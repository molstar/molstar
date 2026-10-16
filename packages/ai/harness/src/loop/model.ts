/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type AssistantMessage,
  digest,
  type Json,
  type ModelConfiguration,
  ModelError,
  type ModelFinish,
  type ProjectedRequest,
  type ProjectInput,
  type ProviderUsage,
  UnrepresentableInput,
} from '@molstar/ai-agent';
import type { AttemptFailure } from '../contracts/records.js';
import { available } from '../tools/plan.js';
import { admit, blob, charge, compact, decide, messagesOf, modelUsage, openCall, type Run, sizeFrom } from './run.js';
import { ZERO } from './state.js';

const failure = (
  code: AttemptFailure['code'],
  message: string,
  transient = false,
  effect: AttemptFailure['effect'] = 'none',
  retryAfterMs?: number,
): AttemptFailure => ({
  code,
  message,
  transient,
  effect,
  ...(retryAfterMs !== undefined ? { retryAfterMs } : {}),
});

/**
 * The model step (harness.md §6): project, admit, invoke, and the decisions on the way. Returns the
 * accepted message; the caller plans its calls, or appends it at a completing boundary.
 */
export async function modelStep(run: Run): Promise<AssistantMessage> {
  const { clock, ids, observations } = run.core;
  const signal = run.abort.signal;
  let configuration: ModelConfiguration = run.model;
  for (;;) {
    let request: ProjectedRequest | undefined;
    let refused: AttemptFailure | undefined;
    let model: Awaited<ReturnType<typeof run.core.model>> | undefined;
    let input: ProjectInput | undefined;
    try {
      model = await run.core.model(configuration);
      const capabilities = model.capabilities;
      input = { messages: messagesOf(run, run.context.entries), tools: available(run), configuration, capabilities };
      request = await model.adapter.project(input);
      const size = sizeFrom(request.inputTokens, capabilities, configuration);
      if (size.inputTokens + size.reservedOutput > size.window)
        refused = failure(
          'capacity',
          `${size.inputTokens} input tokens and ${size.reservedOutput} reserved exceed the window of ${size.window}`,
        );
    } catch (error) {
      if (signal.aborted) throw signal.reason;
      refused = failure(error instanceof UnrepresentableInput ? 'unrepresentable' : 'provider', String(error));
    }
    if (refused !== undefined || request === undefined || model === undefined || input === undefined) {
      const answer = await decide(run, 'failure', { failure: refused ?? failure('provider', 'no request') });
      if ('retry' in answer) await clock.sleep(answer.retry.afterMs ?? 0, signal);
      else if ('replace' in answer) configuration = answer.replace;
      else if ('compact' in answer) await compact(run, answer.compact);
      continue;
    }
    const adapter = model.adapter;
    const evidence: Json = {
      head: run.context.from ?? null,
      entries: digest(run.context.entries),
      configuration: configuration as unknown as Json,
      tools: input.tools.map((tool) => tool.tool.id),
      adapterRevision: request.adapterRevision,
      digest: request.digest,
      inputTokens: request.inputTokens,
      ...(request.segments !== undefined ? { segments: [...request.segments] } : {}),
    };
    const call = openCall(run, 'model', evidence);
    for (let attempt = 1; ; attempt += 1) {
      let fault: AttemptFailure | undefined;
      if (attempt > 1) {
        const again = await adapter.project(input).catch(() => undefined);
        if (again?.digest !== request.digest || again.adapterRevision !== request.adapterRevision)
          fault = failure('integrity', 'a retry projected another request');
      }
      const first = attempt === 1;
      const invoked = fault === undefined;
      if (invoked) {
        const need = { ...ZERO, modelCalls: first ? 1 : 0, uncachedInput: request.inputTokens };
        const cost = run.core.price?.(configuration, need);
        await admit(run, {
          modelCalls: need.modelCalls,
          uncachedInput: need.uncachedInput,
          endBy: 0,
          ...(cost !== undefined ? { cost } : {}),
        });
      }
      const attemptId = `${call.id}:${ids.next('attempt')}`;
      const scope = { path: run.path, callId: call.id, attemptId };
      const startedAt = clock.now();
      let result: { message: AssistantMessage; finish: ModelFinish } | undefined;
      let usage: ProviderUsage | undefined;
      if (invoked) {
        observations.next({ kind: 'attempt.started', ...scope });
        try {
          for await (const frame of adapter.invoke({ request, signal })) {
            if (result !== undefined) {
              fault = failure('stream-protocol', 'a frame after the result', false, 'unknown');
              break;
            }
            if (frame.kind === 'text-delta')
              observations.next({ kind: 'model.delta', ...scope, delta: { text: frame.text } });
            else if (frame.kind === 'reasoning-delta')
              observations.next({ kind: 'model.delta', ...scope, delta: { reasoning: frame.text } });
            else if (frame.kind === 'usage') usage = frame.usage;
            else {
              result = { message: frame.message, finish: frame.finish };
              usage = frame.usage;
            }
          }
          if (result === undefined) fault ??= failure('stream-protocol', 'no result', false, 'unknown');
        } catch (error) {
          if (!signal.aborted)
            fault =
              error instanceof ModelError
                ? failure('provider', error.message, error.transient, error.effect, error.retryAfterMs)
                : failure('provider', String(error), false, 'unknown');
        }
      }
      const output = result === undefined ? 0 : Math.ceil(JSON.stringify(result.message).length / 4);
      let used = invoked ? modelUsage(usage, { input: request.inputTokens, output }, first) : ZERO;
      const cost = run.core.price?.(configuration, used);
      if (cost !== undefined) used = { ...used, cost };
      const outcome = signal.aborted ? 'aborted' : fault !== undefined ? 'failed' : 'succeeded';
      call.attempts.push({
        id: attemptId,
        startedAt,
        endedAt: clock.now(),
        outcome,
        ...(fault !== undefined && !signal.aborted ? { failure: fault } : {}),
        usage: used,
      });
      charge(run, used);
      observations.next({ kind: 'attempt.ended', ...scope, outcome });
      if (signal.aborted) throw signal.reason;
      if (fault !== undefined) {
        const answer = await decide(run, 'failure', { attempt, failure: fault });
        if ('retry' in answer) {
          await clock.sleep(answer.retry.afterMs ?? 0, signal);
          continue;
        }
        if ('replace' in answer) configuration = answer.replace;
        break;
      }
      if (result === undefined) break;
      call.result = blob(run, result.message as unknown as Json);
      if (result.finish === 'output-limit' || result.finish === 'filtered') {
        const answer = await decide(run, 'model-result', { finish: result.finish, message: result.message });
        if (answer !== 'accept') {
          configuration = answer.replace;
          break;
        }
      }
      return result.message;
    }
  }
}
