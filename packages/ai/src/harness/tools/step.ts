/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssistantMessage, Json, Message, ModelVisiblePart, UserMessage } from '../../messages/index.js';
import type { SkipAnswer } from '../contracts/drive.js';
import type { AttemptFailure, CallId, Usage } from '../contracts/records.js';
import type { RunContext, ToolFrame } from '../contracts/runner.js';
import { add, admit, blob, charge, decide, linked, openCall, type Run, raise } from '../loop/run.js';
import { Exit, ZERO } from '../loop/state.js';
import { Journal } from './journal.js';
import { type Planned, projectable, viewOf } from './plan.js';

interface Result {
  readonly status: 'succeeded' | 'failed';
  readonly parts: readonly ModelVisiblePart[];
  readonly contributions?: readonly Message[];
}
type ResultFrame = Extract<ToolFrame, { kind: 'result' }>;

const text = (value: string): ModelVisiblePart[] => [{ kind: 'text', text: value }];
const failed = (parts: readonly ModelVisiblePart[], contributions?: readonly Message[]): Result => ({
  status: 'failed',
  parts,
  ...(contributions !== undefined && contributions.length > 0 ? { contributions } : {}),
});
const textOf = (parts: readonly ModelVisiblePart[]) => parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');

/**
 * The step's skips (harness.md §7): a `skip` without calls (every unstarted call), or one naming calls (only
 * those); and which calls have started.
 */
interface Skip {
  text?: string;
  readonly calls: Map<CallId, string>;
  readonly plan: readonly Planned[];
  readonly started: Set<CallId>;
  readonly results: Map<CallId, Result>;
}
const skipping = (skip: Skip, callId: CallId) => skip.calls.get(callId) ?? skip.text;
/** The step's calls that have neither started nor got a result, but `except`. */
const unstarted = (skip: Skip, except: CallId) =>
  skip.plan
    .filter((c) => c.callId !== except && !skip.started.has(c.callId) && !skip.results.has(c.callId))
    .map(viewOf);
/** A decision's `skip`: its calls, or every unstarted one. */
function skipWith(skip: Skip, answer: SkipAnswer | undefined): void {
  if (answer === undefined) return;
  if (answer.calls === undefined) skip.text ??= answer.text;
  else for (const id of answer.calls) if (!skip.calls.has(id)) skip.calls.set(id, answer.text);
}

/** One attempt of a handler under the frame protocol (harness.md §8). */
async function attemptCall(
  run: Run,
  call: Planned,
  attempt: number,
  attemptId: string,
  signal: AbortSignal,
  journal: Journal,
): Promise<{ result?: ResultFrame; fault?: AttemptFailure; usage: Partial<Usage> }> {
  const handler = call.bound === undefined ? undefined : run.core.tools.resolve(call.bound);
  if (handler === undefined || call.bound === undefined)
    return { fault: { code: 'handler', message: 'no handler', transient: false, effect: 'none' }, usage: {} };
  const context: RunContext = {
    harnessId: run.core.harnessId,
    turnId: run.turnId ?? '',
    callId: call.callId,
    path: run.path,
    key: call.key,
    prior: attempt > 1 ? 'possible' : 'none',
    attempt,
    definition: call.bound,
    arguments: call.arguments,
    signal,
    journal,
    run: (spec) => run.core.nest(run, spec, { callId: call.callId, key: call.key, signal, journal }),
  };
  let result: ResultFrame | undefined;
  let usage: Partial<Usage> = {};
  const protocol = (message: string): AttemptFailure => ({
    code: 'stream-protocol',
    message,
    transient: false,
    effect: 'unknown',
  });
  try {
    for await (const frame of handler.invoke(context)) {
      if (result !== undefined) return { fault: protocol('a frame after the result'), usage };
      if (frame.kind === 'progress') {
        run.core.observations.next({
          kind: 'tool.progress',
          path: run.path,
          callId: call.callId,
          attemptId,
          value: frame.value,
        });
      } else if (frame.kind === 'usage') usage = frame.usage;
      else result = frame;
    }
    await journal.settled();
  } catch (error) {
    return { fault: { code: 'handler', message: String(error), transient: false, effect: 'unknown' }, usage };
  }
  if (result === undefined) return { fault: protocol('no result'), usage };
  if (!projectable(result.contributions ?? []))
    return {
      fault: {
        code: 'invalid-contribution',
        message: 'the contributions are not projectable',
        transient: false,
        effect: result.effect ?? 'unknown',
      },
      usage,
    };
  if (result.status === 'failed')
    return {
      result,
      fault: {
        code: 'handler',
        message: textOf(result.parts),
        transient: result.transient ?? false,
        effect: result.effect ?? 'unknown',
      },
      usage,
    };
  return { result, usage };
}

/** A bound call's attempts until it has its one result (harness.md §6). */
async function runCall(run: Run, call: Planned, skip: Skip): Promise<Result> {
  const { clock, ids, observations } = run.core;
  const record = openCall(
    run,
    'tool',
    { toolId: call.toolId, key: call.key, arguments: call.arguments },
    call.callId,
    call.key,
  );
  for (let attempt = 1; ; attempt += 1) {
    await admit(run, { toolCalls: attempt === 1 ? 1 : 0, endBy: 0 });
    const controller = linked(run.abort.signal);
    const journal = new Journal(run, call.callId);
    const attemptId = `${call.callId}:${ids.next('attempt')}`;
    const scope = { path: run.path, callId: call.callId, attemptId };
    const startedAt = clock.now();
    observations.next({ kind: 'attempt.started', ...scope });
    const { result, fault, usage } = await attemptCall(run, call, attempt, attemptId, controller.signal, journal);
    const aborted = controller.signal.aborted;
    const used: Usage = { ...ZERO, ...usage, toolCalls: attempt === 1 ? 1 : 0 };
    const outcome = aborted ? 'aborted' : fault !== undefined ? 'failed' : 'succeeded';
    record.attempts.push({
      id: attemptId,
      startedAt,
      endedAt: clock.now(),
      outcome,
      ...(fault !== undefined && !aborted ? { failure: fault } : {}),
      usage: used,
    });
    charge(run, used);
    observations.next({ kind: 'attempt.ended', ...scope, outcome });
    if (aborted) {
      journal.drop();
      throw controller.signal.reason;
    }
    let reported: Result | undefined;
    if (fault === undefined && result !== undefined) reported = { ...result };
    else {
      journal.drop();
      const answer = await decide(run, 'failure', {
        call: viewOf(call),
        attempt,
        failure: fault as AttemptFailure,
        unstarted: unstarted(skip, call.callId),
      });
      if ('retry' in answer) {
        await clock.sleep(answer.retry.afterMs ?? 0, run.abort.signal);
        continue;
      }
      if ('report' in answer) {
        skipWith(skip, answer.report.skip);
        reported =
          result !== undefined && fault?.code === 'handler'
            ? failed(result.parts, result.contributions)
            : failed(text(fault?.message ?? 'failed'));
      }
    }
    if (reported === undefined) raise(run, new Exit('fail', { reason: { reason: 'invalid-answer', kind: 'failure' } }));
    const body: Result = {
      status: reported.status,
      parts: reported.parts,
      ...(reported.contributions ? { contributions: reported.contributions } : {}),
    };
    record.result = blob(run, body as unknown as Json);
    return body;
  }
}

/** Runs a stage's calls, up to `limit` at once; an exit waits for every started call to end (harness.md §6). */
async function runStage(run: Run, calls: readonly Planned[], limit: number, results: Map<CallId, Result>, skip: Skip) {
  const queue = [...calls];
  const errors: unknown[] = [];
  const worker = async () => {
    for (let call = queue.shift(); call !== undefined; call = queue.shift()) {
      const skipped = skipping(skip, call.callId);
      if (skipped !== undefined) {
        results.set(call.callId, failed(text(skipped)));
        continue;
      }
      skip.started.add(call.callId);
      try {
        results.set(call.callId, await runCall(run, call, skip));
      } catch (error) {
        errors.push(error);
        queue.length = 0;
      }
    }
  };
  await Promise.all(Array.from({ length: Math.max(1, Math.min(limit, calls.length)) }, worker));
  if (run.abort.signal.aborted) throw run.abort.signal.reason;
  if (errors.length > 0) throw errors[0];
}

/**
 * The tools step (harness.md §6): repair, permission, schedule, run. Returns the round: the assistant
 * message, the results in call order, then the contributions not already in the head.
 */
export async function toolsStep(run: Run, message: AssistantMessage, plan: readonly Planned[]) {
  const results = new Map<CallId, Result>();
  const skip: Skip = { calls: new Map(), plan, started: new Set(), results };
  for (const call of plan) {
    if (call.unbound === undefined) continue;
    const { report } = await decide(run, 'repair', {
      call: viewOf(call),
      reason: call.unbound,
      unstarted: unstarted(skip, call.callId),
    });
    if (!projectable(report.contributions ?? []))
      raise(run, new Exit('fail', { reason: { reason: 'invalid-answer', kind: 'repair' }, path: run.path }));
    results.set(
      call.callId,
      failed(report.parts ?? text(`Cannot call ${call.toolId}: ${call.unbound}`), report.contributions),
    );
    skipWith(skip, report.skip);
  }
  let park: Exit | undefined;
  for (const call of plan) {
    if (results.has(call.callId) || run.tools.preauthorized.includes(call.toolId)) continue;
    const answer = await decide(run, 'permission', { call: viewOf(call) });
    if (typeof answer === 'object' && 'deny' in answer)
      results.set(call.callId, failed(text(`Denied: ${answer.deny}`)));
    else if (typeof answer === 'object' && 'park' in answer)
      park ??= new Exit('park', {
        path: run.path,
        ...(answer.park.note !== undefined ? { note: answer.park.note } : {}),
      });
  }
  if (park !== undefined) raise(run, park);
  const runnable = plan.filter((call) => !results.has(call.callId));
  if (runnable.length > 0) {
    const schedule = await decide(run, 'schedule', { calls: runnable.map(viewOf) });
    const byId = new Map(runnable.map((call) => [call.callId, call]));
    for (const stage of schedule.stages)
      await runStage(
        run,
        stage.flatMap((id) => byId.get(id) ?? []),
        schedule.maxConcurrency ?? stage.length,
        results,
        skip,
      );
  }
  if (run.abort.signal.aborted) throw run.abort.signal.reason;
  const user: UserMessage = {
    role: 'user',
    parts: plan.map((call) => {
      const result = results.get(call.callId) ?? failed(text('Skipped'));
      return {
        kind: 'tool-result',
        requestId: call.requestId,
        tool: { id: call.toolId },
        status: result.status,
        parts: result.parts,
      };
    }),
  };
  const round: Message[] = [message, user];
  const ids = add(run, round);
  let present: Set<string> | undefined;
  for (const call of plan) {
    const contributions = results.get(call.callId)?.contributions ?? [];
    if (contributions.length === 0) continue;
    present ??= new Set(run.context.entries);
    const seen = present;
    const hashed = add(run, contributions);
    for (const [i, contribution] of contributions.entries()) {
      const id = hashed[i] as string;
      if (seen.has(id)) continue;
      seen.add(id);
      round.push(contribution);
      ids.push(id);
    }
  }
  return { round, ids };
}
