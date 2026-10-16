/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  digest,
  hashMessage,
  type Json,
  type Message,
  type ModelAdapter,
  type ModelBindings,
  type ModelCapabilities,
  type ModelConfiguration,
  type ProviderUsage,
  type ToolDefinition,
} from '../../messages/index.js';
import type {
  Answer,
  Decision,
  DecisionKind,
  Decisions,
  Drive,
  DriveContext,
  LastCompaction,
} from '../contracts/drive.js';
import type { HeadEdit } from '../contracts/harness.js';
import type { Observation } from '../contracts/observe.js';
import {
  type BlobRef,
  type CallId,
  type ContextSize,
  type Dimension,
  type EntryId,
  type HarnessId,
  type HeadId,
  type Limits,
  MAIN,
  type ScopeId,
  type StepId,
  type ToolSet,
  type TurnConfiguration,
  type TurnId,
  type Usage,
} from '../contracts/records.js';
import type { Clock, IdSource, RunResult, RunSpec, ToolBindings } from '../contracts/runner.js';
import { legal } from '../drive/index.js';
import type { Subject } from '../observe/index.js';
import { projectable } from '../tools/plan.js';
import {
  addUsage,
  type Collector,
  Exit,
  type OpenCall,
  type Staged,
  type StepState,
  type Store,
  ZERO,
} from './state.js';

/** A nested run's way into its parent call: the call, its key, its signal, and the attempt's journal. */
export interface Via {
  readonly callId: CallId;
  readonly key: string;
  readonly signal: AbortSignal;
  readonly journal: Writer;
}

/** The pointer side of a journal: a nested run holds its pointer through it (harness.md §8, §9). */
export interface Writer {
  open(scope: ScopeId): Staged;
}

/** What every run of one harness shares. */
export interface Core {
  readonly harnessId: HarnessId;
  readonly models: ModelBindings;
  readonly tools: ToolBindings;
  readonly clock: Clock;
  readonly ids: IdSource;
  readonly store: Store;
  readonly observations: Subject<Observation>;
  readonly maxDepth: number;
  price?(configuration: ModelConfiguration, usage: Usage): number;
  model(configuration: ModelConfiguration): Promise<{ adapter: ModelAdapter; capabilities: ModelCapabilities }>;
  nest(parent: Run, spec: RunSpec, via: Via): Promise<RunResult>;
  /** A journal for a drive's nested runs (never dropped). */
  driveJournal(run: Run): Writer;
}

/** A top-level turn's step, or a nested run inside it (harness.md §9): one engine for both. */
export interface Run {
  readonly core: Core;
  readonly turnId?: TurnId;
  /** The top-level step this run belongs to; call ids are prefixed with it. */
  readonly step: StepId;
  readonly depth: number;
  readonly path: readonly CallId[];
  /** Nested runs' calls carry it: the tool call (or a drive's step) that started them. */
  readonly parentCallId?: CallId;
  readonly parentKey?: string;
  readonly model: ModelConfiguration;
  readonly tools: ToolSet;
  /** A nested run's own definitions, loaded for it without being in its head. */
  readonly definitions?: readonly ToolDefinition[];
  readonly budget: { limits: Limits; consumed: Usage };
  readonly drives: readonly Drive[];
  /** The transient head. */
  readonly context: Staged;
  readonly collector: Collector;
  readonly state: StepState;
  /** This run's own abort (a nested run's fail aborts its own work only). */
  readonly abort: AbortController;
  readonly configuration?: TurnConfiguration;
  /** Counts identical calls in a nested run, for content-derived keys (harness.md §8). */
  readonly seen: Map<string, number>;
  readonly lastCompaction?: LastCompaction;
}

/** A controller that aborts with `parent`'s reason. */
export function linked(parent: AbortSignal): AbortController {
  const controller = new AbortController();
  if (parent.aborted) controller.abort(parent.reason);
  else parent.addEventListener('abort', () => controller.abort(parent.reason), { once: true });
  return controller;
}

/** Ends work with an exit: a nested run's `fail` ends that run; any other exit ends the step (harness.md §6, §9). */
export function raise(run: Run, exit: Exit): never {
  if (exit.kind === 'fail' && run.depth > 0) {
    if (!run.abort.signal.aborted) run.abort.abort(exit);
    throw exit;
  }
  throw run.state.exit(exit);
}

/** Waits for `work` unless `signal` aborts first. */
export function unlessAborted<T>(signal: AbortSignal, work: Promise<T>): Promise<T> {
  if (signal.aborted) return Promise.reject(signal.reason);
  return new Promise<T>((resolve, reject) => {
    const abort = () => reject(signal.reason);
    signal.addEventListener('abort', abort, { once: true });
    work.then(resolve, reject).finally(() => signal.removeEventListener('abort', abort));
  });
}

const idsOf = (payload: unknown) => {
  const p = payload as { call?: { callId: string }; calls?: { callId: string }[] };
  const calls = p.call !== undefined ? [p.call.callId] : p.calls?.map((call) => call.callId);
  return calls !== undefined ? { calls } : {};
};

const answerKind = (answer: unknown): string =>
  typeof answer === 'string' ? answer : (Object.keys(answer as object)[0] ?? 'unknown');

/** An entry as this step sees it: one it added, else a stored one. */
export const messageOf = (run: Run, id: EntryId): Message | undefined =>
  run.collector.messages.get(id) ?? run.core.store.entries.get(id);

export function messagesOf(run: Run, ids: readonly EntryId[]): Message[] {
  return ids.map((id) => {
    const message = messageOf(run, id);
    if (message === undefined) throw new Error(`no entry ${id}`);
    return message;
  });
}

/** Hashes messages once each and keeps them with the step's writes (harness.md §2). */
export function add(run: Run, messages: readonly Message[]): EntryId[] {
  return messages.map((message) => {
    const id = hashMessage(message);
    run.collector.messages.set(id, message);
    return id;
  });
}

/** A pointer's transient head in this step, staged once. */
export function openPointer(run: Run, scope: ScopeId): Staged {
  const known = run.collector.pointers.get(scope);
  if (known !== undefined) return known;
  const { store } = run.core;
  const head = store.pointers.get(scope);
  return run.collector.stage(scope, head, head === undefined ? [] : (store.heads.get(head) ?? []));
}

/** A pointer's head as this step sees it, its staged moves included (harness.md §8). */
export function readPointer(run: Run, scope: ScopeId): { head?: HeadId; messages: Message[] } {
  const staged = scope === MAIN ? run.collector.pointers.get(MAIN) : openPointer(run, scope);
  if (staged === undefined) return { messages: [] };
  const unchanged = staged.appended !== undefined && staged.appended.length === 0;
  return {
    ...(unchanged && staged.from !== undefined ? { head: staged.from } : {}),
    messages: messagesOf(run, staged.entries),
  };
}

/** Asks the drive stack (harness.md §7); stop and fail answers become exits (`failAs` replaces a fail's reason). */
export async function decide<K extends DecisionKind>(
  run: Run,
  kind: K,
  payload: Decisions[K][0],
  failAs?: Json,
): Promise<Decisions[K][1]> {
  const decision = { kind, path: run.path, payload } as Decision<K>;
  const ids = idsOf(payload);
  const signal = run.abort.signal;
  run.core.observations.next({ kind: 'decision.opened', decisionKind: kind, path: run.path, ids });
  const context: DriveContext = {
    harnessId: run.core.harnessId,
    ...(run.turnId !== undefined ? { turnId: run.turnId } : {}),
    depth: run.depth,
    budget: { limits: run.budget.limits, consumed: run.budget.consumed },
    head: { id: run.context.from ?? '', entries: run.context.entries },
    ...(run.configuration !== undefined ? { configuration: run.configuration } : {}),
    signal,
    journal: { read: async (scope) => readPointer(run, scope) },
    run: (spec) =>
      run.core.nest(run, spec, {
        callId: run.parentCallId ?? `drive:${run.step}`,
        key: run.parentKey ?? run.step,
        signal,
        journal: run.core.driveJournal(run),
      }),
  };
  let answer: Answer<K> | 'pass' = 'pass';
  for (const drive of run.drives) {
    try {
      answer = await unlessAborted(signal, Promise.resolve(drive.decide(decision, context)));
    } catch (error) {
      if (signal.aborted) throw signal.reason;
      raise(run, new Exit('fail', { reason: { reason: 'drive-error', message: String(error) }, path: run.path }));
    }
    if (answer !== 'pass') break;
  }
  if (answer === 'pass' || !legal(decision as Decision, answer))
    raise(run, new Exit('fail', { reason: { reason: 'invalid-answer', kind }, path: run.path }));
  run.core.observations.next({
    kind: 'decision.answered',
    decisionKind: kind,
    path: run.path,
    answer: answerKind(answer),
    ids,
  });
  if (answer === 'stop') raise(run, new Exit('stop', { path: run.path }));
  if (typeof answer === 'object' && answer !== null && 'fail' in answer)
    raise(run, new Exit('fail', { reason: failAs ?? (answer as { fail: Json }).fail, path: run.path }));
  return answer as Decisions[K][1];
}

/** A model attempt's usage (harness.md §11): `inputTokens` is total input, so uncached = input − cache. */
export function modelUsage(
  usage: ProviderUsage | undefined,
  estimate: { input: number; output: number },
  firstAttempt: boolean,
): Usage {
  const reported = usage !== undefined && (usage.inputTokens !== undefined || usage.outputTokens !== undefined);
  const input = reported ? (usage?.inputTokens ?? 0) : estimate.input;
  const cacheRead = usage?.cacheReadInputTokens ?? 0;
  const cacheWrite = usage?.cacheWriteInputTokens ?? 0;
  return {
    ...ZERO,
    modelCalls: firstAttempt ? 1 : 0,
    uncachedInput: Math.max(0, input - cacheRead - cacheWrite),
    cacheRead,
    cacheWrite,
    output: reported ? (usage?.outputTokens ?? 0) : estimate.output,
    reasoning: usage?.reasoningOutputTokens ?? 0,
    ...(reported ? {} : { estimated: true }),
  };
}

/** Admission before an attempt (harness.md §11): a shortfall opens `budget` until it is extended. */
export async function admit(run: Run, need: Partial<Record<Dimension, number>>): Promise<void> {
  for (;;) {
    const { limits, consumed } = run.budget;
    let short: { dimension: Dimension; need: number } | undefined;
    for (const [dimension, amount] of Object.entries(need) as [Dimension, number][]) {
      const limit = limits[dimension];
      if (limit === undefined) continue;
      const over =
        dimension === 'endBy'
          ? run.core.clock.now() >= limit
          : ((consumed[dimension as keyof Usage] as number | undefined) ?? 0) + amount > limit;
      if (over) {
        short = { dimension, need: amount };
        break;
      }
    }
    if (short === undefined) return;
    const answer = await decide(
      run,
      'budget',
      { dimension: short.dimension, budget: { limits, consumed }, need: short.need },
      short.dimension === 'endBy' ? { reason: 'time' } : undefined,
    );
    run.budget.limits = { ...run.budget.limits, ...answer.extend };
  }
}

export const charge = (run: Run, usage: Usage): void => {
  run.budget.consumed = addUsage(run.budget.consumed, usage);
};

/** A body kept by reference (`digest(body)`), written with the step (harness.md §4). */
export function blob(run: Run, body: Json): BlobRef {
  const ref = digest(body);
  run.collector.blobs.set(ref, body);
  return ref;
}

/** Opens a call record in the step: it joins the step event however the step ends (harness.md §6). */
export function openCall(run: Run, kind: 'model' | 'tool', request: Json, id?: CallId, key?: string): OpenCall {
  const callId = id ?? `${run.step}:${run.core.ids.next('call')}`;
  const call: OpenCall = {
    id: callId,
    key: key ?? callId,
    ...(run.parentCallId !== undefined ? { parentCallId: run.parentCallId } : {}),
    kind,
    request: blob(run, request),
    attempts: [],
  };
  run.collector.calls.push(call);
  return call;
}

/** Appends messages to a context, hashing each once (harness.md §2). */
export function append(run: Run, staged: Staged, messages: readonly Message[], ids?: readonly string[]): void {
  run.collector.append(staged, ids ?? add(run, messages), messages);
}

export function sizeFrom(inputTokens: number, capabilities: ModelCapabilities, model: ModelConfiguration): ContextSize {
  const reservedOutput =
    model.parameters.maxOutputTokens ?? capabilities.defaultMaxOutputTokens ?? capabilities.maxOutputTokens ?? 0;
  const window =
    capabilities.maxInputTokens !== undefined
      ? capabilities.maxInputTokens + reservedOutput
      : capabilities.contextWindowTokens;
  return { inputTokens, reservedOutput, window };
}

export async function sizeOf(run: Run, staged: Staged = run.context): Promise<ContextSize> {
  const inputTokens = run.core.store.tokensOf(staged.entries, (id) => messageOf(run, id));
  const model = await run.core.model(run.model).catch(() => undefined);
  // Without an adapter the window is unknown; the model step reports that failure.
  return model === undefined
    ? { inputTokens, reservedOutput: 0, window: 0 }
    : sizeFrom(inputTokens, model.capabilities, run.model);
}

/** A head edit (compaction's replace-selection); false if it would leave an unprojectable head (harness.md §7). */
export async function applyEdit(run: Run, staged: Staged, edit: HeadEdit): Promise<boolean> {
  const fresh = edit.entries.filter((entry): entry is Message => typeof entry !== 'string');
  const ids = add(run, fresh);
  let next = 0;
  const entries = edit.entries.map((entry) => (typeof entry === 'string' ? entry : (ids[next++] as string)));
  if (entries.some((id) => messageOf(run, id) === undefined) || !projectable(messagesOf(run, entries))) return false;
  run.collector.edit(staged, entries);
  if (staged.scope !== undefined)
    run.collector.compaction.set(staged.scope, { pointer: staged.scope, sizeAfter: await sizeOf(run, staged) });
  return true;
}

/** A drive's edit: one that leaves an unprojectable head is an invalid answer (harness.md §7). */
export async function compact(run: Run, edit: HeadEdit): Promise<void> {
  if (!(await applyEdit(run, run.context, edit)))
    raise(run, new Exit('fail', { reason: { reason: 'invalid-answer', kind: 'edit' }, path: run.path }));
}

/** A round boundary (harness.md §7): a compaction loops back; the answer that remains is returned. */
export async function boundary(run: Run, completing: boolean): Promise<'proceed' | 'end'> {
  for (;;) {
    const answer = await decide(run, 'boundary', {
      size: await sizeOf(run),
      completing,
      ...(run.lastCompaction !== undefined ? { lastCompaction: run.lastCompaction } : {}),
    });
    if (typeof answer === 'object' && 'compact' in answer) await compact(run, answer.compact);
    else return answer;
  }
}
