/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message, ModelAdapter, ModelCapabilities } from '@molstar/ai-agent';
import { canonicalJson } from '@molstar/ai-agent';
import {
  type Harness,
  HarnessError,
  type HarnessOptions,
  type HarnessReader,
  type TurnHandle,
  type TurnOutcome,
} from '../contracts/harness.js';
import type { HarnessObservable, Observation } from '../contracts/observe.js';
import {
  type ActiveTurn,
  type CallRecord,
  type Compaction,
  type ControlEvent,
  type HarnessEvent,
  type Head,
  type HeadId,
  MAIN,
  type PendingStep,
  type ScopeId,
  type StepEvent,
  type StepId,
  type TurnConfiguration,
  type TurnEnd,
  type TurnId,
  type TurnSummary,
} from '../contracts/records.js';
import type { Clock } from '../contracts/runner.js';
import { nest } from '../nested/index.js';
import { Subject } from '../observe/index.js';
import { Journal } from '../tools/journal.js';
import { bind, planCalls, stored, toolCalls } from '../tools/plan.js';
import { toolsStep } from '../tools/step.js';
import { modelStep } from './model.js';
import { append, blob, boundary, type Core, compact, decide, type Run, sizeOf } from './run.js';
import {
  addUsage,
  Collector,
  Exit,
  firstExit,
  Halt,
  headIdOf,
  type Staged,
  StepState,
  Store,
  usageOf,
  ZERO,
} from './state.js';

const realClock: Clock = {
  now: () => Date.now(),
  sleep: (ms, signal) =>
    new Promise<void>((resolve, reject) => {
      if (signal?.aborted) return reject(signal.reason);
      const timer = setTimeout(resolve, ms);
      signal?.addEventListener(
        'abort',
        () => {
          clearTimeout(timer);
          reject(signal.reason);
        },
        { once: true },
      );
    }),
};
let counter = 0;
const randomIds = {
  next: (kind: string) =>
    `${kind}-${Date.now().toString(36)}-${(++counter).toString(36)}-${Math.random().toString(36).slice(2, 8)}`,
};

/** The empty head (a fresh conversation): its id is the same in every harness. */
export const emptyHead = (): Head => ({ id: headIdOf([]), entries: [] });

/** The turn the loop is on: its id, how it is to stop, and the step at work. */
interface Active {
  readonly turnId: TurnId;
  /** Stop was asked (harness.md §5), with its reason. */
  stop?: { reason?: Json };
  /** The step at work (the begin step included). */
  run?: Run;
}

/**
 * The harness (harness.md §4–§11): an in-memory loop that runs one turn at a time. A turn is a begin step, then
 * rounds of a model step and a tools step; each step's writes (entries, heads, the event, blobs) apply at its end, in
 * one synchronous assignment, so a stop or a failure mid-step leaves the journal as it was.
 */
export function createHarness(options: HarnessOptions): Harness {
  const clock = options.clock ?? realClock;
  const ids = options.ids ?? randomIds;
  const id = options.id ?? `harness-${randomIds.next('h')}`;
  const store = new Store();
  const observations = new Subject<Observation>();
  const seqs = new Subject<number>();
  const models = new Map<string, Promise<{ adapter: ModelAdapter; capabilities: ModelCapabilities }>>();
  let active: Active | undefined;

  const core: Core = {
    harnessId: id,
    models: options.models,
    tools: options.tools,
    clock,
    ids,
    store,
    observations,
    maxDepth: options.limits?.depth ?? 4,
    ...(options.price !== undefined ? { price: options.price } : {}),
    model(configuration) {
      const key = canonicalJson(configuration as unknown as Json);
      let found = models.get(key);
      if (found === undefined) {
        found = options.models.resolve(configuration.target).then(async (adapter) => ({
          adapter,
          capabilities: await adapter.capabilitiesFor(configuration),
        }));
        models.set(key, found);
        found.catch(() => models.delete(key));
      }
      return found;
    },
    nest,
    driveJournal: (run) => new Journal(run, `drive:${run.path.join('/')}`),
  };

  // --- Writes: one event at a time, observed as it is written ---

  function write<E extends HarnessEvent>(event: Omit<E, 'seq' | 'at'>): E {
    const full = { ...event, seq: store.seq + 1, at: clock.now() } as E;
    store.events.push(full);
    observations.next({ kind: 'event', event: full });
    seqs.next(full.seq);
    return full;
  }

  const summaryOf = (turn: ActiveTurn, status: TurnEnd['status'], lastStep?: StepId): TurnSummary => ({
    id: turn.id,
    status,
    firstStep: store.firstStep ?? turn.pending.id,
    ...(lastStep !== undefined ? { lastStep } : {}),
    usage: turn.usage,
    endedAt: clock.now(),
  });

  const frozen = (calls: Run['collector']['calls']): CallRecord[] =>
    calls.map((call) => ({ ...call, attempts: [...call.attempts] }));

  /** The step's staged pointers become heads; the messages they name and the step's blobs are stored. */
  function applyMoves(collector: Collector): StepEvent['moves'] {
    const moves: { pointer: ScopeId; from?: HeadId; to: HeadId }[] = [];
    for (const blob of collector.blobs) store.blobs.set(...blob);
    for (const [pointer, staged] of collector.pointers) {
      if (staged.appended !== undefined && staged.appended.length === 0) continue;
      for (const entry of staged.appended ?? staged.entries) {
        const message = collector.messages.get(entry);
        if (message !== undefined && !store.entries.has(entry)) store.entries.set(entry, message);
      }
      const to = store.head(staged.entries);
      if (pointer === MAIN) store.main = to;
      else store.pointers.set(pointer, to);
      moves.push({ pointer, ...(staged.from !== undefined ? { from: staged.from } : {}), to });
    }
    return moves;
  }

  const compactionOf = (collector: Collector, moves: StepEvent['moves']): Compaction[] => {
    const moved = new Set(moves.map((move) => move.pointer));
    return [...collector.compaction.values()].filter((c) => moved.has(c.pointer));
  };

  /** Ends the active turn: its end, and `main` remembered if it completed. */
  function ended(turnId: TurnId, end: TurnEnd): void {
    store.ends.set(turnId, end);
    if (end.status === 'completed') store.completedMain = store.main;
    store.turn = undefined;
    store.firstStep = undefined;
  }

  /** How a step ended (harness.md §6): its writes, its event, and the turn after it. */
  function commit(
    run: Run,
    step: { ending: StepEvent['ending']; next?: PendingStep; end?: TurnEnd; park?: ActiveTurn['park'] },
  ): void {
    const turn = store.turn;
    if (run.state.halted || turn === undefined || turn.id !== run.turnId) return;
    // Stop wins over a settlement in flight (harness.md §5).
    const exit = firstExit(run.state.exits);
    if (step.ending === 'settled' && exit !== undefined) throw exit;
    const calls = frozen(run.collector.calls);
    const moves = step.ending === 'settled' ? applyMoves(run.collector) : [];
    if (step.ending !== 'settled') for (const blob of run.collector.blobs) store.blobs.set(...blob);
    const compaction = compactionOf(run.collector, moves);
    const usage = usageOf(calls);
    const updated: ActiveTurn = {
      ...turn,
      usage: addUsage(turn.usage, usage.all),
      budget: { limits: run.budget.limits, consumed: addUsage(turn.budget.consumed, usage.charged) },
    };
    const summary = step.end !== undefined ? summaryOf(updated, step.end.status, run.step) : undefined;
    const { park: _, ...rest } = updated;
    store.turn =
      step.end !== undefined
        ? undefined
        : { ...rest, pending: step.next ?? turn.pending, ...(step.park !== undefined ? { park: step.park } : {}) };
    if (step.end !== undefined) ended(turn.id, step.end);
    write<StepEvent>({
      kind: 'step',
      turnId: turn.id,
      step: run.step,
      stepKind: turn.pending.kind,
      ending: step.ending,
      calls,
      moves,
      ...(compaction.length > 0 ? { compaction } : {}),
      ...(step.next !== undefined ? { next: step.next } : {}),
      ...(step.end !== undefined ? { end: step.end } : {}),
      ...(summary !== undefined ? { summary } : {}),
      ...(run.budget.limits !== run.state.limits ? { limits: run.budget.limits } : {}),
    });
  }

  /** A step that ends early (harness.md §6): parked, or abandoned with the turn's end. */
  function exitStep(run: Run, exit: Exit): void {
    if (exit.kind === 'park')
      return commit(run, {
        ending: 'parked',
        park: { path: exit.path, ...(exit.note !== undefined ? { note: exit.note } : {}) },
      });
    const stopped = exit.kind === 'stop';
    const reason = stopped ? active?.stop?.reason : exit.reason;
    commit(run, {
      ending: 'abandoned',
      end: { status: stopped ? 'stopped' : 'failed', ...(reason !== undefined ? { reason } : {}) },
    });
  }

  function runOf(
    turn: Pick<ActiveTurn, 'id' | 'configuration' | 'budget'>,
    step: StepId,
    collector: Collector,
    context: Staged,
  ): Run {
    const state = new StepState();
    const lastCompaction = store.lastCompaction(MAIN);
    return {
      core,
      turnId: turn.id,
      configuration: turn.configuration,
      step,
      depth: 0,
      path: [],
      model: turn.configuration.model,
      tools: turn.configuration.tools,
      budget: { limits: { ...turn.budget.limits }, consumed: turn.budget.consumed },
      drives: [options.drive],
      context,
      collector,
      state,
      abort: state.controller,
      seen: new Map(),
      ...(lastCompaction !== undefined ? { lastCompaction } : {}),
    };
  }

  /** The step at work, for stop and reset; a stop asked between steps reaches the next one at once. */
  function open(run: Run): Run {
    if (active !== undefined) {
      active.run = run;
      if (active.stop !== undefined) run.state.exit(new Exit('stop'));
    }
    run.state.limits = run.budget.limits;
    return run;
  }

  /** What went wrong in a step that is not an exit: the turn fails with it rather than the loop. */
  const exitOf = (run: Run, error: unknown): Exit =>
    firstExit(run.state.exits) ??
    (error instanceof Exit ? error : new Exit('fail', { reason: { reason: 'error', message: String(error) } }));

  // --- Steps ---

  /** A round's end (harness.md §6): the boundary, then the commit naming the next step or the turn's end. */
  async function round(run: Run, completing: boolean): Promise<void> {
    await boundary(run, completing);
    commit(run, { ending: 'settled', ...(completing ? { end: { status: 'completed' } } : { next: model() }) });
  }

  const model = (): PendingStep => ({ id: ids.next('step'), kind: 'model' });

  async function modelTop(run: Run): Promise<void> {
    const message = await modelStep(run);
    if (run.abort.signal.aborted) throw run.abort.signal.reason;
    if (toolCalls(message).length > 0) {
      const next = ids.next('step');
      const plan = planCalls(run, message, next);
      const ref = blob(run, message as unknown as Json);
      return commit(run, { ending: 'settled', next: { id: next, kind: 'tools', message: ref, plan: stored(plan) } });
    }
    append(run, run.context, [message]);
    await round(run, true);
  }

  async function toolsTop(run: Run, pending: Extract<PendingStep, { kind: 'tools' }>): Promise<void> {
    const message = store.blobs.get(pending.message) as unknown as Parameters<typeof toolsStep>[1];
    const { round: messages, ids: hashed } = await toolsStep(run, message, bind(run, pending.plan));
    append(run, run.context, messages, hashed);
    await round(run, false);
  }

  /** `endBy` passing during a step opens `budget`; without an extend the step fails (harness.md §11). */
  function endBy(run: Run): AbortController {
    const timer = new AbortController();
    (async () => {
      for (;;) {
        const at = run.budget.limits.endBy;
        if (at === undefined) return;
        if (clock.now() < at) {
          await clock.sleep(at - clock.now(), timer.signal);
          continue;
        }
        const answer = await decide(
          run,
          'budget',
          { dimension: 'endBy', budget: run.budget, need: 0 },
          { reason: 'time' },
        );
        run.budget.limits = { ...run.budget.limits, ...answer.extend };
        if ((run.budget.limits.endBy ?? Number.POSITIVE_INFINITY) <= clock.now())
          run.state.exit(new Exit('fail', { reason: { reason: 'time' } }));
      }
    })().catch(() => undefined);
    return timer;
  }

  /** Runs the pending step and writes how it ended. */
  async function step(turn: ActiveTurn): Promise<void> {
    const collector = new Collector();
    const run = open(
      runOf(turn, turn.pending.id, collector, collector.stage(MAIN, store.main, store.heads.get(store.main) ?? [])),
    );
    const timer = endBy(run);
    try {
      if (run.abort.signal.aborted) throw run.abort.signal.reason;
      if (turn.pending.kind === 'model') await modelTop(run);
      else await toolsTop(run, turn.pending);
    } catch (error) {
      if (run.state.halted || error instanceof Halt) return;
      exitStep(run, exitOf(run, error));
    } finally {
      timer.abort();
    }
  }

  /**
   * The begin step (harness.md §5): `turn-start` (where a drive compacts first), then the turn's messages appended to
   * `main`. A stop or a fail here ends the turn before any model call.
   */
  async function begin(turnId: TurnId, configuration: TurnConfiguration, messages: readonly Message[]) {
    const collector = new Collector();
    const step = ids.next('step');
    const turn = { id: turnId, configuration, budget: { limits: configuration.limits, consumed: ZERO } };
    const run = open(
      runOf(turn, step, collector, collector.stage(MAIN, store.main, store.heads.get(store.main) ?? [])),
    );
    /** The begin step's event: its calls are a drive's nested runs (a compaction's model call). */
    const record = (fields: Pick<StepEvent, 'ending' | 'moves'> & Partial<StepEvent>) => {
      for (const blob of collector.blobs) store.blobs.set(...blob);
      const calls = frozen(collector.calls);
      write<StepEvent>({ kind: 'step', turnId, step, stepKind: 'begin', calls, begun: { configuration }, ...fields });
    };
    try {
      if (run.abort.signal.aborted) throw run.abort.signal.reason;
      const lastCompaction = store.lastCompaction(MAIN);
      const answer = await decide(run, 'turn-start', {
        size: await sizeOf(run),
        ...(lastCompaction !== undefined ? { lastCompaction } : {}),
      });
      if (typeof answer === 'object') await compact(run, answer.compact);
      append(run, run.context, messages);
      if (run.state.halted) return;
      const exit = firstExit(run.state.exits);
      if (exit !== undefined) throw exit;
      const next = model();
      const usage = usageOf(frozen(collector.calls)).all;
      const moves = applyMoves(collector);
      const compaction = compactionOf(collector, moves);
      store.turn = { id: turnId, configuration, budget: turn.budget, pending: next, usage };
      store.firstStep = step;
      record({ ending: 'settled', moves, next, ...(compaction.length > 0 ? { compaction } : {}) });
    } catch (error) {
      if (run.state.halted || error instanceof Halt) return;
      const exit = exitOf(run, error);
      const stopped = exit.kind === 'stop';
      const reason = stopped ? active?.stop?.reason : exit.reason;
      const end: TurnEnd = { status: stopped ? 'stopped' : 'failed', ...(reason !== undefined ? { reason } : {}) };
      const usage = usageOf(frozen(collector.calls)).all;
      const summary: TurnSummary = {
        id: turnId,
        status: end.status,
        firstStep: step,
        lastStep: step,
        usage,
        endedAt: clock.now(),
      };
      store.ends.set(turnId, end);
      record({ ending: 'abandoned', moves: [], end, summary });
    }
  }

  /** The turn's steps until it ends or parks. */
  async function loop(turnId: TurnId, start?: () => Promise<void>): Promise<TurnOutcome> {
    try {
      await Promise.resolve();
      if (start !== undefined) await start();
      for (;;) {
        const end = store.ends.get(turnId);
        if (end !== undefined) return { turnId, end };
        const turn = store.turn;
        if (turn === undefined || turn.id !== turnId) return { turnId, end: { status: 'reset' } };
        if (turn.park !== undefined) {
          // A stop asked as the step parked ends the turn here.
          if (active?.stop === undefined) return { turnId, parked: turn.park };
          endByControl('stopped', active.stop.reason);
          continue;
        }
        await step(turn);
      }
    } finally {
      if (active?.turnId === turnId) active = undefined;
    }
  }

  /** Ends a turn no step is running for (parked, or reset), by a control event. */
  function endByControl(control: 'stopped' | 'reset', reason?: Json, detail: Json = {}): void {
    const turn = store.turn;
    if (turn === undefined) return;
    const summary = summaryOf(turn, control);
    ended(turn.id, { status: control, ...(reason !== undefined ? { reason } : {}) });
    write<ControlEvent>({
      kind: 'control',
      control,
      turnId: turn.id,
      summary,
      detail: { ...(detail as object), ...(reason !== undefined ? { reason } : {}) },
    });
  }

  // --- Reads ---

  const reader: HarnessReader = {
    id,
    execution: async () => ({
      harnessId: id,
      seq: store.seq,
      main: store.main,
      ...(store.turn !== undefined ? { turn: store.turn } : {}),
    }),
    heads: async (wanted) =>
      wanted.map((head) => {
        const entries = store.heads.get(head);
        if (entries === undefined) throw new HarnessError('invalid', `no head ${head}`);
        return { id: head, entries };
      }),
    pointers: async (scopes) =>
      scopes.map((scope) => {
        const head = scope === MAIN ? store.main : store.pointers.get(scope);
        return head === undefined ? undefined : { scope, head };
      }),
    entries: async (wanted) =>
      wanted.map((entry) => {
        const message = store.entries.get(entry);
        if (message === undefined) throw new HarnessError('invalid', `no entry ${entry}`);
        return message;
      }),
    events: async ({ after = 0, limit = 100 }) => store.events.slice(after, after + limit),
    steps: async (wanted) => store.steps(wanted),
    turns: async ({ before, limit = 20 }) => {
      const out: TurnSummary[] = [];
      const turn = store.turn;
      if (before === undefined && turn !== undefined)
        out.push({ id: turn.id, status: 'active', firstStep: store.firstStep ?? turn.pending.id, usage: turn.usage });
      for (
        let i = Math.min(before ?? Number.POSITIVE_INFINITY, store.seq + 1) - 2;
        i >= 0 && out.length < limit;
        i -= 1
      ) {
        const summary = store.events[i]?.summary;
        if (summary !== undefined) out.push(summary);
      }
      return out;
    },
    blobs: async (refs) =>
      refs.map((ref) => {
        const body = store.blobs.get(ref);
        if (body === undefined) throw new HarnessError('invalid', `no blob ${ref}`);
        return body;
      }),
    watch: (): HarnessObservable<number> => ({
      subscribe(onNext, onError, onComplete) {
        const unsubscribe = seqs.subscribe(onNext, onError, onComplete);
        onNext(store.seq);
        return unsubscribe;
      },
    }),
  };

  return {
    ...reader,
    turn({ messages, configuration }): TurnHandle {
      if (active !== undefined || store.turn !== undefined) throw new HarnessError('invalid', 'a turn is active');
      const configured = configuration ?? options.configure?.(messages);
      if (configured === undefined) throw new HarnessError('invalid', 'a turn needs a configuration');
      const turnId = ids.next('turn');
      active = { turnId };
      return { turnId, done: loop(turnId, () => begin(turnId, configured, messages)) };
    },
    stop(reason) {
      if (active !== undefined) {
        active.stop ??= reason !== undefined ? { reason } : {};
        active.run?.state.exit(new Exit('stop'));
      } else if (store.turn !== undefined) endByControl('stopped', reason);
    },
    resume(): TurnHandle {
      const turn = store.turn;
      if (active !== undefined || turn?.park === undefined) throw new HarnessError('invalid', 'no turn is parked');
      const { park: _, ...rest } = turn;
      store.turn = rest;
      write<ControlEvent>({ kind: 'control', control: 'resumed', turnId: turn.id, detail: {} });
      active = { turnId: turn.id };
      return { turnId: turn.id, done: loop(turn.id) };
    },
    reset(heads, reason) {
      const targets: Record<ScopeId, HeadId> = { ...heads };
      if (heads === undefined && store.completedMain !== undefined) targets[MAIN] = store.completedMain;
      for (const head of Object.values(targets))
        if (!store.heads.has(head)) throw new HarnessError('invalid', `no head ${head}`);
      const run = active?.run;
      if (run !== undefined) {
        run.state.halted = true;
        if (!run.abort.signal.aborted) run.abort.abort(new Halt());
      }
      for (const [scope, head] of Object.entries(targets))
        if (scope === MAIN) store.main = head;
        else store.pointers.set(scope, head);
      const detail = { pointers: targets };
      if (store.turn !== undefined) endByControl('reset', reason, detail);
      else
        write<ControlEvent>({
          kind: 'control',
          control: 'reset',
          ...(active !== undefined ? { turnId: active.turnId } : {}),
          detail: { ...detail, ...(reason !== undefined ? { reason } : {}) },
        });
      if (active !== undefined && !store.ends.has(active.turnId)) store.ends.set(active.turnId, { status: 'reset' });
    },
    observe: () => observations,
  };
}
