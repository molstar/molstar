/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message, ModelAdapter, ModelCapabilities } from '@molstar/ai-agent';
import { canonicalJson } from '@molstar/ai-agent';
import type { ItemView, LastCompaction } from '../contracts/drive.js';
import type { Inspectable } from '../contracts/inspect.js';
import type { Observation } from '../contracts/observe.js';
import { type CommitBatch, HarnessError, type Lease, type Move } from '../contracts/ports.js';
import {
  type ActiveTurn,
  type Execution,
  type Item,
  MAIN,
  type PendingStep,
  type Progress,
  SCHEMA_VERSION,
  type Seq,
  type StepId,
} from '../contracts/records.js';
import type { Runner, RunnerExit, RunnerOptions } from '../contracts/runner.js';
import { ZERO } from '../host/writer.js';
import { nest } from '../nested/index.js';
import { Subject } from '../observe/index.js';
import { Journal } from '../tools/journal.js';
import { bind, planCalls, stored, toolCalls } from '../tools/plan.js';
import { toolsStep } from '../tools/step.js';
import { modelStep } from './model.js';
import { append, blob, boundary, type Core, compact, decide, type Run, sizeOf } from './run.js';
import { Collector, Exit, firstExit, Halt, type local, StepState, WorkingSet } from './state.js';

const realClock = {
  now: () => Date.now(),
  sleep: (ms: number, signal?: AbortSignal) =>
    new Promise<void>((resolve, reject) => {
      if (signal?.aborted) return reject(signal.reason);
      const timer = setTimeout(resolve, ms);
      const abort = () => {
        clearTimeout(timer);
        reject(signal?.reason);
      };
      signal?.addEventListener('abort', abort, { once: true });
    }),
};
let counter = 0;
const randomIds = {
  next: (kind: string) =>
    `${kind}-${Date.now().toString(36)}-${(++counter).toString(36)}-${Math.random().toString(36).slice(2, 8)}`,
};

const view = (item: Item): ItemView => ({
  key: item.key,
  messages: item.messages,
  ...(item.options !== undefined ? { options: item.options } : {}),
  submittedAt: item.submittedAt,
  interrupt: item.interrupt !== undefined,
});

/** The engine (harness.md §4–§11): one runner serves turns for one claim, then releases or hands off. */
export function createRunner(options: RunnerOptions): Runner & Inspectable {
  const { port, drive } = options;
  const clock = options.clock ?? realClock;
  const ids = options.ids ?? randomIds;
  const cache = new WorkingSet(port, options.limits?.cacheBytes ?? 0, read);
  const observations = new Subject<Observation>();
  const models = new Map<string, Promise<{ adapter: ModelAdapter; capabilities: ModelCapabilities }>>();
  const own = new Set<StepId>();
  const handled = new Set<string>();
  const life = new AbortController();
  let exec: Execution;
  let lease: Lease;
  let knownSeq: Seq = 0;
  let latestSeq: Seq = 0;
  let committing = false;
  let refreshing: Promise<void> | undefined;
  let lastSent = 0;
  let claimedAt = 0;
  let progress: Progress = { phase: 'idle', since: 0 };
  let current: Run | undefined;
  /** The interrupting item and its turn, until a commit consumes it (harness.md §10). */
  let interrupting: { item: ItemView; turnId: string } | undefined;
  let mainCompaction: LastCompaction | undefined | null = null;

  const halt = (kind: Halt['kind']) => {
    const reason = new Halt(kind);
    if (!life.signal.aborted) life.abort(reason);
  };

  async function read<T>(work: () => Promise<T>): Promise<T> {
    for (let attempt = 1; ; attempt += 1) {
      try {
        return await work();
      } catch (error) {
        if (error instanceof HarnessError || attempt >= 5) throw error;
        await clock.sleep(250 * attempt, life.signal);
      }
    }
  }
  /** A runner operation: retried with its key on an unknown outcome; it counts as a signal. */
  async function send<T>(work: () => Promise<T>): Promise<T> {
    if (life.signal.aborted) throw life.signal.reason;
    for (let attempt = 1; ; attempt += 1) {
      try {
        const result = await work();
        lastSent = clock.now();
        return result;
      } catch (error) {
        const unknown = !(error instanceof HarnessError) || error.code === 'unknown-outcome';
        if (!unknown || attempt >= 5) throw error;
        await clock.sleep(250 * attempt, life.signal);
      }
    }
  }

  /** A commit or begin: the watch holds its refresh until the runner has applied its own change. */
  async function commitOp<T>(work: () => Promise<T>): Promise<T> {
    committing = true;
    try {
      return await send(work);
    } finally {
      committing = false;
    }
  }

  const core: Core = {
    port,
    harnessId: port.id,
    models: options.models,
    tools: options.tools,
    clock,
    ids,
    cache,
    observations,
    maxDepth: options.limits?.depth ?? 4,
    ...(options.price !== undefined ? { price: options.price } : {}),
    inbox: () => exec.inbox.map(view),
    progress(phase, touch = false) {
      if (phase !== progress.phase || touch)
        progress = { ...(current !== undefined ? { step: current.step } : {}), phase, since: clock.now() };
    },
    read,
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

  // --- The Execution: refreshed when seq moves, examined for Stop, reset and interrupts ---

  function examine(): void {
    if (exec.runner.phase !== 'live' || exec.runner.runnerId !== options.runnerId) {
      halt('lost');
      return;
    }
    if (
      interrupting !== undefined &&
      (exec.turn?.id !== interrupting.turnId || !exec.inbox.some((item) => item.key === interrupting?.item.key))
    )
      interrupting = undefined;
    const run = current;
    if (run === undefined || run.state.closed) return;
    if (exec.turn?.id !== run.turnId) {
      run.state.controller.abort(new Halt('reset'));
      return;
    }
    if (exec.turn?.stopRequested !== undefined) {
      run.state.exit(new Exit('stop'));
      return;
    }
    for (const item of exec.inbox) {
      if (item.interrupt?.turnId !== run.turnId || handled.has(item.key) || run.turnId === undefined) continue;
      const turnId = run.turnId;
      handled.add(item.key);
      const seen = view(item);
      const stepKind = exec.turn?.pending.kind === 'tools' ? 'tools' : 'model';
      const asked = decide(run, 'interrupt', { item: seen, stepKind }).then(
        (answer) => {
          interrupting = { item: seen, turnId };
          if (answer === 'abandon') run.state.exit(new Exit('interrupt', { item: seen }));
          else run.state.skip ??= 'Skipped: interrupted';
        },
        () => undefined,
      );
      run.state.pending.push(asked);
    }
  }

  async function refresh(): Promise<void> {
    if (refreshing !== undefined) return await refreshing;
    refreshing = (async () => {
      try {
        exec = await read(() => port.execution());
        knownSeq = exec.seq;
        latestSeq = Math.max(latestSeq, knownSeq);
        examine();
      } finally {
        refreshing = undefined;
      }
    })();
    return await refreshing;
  }

  function watch(): () => void {
    let unsubscribe: (() => void) | undefined;
    const subscribe = () => {
      if (life.signal.aborted) return;
      unsubscribe = port.watch().subscribe(
        (seq) => {
          if (life.signal.aborted) return;
          latestSeq = Math.max(latestSeq, seq);
          if (!committing && latestSeq > knownSeq) refresh().catch(() => undefined);
        },
        () => {
          clock.sleep(1_000, life.signal).then(subscribe, () => undefined);
        },
      );
    };
    subscribe();
    return () => unsubscribe?.();
  }

  async function heartbeat(): Promise<void> {
    const every = options.limits?.signalEveryMs ?? Math.max(1, Math.floor((lease.until - clock.now()) / 3));
    while (!life.signal.aborted) {
      const wait = lastSent + every - clock.now();
      if (wait > 0) {
        await clock.sleep(wait, life.signal).catch(() => undefined);
        continue;
      }
      try {
        lease = await port.renew(lease, { progress });
      } catch (error) {
        if (error instanceof HarnessError && error.code === 'stale') {
          halt('lost');
          return;
        }
      }
      lastSent = clock.now();
    }
  }

  // --- Commits ---

  async function observeEvent(seq: Seq): Promise<void> {
    const events = await read(() => port.events({ after: seq - 1, limit: 1 })).catch(() => []);
    for (const event of events) observations.next({ kind: 'event', event });
  }

  /** After an acknowledged commit: the Execution updated locally (or re-read), eviction, the event. */
  async function committed(
    seq: Seq,
    change: (turn: ActiveTurn) => ActiveTurn | undefined,
    consumed: readonly string[],
    main?: Move,
  ) {
    if (seq === knownSeq + 1) {
      const turn = exec.turn === undefined ? undefined : change(exec.turn);
      const { turn: _, ...rest } = exec;
      exec = {
        ...rest,
        seq,
        main: main?.head.id ?? exec.main,
        inbox: exec.inbox.filter((item) => !consumed.includes(item.key)),
        ...(turn !== undefined ? { turn } : {}),
      };
      knownSeq = seq;
    } else await refresh();
    if (interrupting !== undefined && (consumed.includes(interrupting.item.key) || exec.turn === undefined))
      interrupting = undefined;
    const head = cache.heads.get(exec.main);
    if (head !== undefined) cache.keep([head]);
    const pendingMessage = exec.turn?.pending.kind === 'tools' ? exec.turn.pending.message : undefined;
    for (const ref of cache.blobs.keys()) if (ref !== pendingMessage) cache.blobs.delete(ref);
    await observeEvent(seq);
    if (latestSeq > knownSeq) await refresh();
  }

  /** Commits how a step ended; `collector` holds the moves (an interrupted step's are its own). */
  async function commit(
    run: Run,
    batch: Omit<CommitBatch, 'key' | 'turnId' | 'step' | 'calls' | 'progress' | 'limits'>,
    collector = run.collector,
  ) {
    if (life.signal.aborted) throw life.signal.reason;
    run.state.closed = true;
    await Promise.all(run.state.pending);
    // An exit answered while the step finished (an interrupt's `abandon`) wins over its settlement.
    const exit = firstExit(run.state.exits);
    if (batch.ending === 'settled' && exit !== undefined) throw exit;
    const blobs = [...run.collector.blobs].map(([ref, body]) => ({ ref, body }));
    if (blobs.length > 0) await send(() => port.putBlobs(lease, blobs));
    const full: CommitBatch = {
      key: `${run.step}/${ids.next('key')}`,
      turnId: run.turnId as string,
      step: run.step,
      calls: run.collector.calls,
      progress,
      ...batch,
      ...(run.budget.limits !== run.state.limits ? { limits: run.budget.limits } : {}),
    };
    const done = await commitOp(() => port.commit(lease, full));
    lease = done.lease;
    for (const move of full.moves) {
      const staged = collector.pointers.get(move.pointer);
      for (const id of staged?.appended ?? staged?.entries ?? []) cache.stored.add(id);
    }
    if (full.next !== undefined) own.add(full.next.id);
    const main = full.moves.find((move) => move.pointer === MAIN);
    const compaction = full.compaction?.find((c) => c.pointer === MAIN);
    if (compaction !== undefined) mainCompaction = { seq: done.seq, sizeAfter: compaction.sizeAfter };
    await committed(
      done.seq,
      (turn) => {
        if (full.end !== undefined) return undefined;
        const { park: _, ...rest } = turn;
        return {
          ...rest,
          pending: full.next ?? turn.pending,
          budget: { limits: full.limits ?? turn.budget.limits, consumed: run.budget.consumed },
          ...(full.park !== undefined ? { park: full.park } : {}),
        };
      },
      full.consumed,
      main,
    );
  }

  /** A settling commit: moves, compactions and consumed items from the step's collector. */
  async function settle(run: Run, consumed: readonly string[], next: PendingStep | undefined) {
    const moves = await run.collector.moves(cache);
    const moved = new Set(moves.map((move) => move.pointer));
    const compaction = [...run.collector.compaction.values()].filter((c) => moved.has(c.pointer));
    await commit(run, {
      ending: 'settled',
      messages: run.collector.referenced(cache),
      moves,
      consumed,
      ...(compaction.length > 0 ? { compaction } : {}),
      ...(next !== undefined ? { next } : { end: { status: 'completed' } }),
    });
  }

  /** The items a boundary steers, the interrupting one first (harness.md §10); undefined if one no longer waits. */
  function steered(run: Run, answer: Awaited<ReturnType<typeof boundary>>): Item[] | undefined {
    const keys = [
      ...(interrupting !== undefined && interrupting.turnId === run.turnId ? [interrupting.item.key] : []),
      ...(typeof answer === 'object' ? answer.steer : []),
    ];
    const items = [...new Set(keys)].map((key) => exec.inbox.find((waiting) => waiting.key === key));
    return items.every((item) => item !== undefined) ? items : undefined;
  }

  /** A round's end at depth 0: boundary, steering, and the commit; asked again if a steered item is gone. */
  async function round(run: Run, completing: boolean): Promise<void> {
    const gone = new Set<string>();
    for (;;) {
      const answer = await boundary(run, completing);
      const before = { entries: [...run.context.entries], appended: run.context.appended && [...run.context.appended] };
      const items = steered(run, answer);
      if (items === undefined) {
        // A steered item is gone (withdrawn meanwhile): ask again; naming it again is an invalid answer.
        const missing = (typeof answer === 'object' ? answer.steer : []).filter(
          (key) => !exec.inbox.some((item) => item.key === key),
        );
        if (missing.some((key) => gone.has(key)))
          throw run.state.exit(new Exit('fail', { reason: { reason: 'invalid-answer', kind: 'boundary' } }));
        for (const key of missing) gone.add(key);
        if (latestSeq > knownSeq) await refresh();
        continue;
      }
      await append(
        run,
        run.context,
        items.flatMap((item) => item.messages),
      );
      const next = !completing || items.length > 0 ? { id: ids.next('step'), kind: 'model' as const } : undefined;
      try {
        return await settle(
          run,
          items.map((item) => item.key),
          next,
        );
      } catch (error) {
        if (!(error instanceof HarnessError && error.code === 'inbox-changed')) throw error;
        run.context.entries = before.entries;
        run.context.appended = before.appended;
        await refresh();
      }
    }
  }

  /** A step that ends early (harness.md §6). */
  async function exitStep(run: Run, exit: Exit): Promise<void> {
    const none = { messages: [] as Message[], moves: [] as Move[], consumed: [] as string[] };
    if (exit.kind === 'park')
      return await commit(run, {
        ending: 'parked',
        ...none,
        park: { path: exit.path, ...(exit.note !== undefined ? { note: exit.note } : {}) },
      });
    if (exit.kind === 'interrupt' && exit.item !== undefined) {
      const next = { id: ids.next('step'), kind: 'model' as const };
      const { record, entries } = await mainHead();
      const collector = new Collector();
      const main = collector.stage(MAIN, record, entries);
      collector.append(main, await cache.add(exit.item.messages), exit.item.messages);
      try {
        const moves = await collector.moves(cache);
        const batch = {
          ending: 'interrupted',
          messages: exit.item.messages,
          moves,
          consumed: [exit.item.key],
          next,
        } as const;
        return await commit(run, batch, collector);
      } catch (error) {
        if (!(error instanceof HarnessError && error.code === 'inbox-changed')) throw error;
        // The interrupting item was withdrawn: the step still ends interrupted, with its attempts.
        return await commit(run, { ending: 'interrupted', ...none, next }, new Collector());
      }
    }
    const stopped = exit.kind === 'stop';
    const reason = stopped ? exec.turn?.stopRequested?.reason : exit.reason;
    await commit(run, {
      ending: 'abandoned',
      ...none,
      end: { status: stopped ? 'stopped' : 'failed', ...(reason !== undefined ? { reason } : {}) },
    });
  }

  async function mainHead() {
    const loaded = await cache.load(exec.main);
    if (mainCompaction === null) mainCompaction = await read(() => port.lastCompaction(MAIN));
    return loaded;
  }

  function runOf(
    turn: ActiveTurn | undefined,
    step: StepId,
    collector: Collector,
    context: ReturnType<typeof local>,
  ): Run {
    const state = new StepState();
    return {
      core,
      ...(turn !== undefined ? { turnId: turn.id, configuration: turn.configuration } : {}),
      step,
      depth: 0,
      path: [],
      model: turn?.configuration.model ?? { target: { provider: '', model: '' }, parameters: {} },
      tools: turn?.configuration.tools ?? { enabled: [], preauthorized: [] },
      budget: { limits: { ...turn?.budget.limits }, consumed: turn?.budget.consumed ?? ZERO },
      drives: [drive],
      context,
      collector,
      state,
      abort: state.controller,
      rerun: turn !== undefined && !own.has(step),
      seen: new Map(),
      ...(mainCompaction ? { lastCompaction: mainCompaction } : {}),
    };
  }

  // --- Steps ---

  /** Runs the pending step and commits how it ended. */
  async function step(turn: ActiveTurn, rerun: boolean): Promise<void> {
    const { record, entries } = await mainHead();
    const collector = new Collector();
    const run = runOf(turn, turn.pending.id, collector, collector.stage(MAIN, record, entries));
    const effective: Run = rerun ? { ...run, rerun: true } : run;
    effective.state.limits = effective.budget.limits;
    const stopLife = () => effective.state.controller.abort(life.signal.reason);
    life.signal.addEventListener('abort', stopLife, { once: true });
    current = effective;
    const timer = endBy(effective);
    try {
      examine();
      if (effective.abort.signal.aborted) throw effective.abort.signal.reason;
      if (turn.pending.kind === 'model') await modelTop(effective);
      else await toolsTop(effective, turn.pending);
    } catch (error) {
      if (error instanceof Halt) throw error;
      if (life.signal.aborted) throw life.signal.reason;
      // Stop wins over any other commit: the refused step is abandoned with its attempts (harness.md §5).
      const stopWins = async (refused: unknown) => {
        if (!(refused instanceof HarnessError && refused.code === 'stale')) return false;
        await refresh();
        return exec.turn?.id === turn.id && exec.turn.stopRequested !== undefined;
      };
      let exit = firstExit(effective.state.exits) ?? (error instanceof Exit ? error : undefined);
      if (exit === undefined) {
        if (!(await stopWins(error))) throw error;
        exit = new Exit('stop');
      }
      for (;;) {
        try {
          return await exitStep(effective, exit);
        } catch (refused) {
          if (exit.kind === 'stop' || !(await stopWins(refused))) throw refused;
          exit = new Exit('stop');
        }
      }
    } finally {
      timer.abort();
      current = undefined;
      life.signal.removeEventListener('abort', stopLife);
    }
  }

  async function modelTop(run: Run): Promise<void> {
    const message = await modelStep(run);
    await Promise.all(run.state.pending);
    if (run.abort.signal.aborted) throw run.abort.signal.reason;
    if (toolCalls(message).length > 0) {
      const next = ids.next('step');
      const plan = await planCalls(run, message, next);
      const ref = await blob(run, message as unknown as Json);
      return await settle(run, [], { id: next, kind: 'tools', message: ref, plan: stored(plan) });
    }
    await append(run, run.context, [message]);
    await round(run, true);
  }

  async function toolsTop(run: Run, pending: Extract<PendingStep, { kind: 'tools' }>) {
    const body = cache.blobs.get(pending.message) ?? (await read(() => port.blobs([pending.message])))[0];
    const message = body as unknown as Parameters<typeof toolsStep>[1];
    const { round: messages, ids: hashed } = await toolsStep(run, message, await bind(run, pending.plan));
    await append(run, run.context, messages, hashed);
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

  /** No active turn: `next-turn`, then the begin step (harness.md §5). Returns whether a turn began. */
  async function nextTurn(): Promise<boolean> {
    for (;;) {
      const { record, entries } = await mainHead();
      const collector = new Collector();
      const run = runOf(undefined, ids.next('step'), collector, collector.stage(MAIN, record, entries));
      current = undefined;
      const stopLife = () => run.state.controller.abort(life.signal.reason);
      life.signal.addEventListener('abort', stopLife, { once: true });
      try {
        const configured = exec.inbox.length > 0 ? await options.configure(exec.inbox) : undefined;
        const sized = configured !== undefined ? { ...run, model: configured.model } : run;
        const size =
          configured !== undefined
            ? await sizeOf(sized)
            : { inputTokens: cache.tokensOf(entries), reservedOutput: 0, window: 0 };
        const answer = await decide(run, 'next-turn', {
          size,
          ...(mainCompaction ? { lastCompaction: mainCompaction } : {}),
        });
        if (answer === 'release') return false;
        const items = answer.begin.keys.map((key) => exec.inbox.find((item) => item.key === key));
        if (items.some((item) => item === undefined) || items.length === 0) {
          // Items that no longer wait: ask again with a fresh inbox, or release if it is already fresh.
          if (latestSeq <= knownSeq) return false;
          await refresh();
          continue;
        }
        const chosen = items as Item[];
        const configuration = await options.configure(chosen);
        const turnRun = { ...run, model: configuration.model };
        if (answer.begin.compact !== undefined) await compact(turnRun, answer.begin.compact);
        await append(
          turnRun,
          run.context,
          chosen.flatMap((item) => item.messages),
        );
        const [move] = await collector.moves(cache);
        if (move === undefined) throw new HarnessError('invalid', "a begin step appends its items' messages");
        const turnId = ids.next('turn');
        const next = { id: ids.next('step'), kind: 'model' as const };
        const compaction = collector.compaction.get(MAIN);
        // A drive's nested run (a compaction's model call) is evidence of the begin step.
        const blobs = [...collector.blobs].map(([ref, body]) => ({ ref, body }));
        if (blobs.length > 0) await send(() => port.putBlobs(lease, blobs));
        const done = await commitOp(() =>
          port.begin(lease, {
            key: `${run.step}/begin`,
            turnId,
            step: run.step,
            items: chosen.map((item) => item.key),
            configuration,
            messages: collector.referenced(cache),
            move,
            ...(compaction !== undefined ? { compaction: [compaction] } : {}),
            calls: collector.calls,
            next,
            progress,
          }),
        );
        lease = done.lease;
        for (const id of run.context.appended ?? run.context.entries) cache.stored.add(id);
        own.add(next.id);
        if (compaction !== undefined) mainCompaction = { seq: done.seq, sizeAfter: compaction.sizeAfter };
        exec = {
          ...exec,
          turn: {
            id: turnId,
            schemaVersion: options.schemaVersion ?? SCHEMA_VERSION,
            configuration,
            budget: { limits: configuration.limits, consumed: ZERO },
            pending: next,
            usage: ZERO,
          },
        };
        await committed(
          done.seq,
          (turn) => turn,
          chosen.map((item) => item.key),
          move,
        );
        return true;
      } catch (error) {
        if (error instanceof HarnessError && error.code === 'inbox-changed') {
          await refresh();
          continue;
        }
        if (error instanceof Exit) return false;
        throw error;
      } finally {
        life.signal.removeEventListener('abort', stopLife);
      }
    }
  }

  async function release(): Promise<RunnerExit | undefined> {
    try {
      const released = await commitOp(() => port.release(lease, { seenSeq: knownSeq }));
      life.abort(new Halt('abort'));
      return { reason: 'released', seq: released.seq };
    } catch (error) {
      if (!(error instanceof HarnessError && error.code === 'inbox-changed')) throw error;
      await refresh();
      return undefined;
    }
  }

  async function serve(): Promise<RunnerExit> {
    const retried = new Set<StepId>();
    for (;;) {
      if (life.signal.aborted) throw life.signal.reason;
      const turn = exec.turn;
      if (turn === undefined && (await nextTurn())) continue;
      if (turn === undefined || turn.park !== undefined || turn.blocked !== undefined) {
        const released = await release();
        if (released !== undefined) return released;
        continue;
      }
      const handOff = options.limits?.handOffAfterMs;
      if (handOff !== undefined && clock.now() - claimedAt >= handOff) {
        const released = await release();
        if (released !== undefined) return { ...released, reason: 'handed-off' };
        continue;
      }
      try {
        await step(turn, retried.has(turn.pending.id));
      } catch (error) {
        if (error instanceof Halt && error.kind === 'reset') {
          await refresh();
          continue;
        }
        if (!(error instanceof HarnessError) || !['stale', 'head-changed', 'inbox-changed'].includes(error.code))
          throw error;
        retried.add(turn.pending.id);
        await refresh();
        if (life.signal.aborted) throw life.signal.reason;
        // `stale` with the turn and step unchanged (and no Stop) is the lease: the hold is over (harness.md §4, §16).
        const runner = exec.runner;
        const same =
          exec.turn?.id === turn.id &&
          exec.turn.pending.id === turn.pending.id &&
          exec.turn.stopRequested === undefined;
        if (
          error.code === 'stale' &&
          (runner.phase !== 'live' || runner.runnerId !== options.runnerId || runner.leaseUntil <= clock.now() || same)
        )
          throw new Halt('lost');
      }
    }
  }

  return {
    observe: () => observations,
    retention: () => cache.retention(),
    abort() {
      halt('abort');
    },
    async run(launch) {
      try {
        const claimed = await send(() =>
          port.claim({ launch, runnerId: options.runnerId, schemaVersion: options.schemaVersion ?? SCHEMA_VERSION }),
        );
        lease = claimed.lease;
        exec = claimed.execution;
        knownSeq = latestSeq = exec.seq;
      } catch (error) {
        if (error instanceof HarnessError && error.code === 'stale') return { reason: 'stale-launch' };
        if (error instanceof HarnessError && error.code === 'stale-schema') return { reason: 'stale-schema' };
        throw error;
      }
      claimedAt = clock.now();
      progress = { phase: 'idle', since: claimedAt };
      const unwatch = watch();
      const beating = heartbeat();
      try {
        return await serve();
      } catch (error) {
        const reason = life.signal.aborted ? life.signal.reason : error;
        if (reason instanceof Halt) return { reason: reason.kind === 'abort' ? 'aborted' : 'lost' };
        throw error;
      } finally {
        if (!life.signal.aborted) life.abort(new Halt('abort'));
        unwatch();
        await beating;
      }
    },
  };
}
