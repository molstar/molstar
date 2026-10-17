/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FrameLane, requestFrame } from './frame-state.js';

// Work wanted, done once (rendering.md §2.7): what a change asks to be computed (a tree flattened again, a list
// filtered, a count) is wanted under a key, and a later want of the same key replaces it (the last write wins), so a
// burst of changes costs one run. A queue runs in its lane: `task` (the runtime's microtask flush, before the watches it
// feeds render: `flush()` runs it), `frame` (the animation frame, after the frame states publish, before the renders) or
// `idle` (when the page is idle). Each run first asks the wanted work whether it is `stillWanted`, orders it by
// `priority` (asked again each run) and, with a `budgetMs`, stops when the budget is spent: the rest is carried to the
// next run. Nothing is allocated a want of a known key, nor a run.

export type WorkLane = 'task' | 'frame' | 'idle';

export interface WantOptions {
  /** Higher runs first; a function is asked at each run (a row's distance from the window, now). */
  priority?: number | (() => number);
  /** Asked before it runs: false drops it (a row scrolled away before its turn). */
  stillWanted?: () => boolean;
}

interface Entry {
  run: () => void;
  priority: number | (() => number) | undefined;
  stillWanted: (() => boolean) | undefined;
  /** This run's priority (`priority`, asked once a run). */
  rank: number;
  wanted: boolean;
}

/** Queues with work waiting, by lane (reused, emptied by count). */
const waiting: Record<WorkLane, (WorkQueue | undefined)[]> = { task: [], frame: [], idle: [] };
const waitingCount: Record<WorkLane, number> = { task: 0, frame: 0, idle: 0 };
let taskScheduled = false;
let idleScheduled = false;

/**
 * How the task lane runs: the runtime sets it (dom/scope.ts: its `flush`, which runs the task lane's work first, then
 * the watches). Without the runtime, the work alone.
 */
export const TaskLane: { flush: (() => void) | undefined } = { flush: undefined };

const now: () => number = typeof performance !== 'undefined' ? () => performance.now() : () => Date.now();

const idle: (fn: (deadline?: { timeRemaining(): number }) => void) => void =
  typeof requestIdleCallback === 'function'
    ? (fn) => requestIdleCallback(fn, { timeout: 200 })
    : (fn) => setTimeout(fn, 1);

function enlist(q: WorkQueue) {
  const lane = q.lane;
  waiting[lane][waitingCount[lane]++] = q;
  if (lane === 'task') {
    if (!taskScheduled) {
      taskScheduled = true;
      queueMicrotask(runTaskLaneQueued);
    }
  } else if (lane === 'frame') requestFrame();
  else if (!idleScheduled) {
    idleScheduled = true;
    idle(runIdleLane);
  }
}

function runTaskLaneQueued() {
  taskScheduled = false;
  if (TaskLane.flush) TaskLane.flush();
  else runWork('task');
}

function runIdleLane(deadline?: { timeRemaining(): number }) {
  idleScheduled = false;
  runWork('idle', deadline ? deadline.timeRemaining() : undefined);
}

/** Whether a lane has work waiting (the runtime's flush loops while the task lane does). */
export function hasWork(lane: WorkLane): boolean {
  return waitingCount[lane] > 0;
}

/**
 * Runs the work waiting in `lane` (the runtime calls it: the task lane in `flush`, the frame lane in the frame's tick);
 * `limitMs` caps the whole lane's run (the idle lane's deadline). What a run wants again waits for the lane's next run.
 */
export function runWork(lane: WorkLane, limitMs?: number): void {
  const n = waitingCount[lane];
  if (n === 0) return;
  // This run's queues; a queue wanting work while they run enlists again, for the next.
  const list = waiting[lane];
  waiting[lane] = spare[lane];
  spare[lane] = list;
  waitingCount[lane] = 0;
  const end = limitMs === undefined ? Number.POSITIVE_INFINITY : now() + limitMs;
  for (let i = 0; i < n; i++) {
    const q = list[i]!;
    list[i] = undefined;
    q.runNow(end);
  }
}
const spare: Record<WorkLane, (WorkQueue | undefined)[]> = { task: [], frame: [], idle: [] };

// The frame lane's work runs in the frame's tick, after the frame states publish, before the renders.
FrameLane.work = () => runWork('frame');

/** Idle entries kept past this many are dropped at the end of a run (a queue keyed by rows that come and go). */
const KeepIdle = 512;

/** A queue of work by key (`workQueue`). */
export class WorkQueue {
  /**
   * Every key's entry, wanted or idle: a key wanted again reuses its entry, and the map is not changed (a map's deletes
   * and sets allocate as it rehashes). Idle ones past `KeepIdle` are dropped after a run.
   */
  readonly #entries = new Map<unknown, Entry>();
  #wanted = 0;
  /** A run's order (reused). */
  readonly #order: (Entry | undefined)[] = [];
  #enlisted = false;
  #active = true;
  #disposed = false;
  #running = false;
  #ranked = 0;
  /** Ranks a wanted entry into this run's order (made once: `forEach` with it allocates nothing, `for…of` an iterator). */
  readonly #rank = (e: Entry) => {
    if (!e.wanted) return;
    const order = this.#order;
    const p = e.priority;
    e.rank = p === undefined ? 0 : typeof p === 'number' ? p : p();
    let j = this.#ranked++;
    while (j > 0 && order[j - 1]!.rank < e.rank) {
      order[j] = order[j - 1];
      j--;
    }
    order[j] = e;
  };
  readonly #dropIdle = (e: Entry, key: unknown) => {
    if (!e.wanted) this.#entries.delete(key);
  };

  constructor(
    readonly lane: WorkLane = 'frame',
    /** Each run's budget (ms): the rest waits for the next run. None: all of it. */
    readonly budgetMs?: number,
  ) {}

  /** Work waiting. */
  get size(): number {
    return this.#wanted;
  }

  has(key: unknown): boolean {
    return !!this.#entries.get(key)?.wanted;
  }

  /** Wants `run` done in the lane's next run; a want of the same key before it runs replaces it. */
  want(key: unknown, run: () => void, options?: WantOptions): void {
    if (this.#disposed) return;
    let e = this.#entries.get(key);
    if (!e) {
      e = { run, priority: undefined, stillWanted: undefined, rank: 0, wanted: false };
      this.#entries.set(key, e);
    }
    if (!e.wanted) {
      e.wanted = true;
      this.#wanted++;
    }
    e.run = run;
    e.priority = options?.priority;
    e.stillWanted = options?.stillWanted;
    this.#enlist();
  }

  /** Drops the work wanted under `key`. */
  cancel(key: unknown): void {
    const e = this.#entries.get(key);
    if (e?.wanted) this.#release(e);
  }

  /** Runs all the work waiting now, regardless of the budget (code that reads what it computes, a test). */
  flush(): void {
    this.runNow(Number.POSITIVE_INFINITY, true);
  }

  /** Paused, it runs nothing (hidden UI: `workQueue` in dom/scope.ts follows its scope); resumed, what waited runs. */
  setActive(active: boolean): void {
    if (this.#active === active) return;
    this.#active = active;
    if (active) this.#enlist();
  }

  dispose(): void {
    this.#disposed = true;
    this.#entries.clear();
    this.#wanted = 0;
  }

  /** @internal The lane's run: until `end` (the lane's limit), or the queue's budget. */
  runNow(end: number, all = false): void {
    this.#enlisted = false;
    if (this.#disposed || this.#running || (!this.#active && !all) || this.#wanted === 0) return;
    this.#running = true;
    try {
      const stop = all
        ? Number.POSITIVE_INFINITY
        : this.budgetMs === undefined
          ? end
          : Math.min(end, now() + this.budgetMs);
      // The work as it stands, ordered by priority (asked once each), highest first; in the order wanted otherwise.
      const order = this.#order;
      this.#ranked = 0;
      this.#entries.forEach(this.#rank);
      const n = this.#ranked;
      for (let i = 0; i < n; i++) {
        const e = order[i]!;
        // Cancelled while this run ran: skip. (Wanted again: it runs what was wanted last.)
        if (!e.wanted) continue;
        // At least one a run, however long it takes.
        if (i > 0 && stop !== Number.POSITIVE_INFINITY && now() >= stop) break;
        const run = e.run;
        const still = e.stillWanted;
        this.#release(e);
        if (still && !still()) continue;
        run();
      }
      for (let k = 0; k < n; k++) order[k] = undefined;
      if (this.#entries.size > KeepIdle + this.#wanted) this.#entries.forEach(this.#dropIdle);
    } finally {
      this.#running = false;
    }
    // Carried over (the budget spent) or wanted while running: the lane's next run.
    this.#enlist();
  }

  #enlist() {
    if (this.#enlisted || this.#running || !this.#active || this.#disposed || this.#wanted === 0) return;
    this.#enlisted = true;
    enlist(this);
  }

  #release(e: Entry) {
    e.wanted = false;
    this.#wanted--;
    e.run = noop;
    e.priority = undefined;
    e.stillWanted = undefined;
  }
}

const noop = () => {};

/** A queue of work by key, in `lane` (`WorkQueue`); `dom/scope.ts`'s `workQueue` also pauses and ends it with its scope. */
export function createWorkQueue(options?: { lane?: WorkLane; budgetMs?: number }): WorkQueue {
  return new WorkQueue(options?.lane ?? 'frame', options?.budgetMs);
}

/**
 * The latest call only (rendering.md §2.7): each call aborts the one before (its `signal`), so a slow answer to an old
 * question (a search, a count, a preview, a fetch) never lands after a newer one. A call resolves to `fn`'s result, or
 * to undefined when a newer call (or `abort`) superseded it.
 */
export interface Latest<A extends unknown[], R> {
  (...args: A): Promise<R | undefined>;
  /** Aborts the call in flight, if any. */
  abort(): void;
}

/** The latest call only (`Latest`); `dom/scoped.ts`'s `latest` also aborts it when its scope ends. */
export function createLatest<A extends unknown[], R>(
  fn: (signal: AbortSignal, ...args: A) => R | Promise<R>,
): Latest<A, R> {
  let current: AbortController | undefined;
  const call = (async (...args: A) => {
    current?.abort();
    const mine = new AbortController();
    current = mine;
    try {
      const r = await fn(mine.signal, ...args);
      return mine.signal.aborted ? undefined : r;
    } catch (e) {
      if (mine.signal.aborted) return undefined;
      throw e;
    } finally {
      if (current === mine) current = undefined;
    }
  }) as Latest<A, R>;
  call.abort = () => {
    current?.abort();
    current = undefined;
  };
  return call;
}
