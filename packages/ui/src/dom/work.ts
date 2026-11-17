/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { reportUnhandledError } from '@molstar/core/reactive/config';
import { FrameLane, requestFrame } from './frame.js';
import { DomHost } from './host.js';

// Work wanted, done once (dom §8): what a change asks to be computed (a tree flattened again, a list filtered) is wanted
// under a key, and a later want of the same key replaces it (the last write wins), so a burst of changes costs one run.
// A queue runs in its lane: `task` (the runtime's microtask flush, before the watches it feeds render: `flush()` runs
// it), `frame` (the animation frame, after the frame states publish, before the renders) or `idle` (when the page is
// idle). Each run first asks the wanted work whether it is `stillWanted`, orders it by `priority` (asked again each
// run) and, with a `budgetMs`, stops when the budget is spent: the rest is carried to the next run. Nothing is
// allocated a want of a known key, nor a run.

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
  /** This run's priority (asked once a run). */
  rank: number;
  wanted: boolean;
}

/** Queues with work waiting, by lane (reused, emptied by count). */
const waiting: Record<WorkLane, (WorkQueue | undefined)[]> = { task: [], frame: [], idle: [] };
const waitingCount: Record<WorkLane, number> = { task: 0, frame: 0, idle: 0 };
const spare: Record<WorkLane, (WorkQueue | undefined)[]> = { task: [], frame: [], idle: [] };
let taskScheduled = false;
let idleScheduled = false;

/** How the task lane runs: the runtime sets it (its `flush`: the task lane's work, then the watches). */
export const TaskLane: { flush: (() => void) | undefined } = { flush: undefined };

function enlist(q: WorkQueue) {
  const lane = q.lane;
  waiting[lane][waitingCount[lane]++] = q;
  if (lane === 'task') {
    if (!taskScheduled) {
      taskScheduled = true;
      DomHost.queueMicrotask(runTaskLaneQueued);
    }
  } else if (lane === 'frame') {
    requestFrame();
  } else if (!idleScheduled) {
    idleScheduled = true;
    DomHost.requestIdleCallback(runIdleLane);
  }
}

function runTaskLaneQueued() {
  taskScheduled = false;
  if (TaskLane.flush) TaskLane.flush();
  else runWork('task');
}

// The task lane yields: it runs in microtasks, so work that keeps wanting work there (a run that wants itself, two queues
// wanting each other, a budget's carry-over) would never let a timer, an event or a frame in. A queue's run whose work
// was wanted during a task-lane run is a link of its chain; past `MaxChain` links, and for a budget's carry-over at
// once, the queue waits for a macrotask, and the lane's other queues go on.
const MaxChain = 100;
let taskRunning = false;
/** Task-lane queues waiting for a macrotask (still marked enlisted: a want adds nothing). */
let deferred: WorkQueue[] = [];
let deferScheduled = false;

function defer(q: WorkQueue) {
  deferred.push(q);
  if (deferScheduled) return;
  deferScheduled = true;
  DomHost.setTimeout(resumeDeferred, 0);
}

function resumeDeferred() {
  deferScheduled = false;
  const list = deferred;
  deferred = [];
  for (const q of list) enlist(q);
}

/**
 * Sends the task lane's waiting work to a macrotask (the runtime's flush, when its renders kept marking: work that
 * feeds them must not start the loop again in the next microtask).
 */
export function deferTaskWork(): void {
  const list = waiting.task;
  const n = waitingCount.task;
  waitingCount.task = 0;
  for (let i = 0; i < n; i++) {
    defer(list[i]!);
    list[i] = undefined;
  }
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
  if (lane === 'task') taskRunning = true;
  const list = waiting[lane];
  waiting[lane] = spare[lane];
  spare[lane] = list;
  waitingCount[lane] = 0;
  const end = limitMs === undefined ? Number.POSITIVE_INFINITY : DomHost.now() + limitMs;
  try {
    for (let i = 0; i < n; i++) {
      const q = list[i]!;
      list[i] = undefined;
      q._runNow(end);
    }
  } finally {
    if (lane === 'task') taskRunning = false;
  }
}

// The frame lane's work runs in the frame's tick, after the frame states publish, before the renders.
FrameLane.work = () => runWork('frame');

/** Idle entries kept past this many are dropped at the end of a run (a queue keyed by rows that come and go). */
const KeepIdle = 512;

function noop() {}

/** A priority function's answer; one that throws is reported and ranks 0. */
function rankOf(p: () => number): number {
  try {
    return p();
  } catch (err) {
    reportUnhandledError(err);
    return 0;
  }
}

/** A queue of work by key (`workQueue` in a component: paused and ended with its scope). */
export class WorkQueue {
  /**
   * Every key's entry, wanted or idle: a key wanted again reuses its entry, and the map is not changed (deletes and sets
   * allocate as it rehashes). Idle ones past `KeepIdle` are dropped after a run.
   */
  private readonly _entries = new Map<unknown, Entry>();
  private _wanted = 0;
  /** A run's order (reused). */
  private readonly _order: (Entry | undefined)[] = [];
  private _enlisted = false;
  private _active = true;
  private _disposed = false;
  private _running = false;
  private _ranked = 0;
  /** Task lane: enlisted by a task-lane run, and the runs in a row that were (`MaxChain`). */
  private _fromRun = false;
  private _chain = 0;
  /** Ranks a wanted entry into this run's order (made once: `forEach` with it allocates no iterator). */
  private readonly _rank = (e: Entry) => {
    if (!e.wanted) return;
    const order = this._order;
    const p = e.priority;
    e.rank = p === undefined ? 0 : typeof p === 'number' ? p : rankOf(p);
    let j = this._ranked++;
    while (j > 0 && order[j - 1]!.rank < e.rank) {
      order[j] = order[j - 1];
      j--;
    }
    order[j] = e;
  };
  private readonly _dropIdle = (e: Entry, key: unknown) => {
    if (!e.wanted) this._entries.delete(key);
  };

  constructor(
    readonly lane: WorkLane = 'frame',
    /** Each run's budget (ms): the rest waits for the next run. None: all of it. */
    readonly budgetMs?: number,
  ) {}

  /** Work waiting. */
  get size(): number {
    return this._wanted;
  }

  has(key: unknown): boolean {
    return !!this._entries.get(key)?.wanted;
  }

  /** Wants `run` done in the lane's next run; a want of the same key before it runs replaces it. */
  want(key: unknown, run: () => void, options?: WantOptions): void {
    if (this._disposed) return;
    let e = this._entries.get(key);
    if (!e) {
      e = { run, priority: undefined, stillWanted: undefined, rank: 0, wanted: false };
      this._entries.set(key, e);
    }
    if (!e.wanted) {
      e.wanted = true;
      this._wanted++;
    }
    e.run = run;
    e.priority = options?.priority;
    e.stillWanted = options?.stillWanted;
    this._enlist();
  }

  /** Drops the work wanted under `key`. */
  cancel(key: unknown): void {
    const e = this._entries.get(key);
    if (e?.wanted) this._release(e);
  }

  /** Runs all the work waiting now, regardless of the budget (code that reads what it computes, a test). */
  flush(): void {
    this._runNow(Number.POSITIVE_INFINITY, true);
  }

  /** Paused, it runs nothing; resumed, what waited runs. */
  setActive(active: boolean): void {
    if (this._active === active) return;
    this._active = active;
    if (active) this._enlist();
  }

  /** Drops all the work; a run in progress stops after the entry running now (dom §8). */
  dispose(): void {
    if (this._disposed) return;
    this._disposed = true;
    this._entries.forEach(releaseEntry);
    this._entries.clear();
    this._wanted = 0;
  }

  /** @internal The lane's run: until `end` (the lane's limit), or the queue's budget. */
  _runNow(end: number, all = false): void {
    if (this.lane === 'task' && !all) {
      this._chain = this._fromRun ? this._chain + 1 : 0;
      this._fromRun = false;
      if (this._chain >= MaxChain) {
        // Still enlisted: it waits for a macrotask.
        this._chain = 0;
        reportUnhandledError(new Error('workQueue: task-lane work kept wanting work; the rest runs after a macrotask'));
        defer(this);
        return;
      }
    }
    this._enlisted = false;
    if (this._disposed || this._running || (!this._active && !all) || this._wanted === 0) return;
    this._running = true;
    let carried = false;
    try {
      const stop = all
        ? Number.POSITIVE_INFINITY
        : this.budgetMs === undefined
          ? end
          : Math.min(end, DomHost.now() + this.budgetMs);
      const order = this._order;
      this._ranked = 0;
      this._entries.forEach(this._rank);
      const n = this._ranked;
      for (let i = 0; i < n; i++) {
        // Disposed by an entry: the rest is gone. Paused by one (not for `flush`): the rest stays wanted, for the resume.
        if (this._disposed || (!this._active && !all)) break;
        const e = order[i]!;
        // Cancelled while this run ran: skip. (Wanted again: it runs what was wanted last.)
        if (!e.wanted) continue;
        // At least one a run, however long it takes.
        if (i > 0 && stop !== Number.POSITIVE_INFINITY && DomHost.now() >= stop) {
          carried = true;
          break;
        }
        const run = e.run;
        const still = e.stillWanted;
        this._release(e);
        try {
          if (still && !still()) continue;
          run();
        } catch (err) {
          reportUnhandledError(err);
        }
      }
      for (let k = 0; k < n; k++) order[k] = undefined;
      if (this._entries.size > KeepIdle + this._wanted) this._entries.forEach(this._dropIdle);
    } finally {
      this._running = false;
    }
    // Carried over (the budget spent) or wanted while running: the lane's next run; a budget's carry-over in the task
    // lane after a macrotask (a budget is there to yield).
    this._enlist(carried && this.lane === 'task');
  }

  private _enlist(yieldFirst = false) {
    if (this._enlisted || this._running || !this._active || this._disposed || this._wanted === 0) return;
    this._enlisted = true;
    if (yieldFirst) {
      defer(this);
    } else {
      this._fromRun = taskRunning && this.lane === 'task';
      enlist(this);
    }
  }

  private _release(e: Entry) {
    this._wanted--;
    releaseEntry(e);
  }
}

function releaseEntry(e: Entry) {
  e.wanted = false;
  e.run = noop;
  e.priority = undefined;
  e.stillWanted = undefined;
}

/** A queue of work by key, in `lane` (`WorkQueue`). */
export function createWorkQueue(options?: { lane?: WorkLane; budgetMs?: number }): WorkQueue {
  return new WorkQueue(options?.lane ?? 'frame', options?.budgetMs);
}

/**
 * The latest call only (dom §8): each call aborts the one before (its `signal`), so a slow answer to an old question (a
 * search, a count, a fetch) never lands after a newer one. A call resolves to `fn`'s result, or to undefined when a
 * newer call (or `abort`) superseded it.
 */
export interface Latest<A extends unknown[], R> {
  (...args: A): Promise<R | undefined>;
  /** Aborts the call in flight, if any. */
  abort(): void;
}

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
