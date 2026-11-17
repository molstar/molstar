/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { reportUnhandledError } from '@molstar/core/reactive/config';
import { Listened } from '@molstar/core/reactive/listened';
import { type Observable, Subscriber } from '@molstar/core/reactive/observable';
import {
  batch,
  follow as followNode,
  graphNode,
  listen,
  Sink,
  type StateNode,
  unfollow,
  versionOf,
} from '@molstar/core/reactive/state';
import type { Listenable, Unsubscribable } from '@molstar/core/reactive/types';
import { FrameLane, FrameState, flushFrame, requestFrame } from './frame.js';
import { DomHost } from './host.js';
import { DomCounts, DomProbe } from './probe.js';
import {
  getScope,
  type Live,
  type MaybeReactive,
  requireScope,
  runInScope,
  toStream,
  type ViewScope,
} from './scope.js';
import { deferTaskWork, hasWork, runWork, TaskLane } from './work.js';

// Updates (dom §2): a `watch` and a binding are sinks of the state graph (reactive §4.1). A sink is a dependent of each
// state it follows (a slot in the state's array, no object of its own) and is run in the state flush, after the change
// has spread; only a source outside the graph (a subject, a foreign stream) costs a subscriber. A marked watch goes into
// this file's dirty array and renders in a microtask (or in the frame); a binding writes in the same state flush, parents
// first (`queueBinding`).

// --- The flush -----------------------------------------------------------------------------------

// Two arrays swapped a pass, each with its count, and a flag on each watch (`_queued`): arrays kept at their size and
// emptied by count allocate nothing a flush.
let dirty: (Watch | undefined)[] = [];
let dirtyCount = 0;
let running: (Watch | undefined)[] = [];
// The frame lane's: the watches with a frame state among their sources, rendered in the frame's tick after the frame
// states publish and the microtask lane settles.
let frameDirty: (Watch | undefined)[] = [];
let frameDirtyCount = 0;
let frameRunning: (Watch | undefined)[] = [];
let frameRendering = false;
let scheduled = false;
let flushing = false;
/** Flushes in a row before a watch that keeps marking itself (a render that sets what it watches) is dropped. */
const MaxPasses = 100;

/** What is ordered by its scope's depth: a watch, a binding. */
interface Deep {
  readonly depth: number;
}

/** A counting sort's buffers (reused, grown as needed): the counts by depth, the sorted watches or bindings. */
let depthCounts = new Int32Array(64);
const sortedScratch: (Deep | undefined)[] = [];

/**
 * Puts `n` marked watches (or queued bindings) in depth order, parents first, siblings in the order marked: nothing when
 * they are already (the common case), else a stable counting sort by depth, O(n + depths).
 */
function sortByDepth<T extends Deep>(list: (T | undefined)[], n: number): void {
  let sorted = true;
  let max = 0;
  for (let i = 0; i < n; i++) {
    const d = list[i]!.depth;
    if (d > max) max = d;
    if (i > 0 && d < list[i - 1]!.depth) sorted = false;
  }
  if (sorted) return;
  if (depthCounts.length <= max + 1) depthCounts = new Int32Array(2 * (max + 1));
  const counts = depthCounts;
  counts.fill(0, 0, max + 2);
  for (let i = 0; i < n; i++) counts[list[i]!.depth + 1]++;
  for (let d = 1; d <= max + 1; d++) counts[d] += counts[d - 1];
  const out = sortedScratch as (T | undefined)[];
  for (let i = 0; i < n; i++) {
    const w = list[i]!;
    out[counts[w.depth]++] = w;
  }
  for (let i = 0; i < n; i++) {
    list[i] = out[i];
    out[i] = undefined;
  }
}

/** Renders `n` marked watches, parents first, skipping those stopped since. */
function renderMarked(list: (Watch | undefined)[], n: number): void {
  sortByDepth(list, n);
  for (let i = 0; i < n; i++) {
    const w = list[i]!;
    list[i] = undefined;
    // Stopped since it was marked (paused, disposed): its entry is stale.
    if (!w._queued) continue;
    w._queued = false;
    w._runMarked();
  }
}

function schedule(w: Watch) {
  if (w._queued) return;
  w._queued = true;
  if (w.frame) {
    frameDirty[frameDirtyCount++] = w;
    requestFrame();
    return;
  }
  dirty[dirtyCount++] = w;
  if (scheduled) return;
  scheduled = true;
  DomHost.queueMicrotask(flushQueued);
}

function flushQueued() {
  scheduled = false;
  flush();
}

/**
 * Renders now what `watch`es are waiting for (they render in a microtask otherwise): for code that sets state and then
 * reads the DOM in the same task. The task lane's work first, then the watches, parents before children; a disposed or
 * paused watch not at all; renders that mark others are rendered in the same flush. Nested in a flush, it does nothing.
 */
export function flush(): void {
  if (flushing) return;
  flushing = true;
  try {
    for (let pass = 0; dirtyCount > 0 || hasWork('task'); pass++) {
      if (pass === MaxPasses) {
        reportUnhandledError(
          new Error(`watch: renders kept marking watches (a render sets what it watches?); dropped ${dirtyCount}`),
        );
        for (let i = 0; i < dirtyCount; i++) {
          dirty[i]!._queued = false;
          dirty[i] = undefined;
        }
        dirtyCount = 0;
        // The work that fed them waits for a macrotask, or the next microtask would start the loop again.
        deferTaskWork();
        break;
      }
      runWork('task');
      if (dirtyCount === 0) continue;
      const list = dirty;
      const n = dirtyCount;
      dirty = running;
      dirtyCount = 0;
      running = list;
      renderMarked(list, n);
    }
  } finally {
    flushing = false;
  }
}

TaskLane.flush = flush;
/** The frame's tick, after the frame states published: the microtask lane, then the frame lane, once each. */
FrameLane.after = () => {
  flush();
  if (frameRendering) return;
  frameRendering = true;
  try {
    // What the frame lane's renders mark in the frame lane waits for the next frame.
    const list = frameDirty;
    const n = frameDirtyCount;
    frameDirty = frameRunning;
    frameDirtyCount = 0;
    frameRunning = list;
    renderMarked(list, n);
  } finally {
    frameRendering = false;
  }
  flush();
};

export { flushFrame };

// --- Sources -------------------------------------------------------------------------------------

/** What a `watch` follows: a state, any stream, anything with `listen` (a keyed value), a writable state. */
export type WatchSource = MaybeReactive<unknown> | Listenable<any>;

/** A source outside the graph, as followed: its change calls `changed` (no replay of a current value). */
class ChangeSubscriber extends Subscriber<unknown> {
  constructor(private readonly changed: () => void) {
    super();
  }
  protected override _next(): void {
    this.changed();
  }
  protected override _error(err: any): void {
    reportUnhandledError(err);
    this.unsubscribe();
  }
  protected override _complete(): void {
    this.unsubscribe();
  }
}

function followChanges(source: unknown, changed: () => void): Unsubscribable {
  if (source instanceof Listened) {
    // Its listener list directly: one object, and no replay.
    const s = new ChangeSubscriber(changed);
    source._add(s);
    return s;
  }
  const stop = listen(source as Listenable<unknown> | Observable<unknown>, changed);
  return { unsubscribe: stop };
}

/** A watch's sources as it follows them: a graph node, or the source itself. */
function normalize(sources: readonly WatchSource[]): unknown[] {
  const out: unknown[] = [];
  for (const s of sources) {
    const stream = toStream(s as MaybeReactive<unknown>);
    const source = stream ?? (s && typeof (s as Listenable).listen === 'function' ? s : undefined);
    if (source === undefined) continue;
    out.push(graphNode(source) ?? source);
  }
  return out;
}

function inFrameLane(sources: readonly unknown[]): boolean {
  for (let i = 0; i < sources.length; i++) if (sources[i] instanceof FrameState) return true;
  return false;
}

// --- watch ---------------------------------------------------------------------------------------

/** A `watch` whose sources can change (a recycled row shown another item). */
export interface WatchHandle {
  /** Follows `sources` instead, and renders (now while its scope is active; else on resume). */
  retarget(sources: readonly WatchSource[]): void;
}

/**
 * A watch (dom §2.2): a sink of its graph sources, a subscriber on each foreign one; renders once when it starts and
 * once per flush after a change; stopped (its scope paused), it follows nothing.
 */
class Watch extends Sink implements Live, Unsubscribable, WatchHandle {
  /** @internal Marked, waiting in a lane's array. */
  _queued = false;
  /** In the frame lane: a frame state is among its sources. */
  frame: boolean;
  private _live = false;
  private _dead = false;
  private _foreign: Unsubscribable[] | undefined = undefined;
  /** A foreign source told a change since the last render (it has no version to compare). */
  private _foreignChanged = false;
  /**
   * The graph sources' versions at the last render: one number for one source, else an array (made once). A marked
   * watch whose versions did not move (a state set and set back in a batch) does not render.
   */
  private _seen = -1;
  private _seenAll: number[] | undefined = undefined;
  private readonly _changed = () => {
    if (!this._live) return;
    this._foreignChanged = true;
    schedule(this);
  };

  constructor(
    readonly depth: number,
    private _sources: unknown[],
    private readonly _render: () => void,
  ) {
    super();
    this.frame = inFrameLane(_sources);
  }

  /** @internal Run in the state flush: a graph source changed. */
  _flush(): void {
    if (this._live) schedule(this);
  }

  retarget(sources: readonly WatchSource[]): void {
    if (this._dead) return;
    const live = this._live;
    this._stop();
    this._sources = normalize(sources);
    this.frame = inFrameLane(this._sources);
    if (live) this._start();
  }

  _start(): void {
    if (this._live || this._dead) return;
    this._live = true;
    const sources = this._sources;
    let nodes = 0;
    for (let i = 0; i < sources.length; i++) {
      const s = sources[i];
      if (s instanceof Object && graphNode(s) === s) {
        followNode(s as StateNode<unknown>, this);
        nodes++;
      } else {
        (this._foreign ??= []).push(followChanges(s, this._changed));
      }
    }
    if (nodes > 1 && (this._seenAll === undefined || this._seenAll.length !== nodes)) this._seenAll = new Array(nodes);
    else if (nodes <= 1) this._seenAll = undefined;
    this._run();
  }

  _stop(): void {
    if (!this._live) return;
    this._live = false;
    this._queued = false;
    const sources = this._sources;
    for (let i = 0; i < sources.length; i++) {
      const s = sources[i];
      if (s instanceof Object && graphNode(s) === s) unfollow(s as StateNode<unknown>, this);
    }
    const foreign = this._foreign;
    if (foreign !== undefined) {
      for (let i = 0; i < foreign.length; i++) foreign[i].unsubscribe();
      foreign.length = 0;
    }
  }

  unsubscribe(): void {
    this._dead = true;
    this._stop();
  }

  /** Whether a graph source's version moved since the last render (and records the versions now). */
  private _versionsMoved(): boolean {
    const sources = this._sources;
    let moved = false;
    let k = 0;
    const many = this._seenAll;
    for (let i = 0; i < sources.length; i++) {
      const s = sources[i];
      if (!(s instanceof Object) || graphNode(s) !== s) continue;
      const v = versionOf(s as StateNode<unknown>);
      if (many === undefined) {
        if (v !== this._seen) moved = true;
        this._seen = v;
      } else {
        if (v !== many[k]) moved = true;
        many[k] = v;
      }
      k++;
    }
    return moved;
  }

  /** @internal A marked watch's render: skipped when nothing it follows moved. */
  _runMarked(): void {
    if (!this._live) return;
    const foreign = this._foreignChanged;
    this._foreignChanged = false;
    if (!this._versionsMoved() && !foreign) return;
    this._draw();
  }

  /** @internal Renders (outside any scope: a render writes, it owns nothing), unless stopped meanwhile. */
  _run(): void {
    if (!this._live) return;
    this._foreignChanged = false;
    this._versionsMoved();
    this._draw();
  }

  private _draw(): void {
    DomCounts.render++;
    try {
      runInScope(undefined, this._render);
    } catch (e) {
      // One failing render must not keep the others in the flush from running; the watch goes on.
      reportUnhandledError(e);
    }
  }
}

/**
 * A view as a function of its state (dom §2.2): `render` reads the current values (`state.value`) and writes what
 * differs from the DOM (`writes.ts`). It runs now, and again after any of `sources` changes: the change marks the watch,
 * and every marked watch renders once in a microtask (`flush`), parents before children, so many changes in a task are
 * one render. While the scope is paused it follows nothing; resumed, it renders once. A render owns nothing; one that
 * throws is reported and the watch goes on.
 */
export function watch(sources: readonly WatchSource[], render: () => void): WatchHandle {
  const scope = requireScope('watch');
  const probe = DomProbe.on;
  const run = DomProbe.run;
  let draw = render;
  if (run) {
    const site = new Error().stack ?? '';
    draw = () => run('render', site, false, render);
  }
  if (probe) {
    const inner = draw;
    draw = () => {
      probe('render', scope);
      inner();
    };
  }
  const w = new Watch(scope.depth, normalize(sources), draw);
  scope._addLive(w);
  return w;
}

// --- Bindings ------------------------------------------------------------------------------------

// The bindings' order (dom §2.1): core runs the sinks in the order they were marked, which is their sources' dependents
// order, and a removal there swaps (reorders) it, so a binding under a `Show` could write before the `Show`'s own
// binding removes it (and read a value its content no longer has: `derive(sel, s => s!.label)` with `sel` gone). A
// binding's `_flush` only queues it; the first one queued in a flush marks `drain`, a sink of its own, which the flush
// runs after the sinks marked so far, and which writes the queued bindings parents first (by depth), skipping those
// stopped meanwhile. Two arrays swapped a drain, emptied by count: nothing allocated.
let bindingQueue: (Binding<any> | undefined)[] = [];
let bindingCount = 0;
let bindingRunning: (Binding<any> | undefined)[] = [];
let draining = false;

class BindingDrain extends Sink {
  _flush(): void {
    drainBindings();
  }
}
const drain = new BindingDrain();

function queueBinding(b: Binding<any>): void {
  bindingQueue[bindingCount++] = b;
  if (draining || drain._pending) return;
  drain._mark();
  // Not queued (outside a flush, or marked at this version already): written now, as a sink would be.
  if (!drain._pending) drainBindings();
}

function drainBindings(): void {
  if (draining) return;
  draining = true;
  try {
    // What the writes queue (a binding run by a nested flush) is drained in the same call.
    while (bindingCount > 0) {
      const list = bindingQueue;
      const n = bindingCount;
      bindingQueue = bindingRunning;
      bindingCount = 0;
      bindingRunning = list;
      sortByDepth(list, n);
      for (let i = 0; i < n; i++) {
        const b = list[i]!;
        list[i] = undefined;
        b._drain();
      }
    }
  } finally {
    draining = false;
  }
}

/** A foreign source's values, written by a binding. */
class WriteSubscriber<T> extends Subscriber<T> {
  constructor(private readonly binding: Binding<T>) {
    super();
  }
  protected override _next(value: T): void {
    this.binding._write(value);
  }
  protected override _error(err: any): void {
    reportUnhandledError(err);
    this.unsubscribe();
  }
  protected override _complete(): void {
    this.unsubscribe();
  }
}

/**
 * A binding (dom §2.1): writes a source's value now and after each change, while its scope is active. On a graph node a
 * sink (written in the state flush, once per change, the value consistent); on a foreign source with a value (a
 * `BehaviorSubject`) a subscriber in its list, then its current value; on any other stream a subscription, which
 * replays its current value if it has one.
 */
class Binding<T> extends Sink implements Live, Unsubscribable {
  private _live = false;
  private _dead = false;
  private _version = -1;
  private _sub: Unsubscribable | undefined = undefined;
  private readonly _node: StateNode<T> | undefined;

  constructor(
    /** Its scope's depth: queued bindings write parents first. */
    readonly depth: number,
    private readonly _source: Observable<T>,
    private readonly _fn: (value: T) => void,
    /** Told when it stops after it ran (its scope paused or gone). */
    private readonly _stopped?: () => void,
  ) {
    super();
    this._node = graphNode(_source) as StateNode<T> | undefined;
  }

  /** @internal */
  _write(value: T): void {
    DomCounts.write++;
    try {
      this._fn(value);
    } catch (e) {
      reportUnhandledError(e);
    }
  }

  /** @internal Run in the state flush: the node changed; written by the drain, parents first. */
  _flush(): void {
    if (this._live && this._node !== undefined) queueBinding(this);
  }

  /** @internal The drain's write: skipped when stopped since, or when the version did not move. */
  _drain(): void {
    const node = this._node;
    if (!this._live || node === undefined) return;
    // The version first (it brings the node up to date), then its cached value: one pull.
    const version = versionOf(node);
    if (version === this._version) return;
    this._version = version;
    this._write(node.value);
  }

  _start(): void {
    if (this._live || this._dead) return;
    this._live = true;
    const node = this._node;
    if (node !== undefined) {
      followNode(node, this);
      // In a batch: a set made by the first write is written after it, not nested.
      batch(() => {
        const value = node.value;
        this._version = versionOf(node);
        this._write(value);
      });
      return;
    }
    const source = this._source as Observable<T> & { readonly value?: T };
    if (source instanceof Listened && 'value' in source) {
      // Listened first, so a value told during the first write is not missed; written again if it moved meanwhile.
      const s = new WriteSubscriber(this);
      source._add(s);
      this._sub = s;
      const v = source.value as T;
      this._write(v);
      if (this._live && !Object.is(source.value, v)) this._write(source.value as T);
      return;
    }
    // The stream replays its current value on subscribe (if it has one), and the write may stop the binding meanwhile.
    const sub = source.subscribe(new WriteSubscriber(this));
    if (this._live) this._sub = sub;
    else sub.unsubscribe();
  }

  _stop(): void {
    const was = this._live;
    this._live = false;
    if (this._node !== undefined) {
      if (was) unfollow(this._node, this);
    } else {
      const sub = this._sub;
      this._sub = undefined;
      sub?.unsubscribe();
    }
    if (was) this._stopped?.();
  }

  unsubscribe(): void {
    this._dead = true;
    this._stop();
  }
}

function makeBinding<T>(
  scope: ViewScope,
  stream: Observable<T>,
  fn: (value: T) => void,
  sideEffect: boolean,
  stopped?: () => void,
) {
  const probe = DomProbe.on;
  const run = sideEffect ? undefined : DomProbe.run;
  let write = fn;
  let onStop = stopped;
  if (run) {
    const site = new Error().stack ?? '';
    // A first write since the binding started: made, or resumed (dom §9).
    let first = true;
    write = (v: T) => {
      const was = first;
      first = false;
      run('write', site, was, () => fn(v));
    };
    onStop = () => {
      first = true;
      stopped?.();
    };
  }
  if (probe) {
    const inner = write;
    write = (v: T) => {
      probe('write', scope);
      inner(v);
    };
  }
  scope._addLive(new Binding(scope.depth, stream, write, onStop));
}

/**
 * Calls `fn` with the value, and again with each new value while the current scope lives and is active. A plain value
 * is applied once. Paused, the binding follows nothing; resumed, it writes the current value.
 */
export function bind<T>(source: MaybeReactive<T>, fn: (value: T) => void): void {
  const stream = toStream(source);
  if (!stream) {
    fn(source as T);
    return;
  }
  makeBinding(requireScope(), stream, fn, false);
}

/**
 * `bind` for what keeps state per start (a region's closing, a switch's cases): `fn` is told whether the value is the
 * first since the binding started (made, or resumed after a pause), and `stopped` runs when it stops (paused or gone).
 * A plain value is applied once, as the first.
 */
export function follow<T>(
  source: MaybeReactive<T>,
  fn: (value: T, first: boolean) => void,
  stopped?: () => void,
): void {
  const stream = toStream(source);
  if (!stream) {
    fn(source as T, true);
    return;
  }
  let first = true;
  makeBinding(
    requireScope('follow'),
    stream,
    (v) => {
      const was = first;
      first = false;
      fn(v, was);
    },
    false,
    () => {
      first = true;
      stopped?.();
    },
  );
}

/**
 * Calls `fn` with each value of `source` while the scope lives and is active, as `bind` does, for a side effect that
 * writes no DOM (bringing a row into view, telling a model): the probe does not count it as a write.
 */
export function effect<T>(source: MaybeReactive<T>, fn: (value: T) => void): void {
  const stream = toStream(source);
  if (!stream) {
    fn(source as T);
    return;
  }
  makeBinding(requireScope('effect'), stream, fn, true);
}

/**
 * A frame state (`FrameState`, dom §3) for input in a component: `offer` the latest value at every event; it publishes
 * once a frame, and its pending offer goes with the current scope.
 */
export function frameState<T>(initial: T, equals?: (a: T, b: T) => boolean): FrameState<T> {
  const s = new FrameState(initial, equals);
  getScope()?.onCleanup(() => s.dispose());
  return s;
}
