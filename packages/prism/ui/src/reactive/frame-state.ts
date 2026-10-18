/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from './readonly-state.js';
import { Listened } from './state.js';

// The frame lane (rendering.md §2.3): input comes faster than frames (a trackpad's wheel and pointer events at 120 Hz
// or more, several a frame), so input-fed state is batched at its source. A frame state's `offer` only keeps the
// latest value; the first offer in a frame asks for one animation frame, which publishes every offered frame state
// once, runs the frame lane's work queues (`FrameLane.work`, reactive/work-queue.ts), then lets the runtime render what
// follows from them (`FrameLane.after`: the microtask lane's watches, then the frame lane's, which are the watches with
// a frame state among their sources). One `requestAnimationFrame` a frame for all of them; nothing is allocated an
// offer.

const raf: (fn: (time: number) => void) => unknown =
  typeof requestAnimationFrame === 'function'
    ? (fn) => requestAnimationFrame(fn)
    : (fn) => setTimeout(() => fn(performance.now()), 16);

/** The page's frame requests and the callbacks run: counted for the idle and interaction budget. */
export const FrameCounts = { requested: 0, callbacks: 0 };

/** The frame states offered in this frame (reused, emptied by count). */
let offered: (FrameState<unknown> | undefined)[] = [];
let offeredCount = 0;
let requested = false;

/** What the runtime runs in the frame, after the frame states publish (dom/scope.ts: both lanes' watches). */
export const FrameLane: { work: (() => void) | undefined; after: (() => void) | undefined } = {
  work: undefined,
  after: undefined,
};

/** Asks for the frame's tick (once a frame). */
export function requestFrame(): void {
  if (requested) return;
  requested = true;
  FrameCounts.requested++;
  raf(frameTick);
}

// Frame callbacks (rendering.md §2.3, one frame tick): what a component would ask `requestAnimationFrame` for (a
// placement, a canvas's draw, a simulation's step) joins the page's one frame request, run in the order asked, before
// the frame states publish. One asked for while they run waits for the next frame. Reused arrays, emptied by count.
let callbacks: (((time: number) => void) | undefined)[] = [];
let callbackIds: number[] = [];
let callbackCount = 0;
let spareCallbacks: (((time: number) => void) | undefined)[] = [];
let spareIds: number[] = [];
let nextId = 1;

/** Runs `fn` in the next frame (the page's one frame tick); returns an id for `cancelFrameCallback`. */
export function requestFrameCallback(fn: (time: number) => void): number {
  const id = nextId++;
  callbacks[callbackCount] = fn;
  callbackIds[callbackCount++] = id;
  requestFrame();
  return id;
}

/** Frame callbacks asked for and not yet run nor cancelled (tests: a view disposed asks for none). */
export function pendingFrameCallbacks(): number {
  let n = 0;
  for (let i = 0; i < callbackCount; i++) if (callbacks[i]) n++;
  return n;
}

/** Drops a frame callback not yet run, also one of the frame running now (cancelled by a callback before it). */
export function cancelFrameCallback(id: number): void {
  for (let i = 0; i < callbackCount; i++)
    if (callbackIds[i] === id) {
      callbacks[i] = undefined;
      // Cancelled at the end (asked for again at the next event, the common case): dropped, so a run of events between
      // two frames keeps one entry, not one an event.
      while (callbackCount > 0 && callbacks[callbackCount - 1] === undefined) callbackCount--;
      return;
    }
  const list = runningCallbacks;
  if (!list) return;
  for (let i = runningAt + 1; i < runningCount; i++)
    if (spareIds[i] === id) {
      list[i] = undefined;
      return;
    }
}

/** The frame's callbacks while they run (their ids are `spareIds`), and the one running: a cancel reaches the rest. */
let runningCallbacks: (((time: number) => void) | undefined)[] | undefined;
let runningCount = 0;
let runningAt = 0;

function runCallbacks(time: number): void {
  const list = callbacks;
  const ids = callbackIds;
  const n = callbackCount;
  if (n === 0) return;
  callbacks = spareCallbacks;
  callbackIds = spareIds;
  callbackCount = 0;
  spareCallbacks = list;
  spareIds = ids;
  runningCallbacks = list;
  runningCount = n;
  for (let i = 0; i < n; i++) {
    runningAt = i;
    const fn = list[i];
    list[i] = undefined;
    if (!fn) continue;
    FrameCounts.callbacks++;
    try {
      fn(time);
    } catch (e) {
      console.error(e);
    }
  }
  runningCallbacks = undefined;
  runningCount = 0;
}

function frameTick(time: number): void {
  requested = false;
  runCallbacks(time);
  publishOffered();
  FrameLane.work?.();
  FrameLane.after?.();
}

function publishOffered(): void {
  // Offers made while publishing (a listener offering again) go to the next frame.
  const list = offered;
  const n = offeredCount;
  offered = spare;
  offeredCount = 0;
  spare = list;
  for (let i = 0; i < n; i++) {
    const s = list[i]!;
    list[i] = undefined;
    s.publish();
  }
}
let spare: (FrameState<unknown> | undefined)[] = [];

/**
 * Publishes the offered frame states and renders what follows, now (tests; code that must read the DOM the input
 * leads to). Not the runtime's `flush()`: that one runs at every delegated event, and would publish per event.
 */
export function flushFrame(): void {
  publishOffered();
  FrameLane.work?.();
  FrameLane.after?.();
}

/**
 * A state fed by input (rendering.md §2.3, the frame lane): `offer` keeps the latest value (no listener told, nothing
 * allocated), and it is published once a frame, told to its listeners when it differs (`compare`). `value` is the
 * published value (what its views show); `latest` the last offered. `dispose` drops a pending offer (`frameState` in
 * dom/scope.ts does it when its scope ends).
 */
export class FrameState<T> extends Listened<T> implements ReadonlyState<T> {
  /** The marker the runtime reads: a watch following one renders in the frame lane. */
  readonly frameLane = true;
  #value: T;
  #next: T;
  #pending = false;
  #disposed = false;
  readonly #compare: (a: T, b: T) => boolean;

  constructor(initial: T, compare: (a: T, b: T) => boolean = Object.is) {
    super();
    this.#value = initial;
    this.#next = initial;
    this.#compare = compare;
  }

  get value(): T {
    return this.#value;
  }

  /** The last offered value (not yet published, maybe). */
  get latest(): T {
    return this.#pending ? this.#next : this.#value;
  }

  /** Keeps `v` to publish in the next frame. */
  offer(v: T): void {
    if (this.#disposed) return;
    this.#next = v;
    if (this.#pending) return;
    this.#pending = true;
    offered[offeredCount++] = this as FrameState<unknown>;
    requestFrame();
  }

  /** @internal Publishes the offered value (the frame's tick). */
  publish(): void {
    if (!this.#pending) return;
    this.#pending = false;
    const v = this.#next;
    if (this.#compare(this.#value, v)) return;
    this.#value = v;
    this.tell(v);
  }

  dispose(): void {
    this.#disposed = true;
    this.#pending = false;
  }
}

/** A frame state (`FrameState`); `dom/scope.ts`'s `frameState` also drops it with its scope. */
export function createFrameState<T>(initial: T, compare?: (a: T, b: T) => boolean): FrameState<T> {
  return new FrameState(initial, compare);
}
