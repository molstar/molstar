/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { reportUnhandledError } from '@molstar/core/reactive/config';
import { batch, State } from '@molstar/core/reactive/state';
import { DomHost } from './host.js';

// The frame lane (dom §3): input comes faster than frames (a trackpad's wheel and pointer events at 120 Hz or more), so
// input-fed state is batched at its source. A frame state's `offer` keeps only the latest value; the first offer in a
// frame asks for the page's one animation frame, whose tick runs the frame callbacks, publishes every offered frame
// state once, runs the frame lane's work (`FrameLane.work`) and then what the runtime renders (`FrameLane.after`). One
// `requestAnimationFrame` a frame for all of them; nothing allocated an offer or a tick.

/** The page's frame requests and the callbacks run (an idle page asks for none). */
export const FrameCounts = { requested: 0, callbacks: 0 };

/** What runs in the frame's tick after the frame states publish: the frame lane's work, then the renders. */
export const FrameLane: { work: (() => void) | undefined; after: (() => void) | undefined } = {
  work: undefined,
  after: undefined,
};

let requested = false;

/** Asks for the frame's tick (once a frame). */
export function requestFrame(): void {
  if (requested) return;
  requested = true;
  FrameCounts.requested++;
  DomHost.requestAnimationFrame(frameTick);
}

// Frame callbacks: what a component would ask `requestAnimationFrame` for joins the one tick, run in the order asked,
// before the frame states publish. One asked for while they run waits for the next frame. Reused arrays, emptied by
// count.
let callbacks: (((time: number) => void) | undefined)[] = [];
let callbackIds: number[] = [];
let callbackCount = 0;
let spareCallbacks: (((time: number) => void) | undefined)[] = [];
let spareIds: number[] = [];
let nextId = 1;
/** The frame's callbacks while they run (their ids are `spareIds`), and the one running: a cancel reaches the rest. */
let runningCallbacks: (((time: number) => void) | undefined)[] | undefined;
let runningCount = 0;
let runningAt = 0;

/** Runs `fn` in the next frame (the page's one tick); returns an id for `cancelFrameCallback`. */
export function requestFrameCallback(fn: (time: number) => void): number {
  const id = nextId++;
  callbacks[callbackCount] = fn;
  callbackIds[callbackCount++] = id;
  requestFrame();
  return id;
}

/** Frame callbacks asked for and not yet run nor cancelled. */
export function pendingFrameCallbacks(): number {
  let n = 0;
  for (let i = 0; i < callbackCount; i++) if (callbacks[i]) n++;
  return n;
}

/** Drops a frame callback not yet run, also one of the frame running now (cancelled by a callback before it). */
export function cancelFrameCallback(id: number): void {
  for (let i = 0; i < callbackCount; i++) {
    if (callbackIds[i] !== id) continue;
    callbacks[i] = undefined;
    // Cancelled at the end (asked for again at the next event, the common case): dropped, so events between two frames
    // keep one entry, not one each.
    while (callbackCount > 0 && callbacks[callbackCount - 1] === undefined) callbackCount--;
    return;
  }
  const list = runningCallbacks;
  if (!list) return;
  for (let i = runningAt + 1; i < runningCount; i++) {
    if (spareIds[i] === id) {
      list[i] = undefined;
      return;
    }
  }
}

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
      reportUnhandledError(e);
    }
  }
  runningCallbacks = undefined;
  runningCount = 0;
}

/** The frame states offered in this frame (reused, emptied by count). */
let offered: (FrameState<any> | undefined)[] = [];
let offeredCount = 0;
let spareOffered: (FrameState<any> | undefined)[] = [];

/** The frame states being published, and how many (the batch's function reads them: no closure a tick). */
let publishing: (FrameState<any> | undefined)[] = [];
let publishingCount = 0;

function publishAll(): void {
  const list = publishing;
  const n = publishingCount;
  publishingCount = 0;
  for (let i = 0; i < n; i++) {
    const s = list[i]!;
    list[i] = undefined;
    s._publish();
  }
}

function publishOffered(): void {
  // Offers made while publishing (a listener offering again) go to the next frame.
  const list = offered;
  const n = offeredCount;
  offered = spareOffered;
  offeredCount = 0;
  spareOffered = list;
  if (n === 0) return;
  publishing = list;
  publishingCount = n;
  // In one batch (dom §3): a derived value or a binding over several frame states sees them all published, never some.
  batch(publishAll);
}

function frameTick(time: number): void {
  requested = false;
  runCallbacks(time);
  publishOffered();
  FrameLane.work?.();
  FrameLane.after?.();
}

/**
 * Publishes the offered frame states and renders what follows, now (tests; code that must read the DOM the input leads
 * to). Not `flush()`: that one runs at every delegated event, and would publish per event.
 */
export function flushFrame(): void {
  publishOffered();
  FrameLane.work?.();
  FrameLane.after?.();
}

/**
 * A state fed by input (dom §3): `offer` keeps the latest value (no listener told, nothing allocated) and it is
 * published once a frame, a `set` then (compared, told to its listeners, marked in the state graph: `derive` follows it
 * as any state). `value` is the published value; `latest` the last offered. A `watch` following one renders in the
 * frame lane. `dispose` drops a pending offer.
 */
export class FrameState<T> extends State<T> {
  /** The marker the runtime reads: a watch following one renders in the frame lane. */
  readonly frameLane = true;
  private _next: T;
  private _offered = false;
  private _disposed = false;

  constructor(initial: T, equals?: (a: T, b: T) => boolean) {
    super(initial, equals);
    this._next = initial;
  }

  /** The last offered value (not yet published, maybe). */
  get latest(): T {
    return this._offered ? this._next : this.value;
  }

  /** Keeps `v` to publish in the next frame. */
  offer(v: T): void {
    if (this._disposed) return;
    this._next = v;
    if (this._offered) return;
    this._offered = true;
    offered[offeredCount++] = this;
    requestFrame();
  }

  /** @internal Publishes the offered value (the frame's tick). */
  _publish(): void {
    if (!this._offered) return;
    this._offered = false;
    const v = this._next;
    this._next = undefined as T;
    this.set(v);
  }

  dispose(): void {
    this._disposed = true;
    this._offered = false;
    this._next = undefined as T;
  }
}
