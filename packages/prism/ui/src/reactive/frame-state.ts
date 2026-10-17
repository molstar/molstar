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

const raf: (fn: () => void) => unknown =
  typeof requestAnimationFrame === 'function' ? (fn) => requestAnimationFrame(fn) : (fn) => setTimeout(fn, 16);

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
  raf(frameTick);
}

function frameTick(): void {
  requested = false;
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
