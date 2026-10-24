/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createLatest, createWorkQueue, type Latest, type WorkLane, type WorkQueue } from '../reactive/work-queue.js';
import { cancelFrameCallback, requestFrameCallback } from '../reactive/frame-state.js';
import { requireScope } from './scope.js';

// What a component attaches imperatively, released by construction (rendering.md §2.4): each helper ties what it
// attaches to the current scope, so it goes when the scope ends, and a listener is also detached while the scope is
// paused (hidden UI does nothing). The only way components attach listeners, observers, frames and timers
// (`_test/scoped-attach.test.ts` holds raw ones to a shrinking list).

/**
 * Listens to `type` on `target` while the current scope lives and is active (detached while it is paused, and for
 * good when it ends). Returns a function that removes it sooner.
 */
export function on<E extends Event = Event>(
  target: EventTarget,
  type: string,
  fn: (e: E) => void,
  options?: AddEventListenerOptions,
): () => void {
  const scope = requireScope('on');
  let attached = false;
  let gone = false;
  const listener = fn as EventListener;
  const set = (a: boolean) => {
    if (gone || a === attached) return;
    attached = a;
    if (a) target.addEventListener(type, listener, options);
    else target.removeEventListener(type, listener, options);
  };
  scope.onActive(set);
  const off = () => {
    set(false);
    gone = true;
  };
  scope.onCleanup(off);
  return off;
}

/**
 * A ResizeObserver calling `fn`, observing `el` if given (more with `.observe`), disconnected when the current scope
 * ends; undefined where there is none (jsdom).
 */
export function observeResize(el: Element | undefined, fn: ResizeObserverCallback): ResizeObserver | undefined {
  if (typeof ResizeObserver === 'undefined') return undefined;
  const observer = new ResizeObserver(fn);
  if (el) observer.observe(el);
  requireScope('observeResize').onCleanup(() => observer.disconnect());
  return observer;
}

/** An IntersectionObserver calling `fn`, observing `el`, disconnected when the current scope ends (as above). */
export function observeIntersection(
  el: Element,
  fn: IntersectionObserverCallback,
  options?: IntersectionObserverInit,
): IntersectionObserver | undefined {
  if (typeof IntersectionObserver === 'undefined') return undefined;
  const observer = new IntersectionObserver(fn, options);
  observer.observe(el);
  requireScope('observeIntersection').onCleanup(() => observer.disconnect());
  return observer;
}

/**
 * Runs `fn` in the next frame (the page's one frame tick, `requestFrameCallback`), unless the current scope ends first;
 * paused then, in the first frame after it resumes (hidden UI does nothing). Returns its cancel.
 */
export function frame(fn: (time: number) => void): () => void {
  const scope = requireScope('frame');
  let done = false;
  let id = 0;
  let waiting = false;
  let hooked = false;
  const run = (time: number) => {
    if (done || scope.disposed) return;
    if (!scope.active) {
      // Paused: wait for the resume (the hook added once, only in this rare case).
      waiting = true;
      if (!hooked) {
        hooked = true;
        scope.onActive((active) => {
          if (!active || !waiting || done) return;
          waiting = false;
          id = requestFrameCallback(run);
        });
      }
      return;
    }
    done = true;
    fn(time);
  };
  id = requestFrameCallback(run);
  const cancel = () => {
    if (done) return;
    done = true;
    cancelFrameCallback(id);
  };
  scope.onCleanup(cancel);
  return cancel;
}

/** Runs `fn` after `ms`, unless the current scope ends first. Returns its cancel. */
export function timeout(fn: () => void, ms: number): () => void {
  const id = setTimeout(fn, ms);
  const cancel = () => clearTimeout(id);
  requireScope('timeout').onCleanup(cancel);
  return cancel;
}

/** A restartable timer of the scope (`restartable`). */
export interface Restartable {
  /** Runs the function after `ms`, replacing a wait in progress. */
  start(ms: number): void;
  cancel(): void;
}

/**
 * A timer that can be started again, for a callback made in an effect (where `timeout` has no scope to bind to once the
 * call is made later, and would be one per start): one wait at a time, held while the scope is paused (it runs, in full,
 * after the scope resumes) and dropped when it ends. Hidden UI writes nothing.
 */
export function restartable(fn: () => void): Restartable {
  const scope = requireScope('restartable');
  let id: ReturnType<typeof setTimeout> | undefined;
  let ms = 0;
  let waiting = false;
  const arm = () => {
    clearTimeout(id);
    id = setTimeout(() => {
      id = undefined;
      waiting = false;
      fn();
    }, ms);
  };
  const cancel = () => {
    clearTimeout(id);
    id = undefined;
    waiting = false;
  };
  scope.onActive((active) => {
    if (active) {
      if (waiting) arm();
    } else {
      clearTimeout(id);
      id = undefined;
    }
  });
  scope.onCleanup(cancel);
  return {
    start(wait) {
      ms = wait;
      waiting = true;
      if (scope.active) arm();
    },
    cancel,
  };
}

/**
 * A queue of work by key (reactive/work-queue.ts) owned by the current scope: paused with it (hidden UI computes
 * nothing; what waited runs on resume), dropped when it ends.
 */
export function workQueue(options?: { lane?: WorkLane; budgetMs?: number }): WorkQueue {
  const q = createWorkQueue(options);
  const scope = requireScope('workQueue');
  scope.onActive((a) => q.setActive(a));
  scope.onCleanup(() => q.dispose());
  return q;
}

/** The latest call only (reactive/work-queue.ts, `Latest`), its call in flight aborted when the current scope ends. */
export function latest<A extends unknown[], R>(fn: (signal: AbortSignal, ...args: A) => R | Promise<R>): Latest<A, R> {
  const l = createLatest(fn);
  requireScope('latest').onCleanup(() => l.abort());
  return l;
}

/** A drag's handlers (`drag`). */
export interface DragHandlers<S> {
  /** A press on the target: the drag's state to start one, undefined to let it be. */
  start(e: PointerEvent): S | undefined;
  /** The pointer moved while the drag is on. */
  move?(e: PointerEvent, state: S): void;
  /**
   * It ended: released (`up`: the drag's own pointer let go, the one way to commit what a drag decides on its release, a
   * row's move, a seek), or cancelled (`pointercancel`, the capture lost, a move with no button, a new press after a lost
   * release, the scope paused or ended; `e` is the event, undefined for the scope): leave things as they were.
   */
  end?(state: S, how: 'up' | 'cancel', e: PointerEvent | undefined): void;
}

/**
 * A pointer drag on `target` (prism AGENTS.md, gestures; the scroll area's shape): its listeners attached once, with
 * the current scope (`on`: detached while it is paused, gone when it ends), acting only while a drag is on; the
 * pointer captured on the press, so the moves come to the target wherever the pointer goes. Chrome can drop the
 * capture just before the release (`lostpointercapture` with no button held, its `pointerup` a moment later to the
 * element under the pointer): a capture lost with no button held waits a moment for that release, read on the document
 * (capture phase, attached once with the scope), which ends the drag as `up`; a press anywhere meanwhile, or no release,
 * cancels it. Outside that wait the document's releases are not the drag's. It follows only the
 * pointer that started it (rendering.md §2.4): another finger's moves and releases are not its own, and its press
 * starts nothing while the first pointer still holds the capture. Every way a drag can end ends it: a release (`up`),
 * and as a cancel, `pointercancel`, the capture lost, a move with no button held (a release that never came), a new
 * press, the scope paused or ended. A press never attaches anything, so a lost release leaves nothing behind.
 */
/** How long a capture lost with no button held waits for its release (Chrome's came 1 ms after it). */
const LostReleaseMs = 100;

export function drag<S>(target: Element, handlers: DragHandlers<S>): { readonly active: boolean; cancel(): void } {
  const scope = requireScope('drag');
  let state: S | undefined;
  let pointer = -1;
  // The capture lost with no button held: the release that follows ends the drag; none after a moment, a cancel.
  let lost: ReturnType<typeof setTimeout> | undefined;
  const finish = (how: 'up' | 'cancel', e?: PointerEvent) => {
    if (lost !== undefined) clearTimeout(lost);
    lost = undefined;
    if (state === undefined) return;
    const s = state;
    state = undefined;
    if (pointer >= 0 && 'releasePointerCapture' in target) {
      try {
        if ((target as HTMLElement).hasPointerCapture?.(pointer))
          (target as HTMLElement).releasePointerCapture(pointer);
      } catch {
        // Gone already.
      }
    }
    pointer = -1;
    handlers.end?.(s, how, e);
  };
  on<PointerEvent>(target, 'pointerdown', (e) => {
    if (state !== undefined) {
      // Another finger while the first one holds the drag: not this drag's.
      if (e.pointerId !== pointer && (target as HTMLElement).hasPointerCapture?.(pointer)) return;
      // A press while one is on (its release lost): that one ends first.
      finish('cancel', e);
    }
    const s = handlers.start(e);
    if (s === undefined) return;
    state = s;
    pointer = e.pointerId;
    try {
      (target as HTMLElement).setPointerCapture?.(e.pointerId);
    } catch {
      // A synthetic or already-ended pointer: the moves still come over the target.
    }
  });
  on<PointerEvent>(target, 'pointermove', (e) => {
    if (state === undefined || e.pointerId !== pointer || lost !== undefined) return;
    // No button held: the release went missing; end here.
    if (e.buttons === 0 && e.pointerType !== 'touch') return finish('cancel', e);
    handlers.move?.(e, state);
  });
  on<PointerEvent>(target, 'pointerup', (e) => {
    if (e.pointerId === pointer) finish('up', e);
  });
  // While a capture lost with no button held waits: its release, sent elsewhere in the page, is the drag's own pointer
  // let go; a press is another gesture (the release never came).
  const doc = target.ownerDocument;
  if (doc) {
    on<PointerEvent>(
      doc,
      'pointerup',
      (e) => {
        if (lost !== undefined && e.pointerId === pointer) finish('up', e);
      },
      { capture: true },
    );
    on<PointerEvent>(
      doc,
      'pointerdown',
      (e) => {
        if (lost !== undefined) finish('cancel', e);
      },
      { capture: true },
    );
  }
  on<PointerEvent>(target, 'pointercancel', (e) => {
    if (e.pointerId === pointer) finish('cancel', e);
  });
  on<PointerEvent>(target, 'lostpointercapture', (e) => {
    if (state === undefined || e.pointerId !== pointer) return;
    // No button held: the release is coming (Chrome drops the capture just before it); a moment for it.
    if (e.buttons === 0 && e.pointerType !== 'touch') {
      lost ??= setTimeout(() => finish('cancel', e), LostReleaseMs);
      return;
    }
    finish('cancel', e);
  });
  scope.onActive((active) => {
    if (!active) finish('cancel');
  });
  scope.onCleanup(() => finish('cancel'));
  return {
    get active() {
      return state !== undefined;
    },
    cancel: () => finish('cancel'),
  };
}

/**
 * A frame callback the current scope owns, asked for from anywhere: an event handler, a binding, the frame itself, where
 * `frame()` cannot be called (no scope there). `request()` asks for `fn` in the next frame on the page's one tick
 * (once a frame, however often asked); paused, a request waits for the resume (one asked for when it paused is put off
 * too); it is cancelled when the scope ends. A loop that runs while something moves requests again from `fn`.
 */
export function frames(fn: (time: number) => void): { request(): void; cancel(): void } {
  const scope = requireScope('frames');
  let id = 0;
  let waiting = false;
  let gone = false;
  const run = (time: number) => {
    id = 0;
    if (gone) return;
    if (!scope.active) {
      waiting = true;
      return;
    }
    fn(time);
  };
  const request = () => {
    if (gone || id) return;
    if (!scope.active) {
      waiting = true;
      return;
    }
    id = requestFrameCallback(run);
  };
  const cancel = () => {
    waiting = false;
    if (!id) return;
    cancelFrameCallback(id);
    id = 0;
  };
  scope.onActive((active) => {
    if (!active) {
      if (id) {
        cancelFrameCallback(id);
        id = 0;
        waiting = true;
      }
    } else if (waiting) {
      waiting = false;
      request();
    }
  });
  scope.onCleanup(() => {
    cancel();
    gone = true;
  });
  return { request, cancel };
}
