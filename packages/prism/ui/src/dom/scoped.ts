/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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

/** Runs `fn` in the next animation frame, unless the current scope ends first. Returns its cancel. */
export function frame(fn: (time: number) => void): () => void {
  const id = requestAnimationFrame(fn);
  let done = false;
  const cancel = () => {
    if (done) return;
    done = true;
    cancelAnimationFrame(id);
  };
  requireScope('frame').onCleanup(cancel);
  return cancel;
}

/** Runs `fn` after `ms`, unless the current scope ends first. Returns its cancel. */
export function timeout(fn: () => void, ms: number): () => void {
  const id = setTimeout(fn, ms);
  const cancel = () => clearTimeout(id);
  requireScope('timeout').onCleanup(cancel);
  return cancel;
}
