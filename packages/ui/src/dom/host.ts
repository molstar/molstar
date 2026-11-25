/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The browser APIs the runtime schedules with (dom §3), read through one object: a test or a host without a browser
 * replaces them. The defaults read the globals at each call, so fake timers that patch the globals drive them too.
 */
export const DomHost: {
  queueMicrotask(fn: () => void): void;
  /** The page's frame. Without `requestAnimationFrame` (Node, a worker), a 16 ms timer. */
  requestAnimationFrame(fn: (time: number) => void): void;
  /** Idle time. Without `requestIdleCallback`, a 1 ms timer. */
  requestIdleCallback(fn: (deadline?: { timeRemaining(): number }) => void): void;
  setTimeout(fn: () => void, ms: number): unknown;
  clearTimeout(id: unknown): void;
  now(): number;
} = {
  queueMicrotask: (fn) => queueMicrotask(fn),
  requestAnimationFrame: (fn) => {
    if (typeof requestAnimationFrame === 'function') requestAnimationFrame(fn);
    else setTimeout(() => fn(DomHost.now()), 16);
  },
  requestIdleCallback: (fn) => {
    if (typeof requestIdleCallback === 'function') requestIdleCallback(fn, { timeout: 200 });
    else setTimeout(fn, 1);
  },
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (id) => clearTimeout(id as ReturnType<typeof setTimeout>),
  now: () => (typeof performance !== 'undefined' ? performance.now() : Date.now()),
};
