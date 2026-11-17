/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  createWorkQueue,
  drag,
  flush,
  frame,
  frames,
  latest,
  observeIntersection,
  observeResize,
  on,
  pendingFrameCallbacks,
  restartable,
  runInScope,
  timeout,
  ViewScope,
  workQueue,
} from '../index.js';
import { Subject } from '@molstar/core/reactive';
import { held } from '../../state.js';
import { DomHost } from '../host.js';

// Scoped attachments, gestures and work queues (dom §6, §8).

// The page's frame through the host object: a test runs the tick itself.
const realFrame = DomHost.requestAnimationFrame;
let tick: ((time: number) => void) | undefined;

beforeAll(async () => {
  // A frame asked for by an earlier file goes first, on the real clock.
  await new Promise((r) => setTimeout(r, 20));
});

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    tick = fn;
  };
});

afterEach(() => {
  DomHost.requestAnimationFrame = realFrame;
  flushFrameTick();
  vi.useRealTimers();
});

function pointer(type: string, init: { pointerId?: number; buttons?: number; pointerType?: string } = {}) {
  const e = new Event(type, { bubbles: true });
  Object.assign(e, { pointerId: 1, buttons: 1, pointerType: 'mouse', ...init });
  return e as PointerEvent;
}

describe('scoped attachments', () => {
  it('on: detached while paused, removed when the scope ends', () => {
    const el = document.createElement('div');
    const seen = vi.fn();
    const scope = new ViewScope();
    runInScope(scope, () => on(el, 'ping', seen));
    el.dispatchEvent(new Event('ping'));
    scope.setActive(false);
    el.dispatchEvent(new Event('ping'));
    scope.setActive(true);
    el.dispatchEvent(new Event('ping'));
    scope.dispose();
    el.dispatchEvent(new Event('ping'));
    expect(seen).toHaveBeenCalledTimes(2);
  });

  it('frame runs in the next tick, after a resume when paused, never after the end', () => {
    const scope = new ViewScope();
    const seen: string[] = [];
    runInScope(scope, () => {
      frame(() => seen.push('a'));
      frame(() => seen.push('b'));
    });
    scope.setActive(false);
    flushFrameTick();
    expect(seen).toEqual([]);
    scope.setActive(true);
    flushFrameTick();
    expect(seen).toEqual(['a', 'b']);
    runInScope(scope, () => frame(() => seen.push('c')));
    scope.dispose();
    flushFrameTick();
    expect(seen).toEqual(['a', 'b']);
    expect(pendingFrameCallbacks()).toBe(0);
  });

  it('frames: requested once a frame, held while paused, cancelled at the end', () => {
    const scope = new ViewScope();
    let runs = 0;
    const f = runInScope(scope, () => frames(() => runs++));
    f.request();
    f.request();
    flushFrameTick();
    expect(runs).toBe(1);
    scope.setActive(false);
    f.request();
    flushFrameTick();
    expect(runs).toBe(1);
    scope.setActive(true);
    flushFrameTick();
    expect(runs).toBe(2);
    f.request();
    scope.dispose();
    flushFrameTick();
    expect(runs).toBe(2);
  });

  it('timeout and restartable: cancelled with the scope; a restartable waits out a pause', () => {
    vi.useFakeTimers();
    const scope = new ViewScope();
    const seen: string[] = [];
    const r = runInScope(scope, () => {
      timeout(() => seen.push('t'), 10);
      return restartable(() => seen.push('r'));
    });
    r.start(5);
    r.start(20);
    vi.advanceTimersByTime(15);
    expect(seen).toEqual(['t']);
    scope.setActive(false);
    vi.advanceTimersByTime(100);
    expect(seen).toEqual(['t']);
    scope.setActive(true);
    vi.advanceTimersByTime(20);
    expect(seen).toEqual(['t', 'r']);
    r.start(5);
    scope.dispose();
    vi.advanceTimersByTime(10);
    expect(seen).toEqual(['t', 'r']);
  });

  it('latest: only the latest call lands; the scope aborts the one in flight', async () => {
    const scope = new ViewScope();
    const signals: AbortSignal[] = [];
    const call = runInScope(scope, () =>
      latest(async (signal, x: number) => {
        signals.push(signal);
        await Promise.resolve();
        return x;
      }),
    );
    const first = call(1);
    const second = call(2);
    expect(await first).toBeUndefined();
    expect(await second).toBe(2);
    const third = call(3);
    scope.dispose();
    expect(await third).toBeUndefined();
    expect(signals[2].aborted).toBe(true);
  });
});

describe('drag', () => {
  function setup() {
    const target = document.createElement('div');
    document.body.appendChild(target);
    const scope = new ViewScope();
    const log: string[] = [];
    const d = runInScope(scope, () =>
      drag(target, {
        start: (e) => (log.push(`start ${e.pointerId}`), { n: 0 }),
        move: (_e, s) => log.push(`move ${++s.n}`),
        end: (_s, how) => log.push(`end ${how}`),
      }),
    );
    return { target, scope, log, d };
  }

  it('a press, moves and the release: up', () => {
    const { target, log, d } = setup();
    target.dispatchEvent(pointer('pointerdown'));
    target.dispatchEvent(pointer('pointermove'));
    target.dispatchEvent(pointer('pointermove', { pointerId: 2 }));
    expect(d.active).toBe(true);
    target.dispatchEvent(pointer('pointerup'));
    expect(log).toEqual(['start 1', 'move 1', 'end up']);
    expect(d.active).toBe(false);
  });

  it('a move with no button held, a lost capture with no release, a pause: cancel', () => {
    vi.useFakeTimers();
    const { target, scope, log } = setup();
    target.dispatchEvent(pointer('pointerdown'));
    target.dispatchEvent(pointer('pointermove', { buttons: 0 }));
    target.dispatchEvent(pointer('pointerdown'));
    target.dispatchEvent(pointer('lostpointercapture', { buttons: 0 }));
    vi.advanceTimersByTime(200);
    target.dispatchEvent(pointer('pointerdown'));
    scope.setActive(false);
    expect(log).toEqual(['start 1', 'end cancel', 'start 1', 'end cancel', 'start 1', 'end cancel']);
  });

  it("a capture lost just before the release: the release elsewhere in the page is the drag's", () => {
    const { target, log } = setup();
    target.dispatchEvent(pointer('pointerdown'));
    target.dispatchEvent(pointer('lostpointercapture', { buttons: 0 }));
    document.body.dispatchEvent(pointer('pointerup'));
    expect(log).toEqual(['start 1', 'end up']);
  });
});

describe('work queues', () => {
  it('runs a key once a burst (the last want), by priority, dropping what is no longer wanted', () => {
    const q = createWorkQueue({ lane: 'task' });
    const seen: string[] = [];
    q.want('a', () => seen.push('a1'));
    q.want('a', () => seen.push('a2'), { priority: 1 });
    q.want('b', () => seen.push('b'), { priority: 5 });
    q.want('c', () => seen.push('c'), { stillWanted: () => false });
    flush();
    expect(seen).toEqual(['b', 'a2']);
    expect(q.size).toBe(0);
  });

  it('a budget carries the rest to the next run; flush runs all now', () => {
    let now = 0;
    const q = createWorkQueue({ lane: 'idle', budgetMs: 5 });
    const seen: number[] = [];
    for (let i = 0; i < 4; i++)
      q.want(i, () => {
        seen.push(i);
        now += 3;
      });
    const spy = vi.spyOn(performance, 'now').mockImplementation(() => now);
    q._runNow(Number.POSITIVE_INFINITY);
    expect(seen).toEqual([0, 1]);
    q.flush();
    expect(seen).toEqual([0, 1, 2, 3]);
    spy.mockRestore();
  });

  it('a scope pauses its queue and ends it', () => {
    const scope = new ViewScope();
    const q = runInScope(scope, () => workQueue({ lane: 'task' }));
    const seen: string[] = [];
    scope.setActive(false);
    q.want('x', () => seen.push('x'));
    flush();
    expect(seen).toEqual([]);
    scope.setActive(true);
    flush();
    expect(seen).toEqual(['x']);
    scope.dispose();
    q.want('y', () => seen.push('y'));
    flush();
    expect(seen).toEqual(['x']);
  });
});

/** The frame's tick, as the page runs it: callbacks, frame states, renders. */
function flushFrameTick() {
  const t = tick;
  tick = undefined;
  t?.(0);
}

describe('review regressions', () => {
  it('a queue disposed by its own work runs nothing more; paused by it, the rest waits for the resume', () => {
    const q = createWorkQueue({ lane: 'task' });
    const seen: string[] = [];
    q.want('a', () => q.dispose(), { priority: 2 });
    q.want('b', () => seen.push('b'));
    flush();
    expect(seen).toEqual([]);
    expect(q.size).toBe(0);

    const p = createWorkQueue({ lane: 'task' });
    p.want('a', () => p.setActive(false), { priority: 2 });
    p.want('b', () => seen.push('b'));
    flush();
    expect(seen).toEqual([]);
    expect(p.size).toBe(1);
    p.setActive(true);
    flush();
    expect(seen).toEqual(['b']);
  });

  it('restartable: a start after its scope ended arms nothing', () => {
    vi.useFakeTimers();
    const fn = vi.fn();
    const scope = new ViewScope();
    const r = runInScope(scope, () => restartable(fn));
    scope.dispose();
    r.start(10);
    vi.advanceTimersByTime(20);
    expect(fn).not.toHaveBeenCalled();
  });

  it('outside a scope, timeout, observers and held throw before they make anything', () => {
    vi.useFakeTimers();
    const fn = vi.fn();
    expect(() => timeout(fn, 10)).toThrow(/outside a render scope/);
    vi.advanceTimersByTime(20);
    expect(fn).not.toHaveBeenCalled();

    let made = 0;
    class Observer {
      constructor() {
        made++;
      }
      observe() {}
      disconnect() {}
    }
    vi.stubGlobal('ResizeObserver', Observer);
    vi.stubGlobal('IntersectionObserver', Observer);
    try {
      const el = document.createElement('div');
      expect(() => observeResize(el, () => {})).toThrow(/outside a render scope/);
      expect(() => observeIntersection(el, () => {})).toThrow(/outside a render scope/);
      expect(made).toBe(0);
    } finally {
      vi.unstubAllGlobals();
    }

    const source = new Subject<number>();
    expect(() => held(source, 0)).toThrow(/outside a render scope/);
    expect(source.observed).toBe(false);
  });
});
