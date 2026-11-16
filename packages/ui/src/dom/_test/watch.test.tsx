/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { ReactiveConfig } from '@molstar/core/reactive/config';
import { BehaviorSubject, batch, derive, KeyedSelector, State, Subject } from '@molstar/core/reactive';
import {
  bind,
  createRoot,
  DomCounts,
  DomProbe,
  effect,
  FrameState,
  flush,
  flushFrame,
  follow,
  frameState,
  onActive,
  onCleanup,
  ViewScope,
  runInScope,
  watch,
} from '../index.js';

// `watch`, bindings and scopes (dom §1–3).

afterEach(() => {
  ReactiveConfig.onUnhandledError = null;
  DomProbe.on = undefined;
  DomProbe.run = undefined;
});

const tick = () => Promise.resolve();

describe('watch', () => {
  it('renders now, then once a microtask however many sources changed', async () => {
    const a = new State(1);
    const b = new State('x');
    const renders: string[] = [];
    const dispose = createRoot((d) => {
      watch([a, b], () => renders.push(`${a.value}${b.value}`));
      return d;
    });
    expect(renders).toEqual(['1x']);
    a.set(2);
    b.set('y');
    a.set(3);
    expect(renders).toEqual(['1x']);
    await tick();
    expect(renders).toEqual(['1x', '3y']);
    dispose();
    a.set(4);
    await tick();
    expect(renders.length).toBe(2);
    expect(a.observed || b.observed).toBe(false);
  });

  it('is one dependent slot per state, no listener object', () => {
    const a = new State(1);
    const b = new State(2);
    createRoot(() => watch([a, b], () => {}));
    expect((a as any)._dependents.length).toBe(1);
    expect((a as any)._subs).toBeUndefined();
    expect((b as any)._dependents.length).toBe(1);
  });

  it('a diamond renders once, with consistent values', () => {
    const a = new State(1);
    const double = derive(a, (x) => x * 2);
    const seen: string[] = [];
    createRoot(() => watch([a, double], () => seen.push(`${a.value}:${double.value}`)));
    a.set(2);
    flush();
    expect(seen).toEqual(['1:2', '2:4']);
  });

  it('renders parents before children, siblings in the order marked', () => {
    const s = new State(0);
    const order: string[] = [];
    const root = new ViewScope();
    const child = root.child();
    const grandchild = child.child();
    runInScope(grandchild, () => watch([s], () => order.push('grandchild')));
    runInScope(root, () => watch([s], () => order.push('root')));
    runInScope(child, () => watch([s], () => order.push('child')));
    order.length = 0;
    s.set(1);
    flush();
    expect(order).toEqual(['root', 'child', 'grandchild']);
  });

  it('a render that sets other state renders them in the same flush; one that keeps marking itself is dropped', () => {
    const errors: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => errors.push(e);
    const a = new State(0);
    const b = new State(0);
    const seen: number[] = [];
    createRoot(() => {
      watch([a], () => b.set(a.value * 10));
      watch([b], () => seen.push(b.value));
    });
    a.set(1);
    flush();
    expect(seen).toEqual([0, 10]);
    const loop = new State(0);
    createRoot(() => watch([loop], () => loop.set(loop.value + 1)));
    loop.set(100);
    flush();
    expect(errors.length).toBe(1);
  });

  it('runs a render outside any scope, reports a throwing render and goes on', () => {
    const errors: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => errors.push(e);
    const s = new State(0);
    let inScope: unknown = 'unset';
    createRoot(() =>
      watch([s], () => {
        try {
          onCleanup(() => {});
        } catch (e) {
          inScope = e;
        }
        if (s.value === 1) throw new Error('render failed');
      }),
    );
    expect(inScope).toBeInstanceOf(Error);
    s.set(1);
    flush();
    s.set(2);
    flush();
    expect(errors).toEqual([new Error('render failed')]);
  });

  it('follows foreign sources: a subject, a BehaviorSubject, a keyed value', async () => {
    const subject = new Subject<number>();
    const behavior = new BehaviorSubject('a');
    const picked = new State('r1');
    const selector = new KeyedSelector(picked);
    let renders = 0;
    const dispose = createRoot((d) => {
      watch([subject, behavior, selector.key('r2')], () => renders++);
      return d;
    });
    subject.next(1);
    await tick();
    behavior.next('b');
    await tick();
    picked.set('r3');
    await tick();
    picked.set('r2');
    await tick();
    expect(renders).toBe(4);
    dispose();
    expect(subject.observed || behavior.observed || picked.observed).toBe(false);
  });

  it('paused, follows nothing; resumed, renders once', () => {
    const s = new State(0);
    const root = new ViewScope();
    let renders = 0;
    runInScope(root, () => watch([s], () => renders++));
    root.setActive(false);
    expect(s.observed).toBe(false);
    s.set(1);
    s.set(2);
    flush();
    expect(renders).toBe(1);
    root.setActive(true);
    expect(renders).toBe(2);
  });

  it('retarget follows other sources and renders', () => {
    const a = new State(1);
    const b = new State(2);
    let source = a;
    const seen: number[] = [];
    const handle = createRoot(() => watch([a], () => seen.push(source.value)));
    source = b;
    handle.retarget([b]);
    a.set(10);
    b.set(20);
    flush();
    expect(seen).toEqual([1, 2, 20]);
    expect(a.observed).toBe(false);
  });

  it('a watch of a frame state renders in the frame, once for many offers', () => {
    const pos = new FrameState(0);
    const model = new State('m');
    const seen: string[] = [];
    createRoot(() => watch([pos, model], () => seen.push(`${pos.value}${model.value}`)));
    for (let i = 1; i <= 10; i++) pos.offer(i);
    model.set('n');
    flush();
    expect(seen).toEqual(['0m']);
    flushFrame();
    expect(seen).toEqual(['0m', '10n']);
  });

  it('a frame state is a state: derive follows it; its pending offer goes with its scope', () => {
    const fs = new FrameState(1);
    const doubled = derive(fs, (x) => x * 2);
    fs.offer(5);
    expect(doubled.value).toBe(2);
    flushFrame();
    expect(doubled.value).toBe(10);
    let local!: FrameState<number>;
    const dispose = createRoot((d) => {
      local = frameState(0);
      return d;
    });
    local.offer(3);
    dispose();
    flushFrame();
    expect(local.value).toBe(0);
  });
});

describe('bindings', () => {
  it('bind writes the value now and after each change, once per batch; a set during the first write comes after it', () => {
    const s = new State(1);
    const seen: number[] = [];
    createRoot(() =>
      bind(s, (v) => {
        seen.push(v);
        if (v === 1) s.set(2);
      }),
    );
    batch(() => {
      s.set(3);
      s.set(4);
    });
    expect(seen).toEqual([1, 2, 4]);
  });

  it('a plain value is applied once; a foreign stream is subscribed', () => {
    const seen: unknown[] = [];
    const subject = new Subject<string>();
    createRoot(() => {
      bind(5, (v) => seen.push(v));
      bind(subject, (v) => seen.push(v));
    });
    subject.next('x');
    expect(seen).toEqual([5, 'x']);
  });

  it('follow is told the first value of each start and when it stops', () => {
    const s = new State('a');
    const root = new ViewScope();
    const log: string[] = [];
    runInScope(root, () =>
      follow(
        s,
        (v, first) => log.push(`${v}${first ? '!' : ''}`),
        () => log.push('stop'),
      ),
    );
    s.set('b');
    root.setActive(false);
    root.setActive(true);
    root.dispose();
    expect(log).toEqual(['a!', 'b', 'stop', 'b!', 'stop']);
  });

  it('effect runs like bind and is not a write for the probe', () => {
    const kinds: string[] = [];
    DomProbe.run = (kind, _site, _first, fn) => {
      kinds.push(kind);
      fn();
    };
    const s = new State(0);
    let effects = 0;
    createRoot(() => {
      effect(s, () => effects++);
      bind(s, () => {});
    });
    s.set(1);
    expect(effects).toBe(2);
    expect(kinds).toEqual(['write', 'write']);
  });

  it('counts writes and renders', () => {
    const before = { ...DomCounts };
    const s = new State(0);
    createRoot(() => {
      bind(s, () => {});
      watch([s], () => {});
    });
    s.set(1);
    flush();
    expect(DomCounts.write - before.write).toBe(2);
    expect(DomCounts.render - before.render).toBe(2);
  });
});

describe('ViewScope', () => {
  it('disposes children first (newest first), then cleanups in reverse; a failing cleanup is reported', () => {
    const errors: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => errors.push(e);
    const order: string[] = [];
    const root = new ViewScope();
    root.onCleanup(() => order.push('root 1'));
    root.onCleanup(() => {
      throw new Error('cleanup failed');
    });
    root.onCleanup(() => order.push('root 3'));
    const a = root.child();
    a.onCleanup(() => order.push('a'));
    const b = root.child();
    b.onCleanup(() => order.push('b'));
    root.dispose();
    expect(order).toEqual(['b', 'a', 'root 3', 'root 1']);
    expect(errors).toEqual([new Error('cleanup failed')]);
    expect(root.disposed && a.disposed && b.disposed).toBe(true);
  });

  it('a child disposed alone leaves its parent; a child of a disposed scope is disposed', () => {
    const root = new ViewScope();
    const kids = Array.from({ length: 3 }, () => root.child());
    kids[1].dispose();
    expect((root as any)._firstChild).toBe(kids[0]);
    expect((kids[0] as any)._next).toBe(kids[2]);
    root.dispose();
    expect(root.child().disposed).toBe(true);
  });

  it('onActive tells the state now and on every change, pausing hooks before bindings', () => {
    const root = new ViewScope();
    const log: string[] = [];
    const s = new State(0);
    runInScope(root, () => {
      onActive((a) => log.push(`hook ${a}`));
      bind(s, () => log.push('write'));
    });
    root.setActive(false);
    root.setActive(true);
    expect(log).toEqual(['hook true', 'write', 'hook false', 'write', 'hook true']);
  });

  it('follow: the scope is active while a state is true', () => {
    const shown = new State(true);
    const root = new ViewScope();
    root.follow(shown);
    shown.set(false);
    expect(root.active).toBe(false);
    shown.set(true);
    expect(root.active).toBe(true);
  });

  it('requires a scope for reactive bindings and watches', () => {
    expect(() => bind(new State(1), () => {})).toThrow(/outside a render scope/);
    expect(() => watch([], () => {})).toThrow(/outside a render scope/);
  });

  it('probe.on is told writes and renders with their scope', () => {
    const on = vi.fn();
    DomProbe.on = on;
    const s = new State(0);
    createRoot(() => {
      bind(s, () => {});
      watch([s], () => {});
    });
    expect(on.mock.calls.map((c) => c[0])).toEqual(['write', 'render']);
  });
});

describe('a set back in a batch', () => {
  it('neither writes a binding nor renders a watch', () => {
    const s = new State(0);
    const other = new Subject<number>();
    let writes = 0;
    let renders = 0;
    createRoot(() => {
      bind(s, () => writes++);
      watch([s], () => renders++);
    });
    batch(() => {
      s.set(1);
      s.set(0);
    });
    flush();
    expect(writes).toBe(1);
    expect(renders).toBe(1);
    // A foreign source has no version: its change always renders.
    createRoot(() => watch([s, other], () => renders++));
    other.next(1);
    flush();
    expect(renders).toBe(3);
  });

  it('a state followed only by a binding is observed', () => {
    const s = new State(0);
    const dispose = createRoot((d) => {
      bind(s, () => {});
      return d;
    });
    expect(s.observed).toBe(true);
    dispose();
    expect(s.observed).toBe(false);
  });
});
