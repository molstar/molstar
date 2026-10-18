/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Observable } from 'rxjs';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { keyedSelector } from '../../reactive/keyed-selector.js';
import { mapState } from '../../reactive/map-state.js';
import { batch, observableState } from '../../reactive/observable-state.js';
import { render } from '../render.js';
import { createRoot, flush, getScope, runInScope, Scope, watch } from '../scope.js';
import { setAttr, setClass, setShown, setStyle, setText } from '../writes.js';

// `watch` (rendering.md §2): a view as render(state). The first render now; changes mark the watch, and marked watches
// render once each in a microtask (parents first), or at `flush()`. Paused, nothing; resumed, one render. Disposed,
// nothing and no listener left.

const tick = () => new Promise<void>((r) => queueMicrotask(r));

afterEach(() => {
  flush();
  vi.restoreAllMocks();
});

describe('watch', () => {
  it('renders now, then once per task however many changes', async () => {
    const a = observableState(1);
    const b = observableState('x');
    const seen: string[] = [];
    const dispose = createRoot((dispose) => {
      watch([a.asReadonly, b], () => seen.push(`${a.value}${b.value}`));
      return dispose;
    });
    expect(seen).toEqual(['1x']);
    a.set(2);
    b.set('y');
    a.set(3);
    expect(seen).toEqual(['1x']);
    await tick();
    expect(seen).toEqual(['1x', '3y']);
    batch(() => {
      a.set(4);
      b.set('z');
    });
    await tick();
    expect(seen).toEqual(['1x', '3y', '4z']);
    // A set to the same value is no change.
    a.set(4);
    await tick();
    expect(seen).toHaveLength(3);
    dispose();
  });

  it('renders parents before children, siblings in the order they changed', () => {
    const a = observableState(0);
    const b = observableState(0);
    const seen: string[] = [];
    const root = new Scope();
    const child = root.child();
    const grandchild = child.child();
    runInScope(grandchild, () => watch([a.asReadonly], () => seen.push('grandchild')));
    runInScope(child, () => watch([b.asReadonly], () => seen.push('child b')));
    runInScope(child, () => watch([a.asReadonly], () => seen.push('child a')));
    runInScope(root, () => watch([a.asReadonly, b.asReadonly], () => seen.push('root')));
    seen.length = 0;
    a.set(1);
    b.set(1);
    flush();
    expect(seen).toEqual(['root', 'child a', 'child b', 'grandchild']);
    root.dispose();
  });

  it('skips a watch disposed before the flush, and one a parent render disposes', () => {
    const a = observableState(0);
    const seen: string[] = [];
    const root = new Scope();
    const child = root.child();
    runInScope(child, () => watch([a.asReadonly], () => seen.push('child')));
    runInScope(root, () =>
      watch([a.asReadonly], () => {
        seen.push('root');
        if (a.value === 1) child.dispose();
      }),
    );
    seen.length = 0;
    a.set(1);
    flush();
    expect(seen).toEqual(['root']);
    const other = new Scope();
    runInScope(other, () => watch([a.asReadonly], () => seen.push('other')));
    seen.length = 0;
    a.set(2);
    other.dispose();
    flush();
    expect(seen).toEqual(['root']);
    root.dispose();
  });

  it('listens to nothing while paused, renders once on resume', () => {
    const a = observableState(0);
    // Counts the listeners it has.
    let listening = 0;
    const counted = Object.assign(new Observable<number>(), {
      value: 0,
      listen: (fn: () => void) => {
        listening++;
        const stop = a.asReadonly.listen!(fn);
        return () => {
          listening--;
          stop();
        };
      },
    });
    const seen: number[] = [];
    const root = new Scope();
    runInScope(root, () => watch([counted], () => seen.push(a.value)));
    expect(listening).toBe(1);
    a.set(1);
    root.setActive(false);
    expect(listening).toBe(0);
    // Marked, then paused before the flush: it does not render.
    flush();
    expect(seen).toEqual([0]);
    a.set(2);
    a.set(3);
    flush();
    expect(seen).toEqual([0]);
    root.setActive(true);
    expect(seen).toEqual([0, 3]);
    expect(listening).toBe(1);
    a.set(4);
    flush();
    expect(seen).toEqual([0, 3, 4]);
    root.dispose();
    expect(listening).toBe(0);
    a.set(5);
    flush();
    expect(seen).toEqual([0, 3, 4]);
  });

  it('made in a paused scope, renders first when it resumes', () => {
    const a = observableState(0);
    const seen: number[] = [];
    const root = new Scope();
    root.setActive(false);
    runInScope(root, () => watch([a.asReadonly], () => seen.push(a.value)));
    expect(seen).toEqual([]);
    a.set(1);
    root.setActive(true);
    expect(seen).toEqual([1]);
    root.dispose();
  });

  it('flush renders now; renders that mark others render in the same flush', () => {
    const a = observableState(0);
    const b = observableState(0);
    const seen: string[] = [];
    const root = new Scope();
    runInScope(root, () => {
      watch([a.asReadonly], () => {
        seen.push(`a${a.value}`);
        b.set(a.value * 10);
      });
      watch([b.asReadonly], () => seen.push(`b${b.value}`));
    });
    seen.length = 0;
    a.set(1);
    flush();
    expect(seen).toEqual(['a1', 'b10']);
    root.dispose();
  });

  it('a render that keeps marking itself is dropped, with an error', () => {
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    const a = observableState(0);
    const root = new Scope();
    runInScope(root, () => watch([a.asReadonly], () => a.set(a.value + 1)));
    flush();
    expect(error).toHaveBeenCalled();
    expect(a.value).toBeLessThan(200);
    root.dispose();
  });

  it('a render that throws is logged; the others render and it goes on', () => {
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    const a = observableState(0);
    const seen: string[] = [];
    const root = new Scope();
    runInScope(root, () => {
      watch([a.asReadonly], () => {
        if (a.value === 1) throw new Error('boom');
        seen.push(`first ${a.value}`);
      });
      watch([a.asReadonly], () => seen.push(`second ${a.value}`));
    });
    a.set(1);
    flush();
    expect(error).toHaveBeenCalledTimes(1);
    expect(seen).toEqual(['first 0', 'second 0', 'second 1']);
    a.set(2);
    flush();
    expect(seen).toEqual(['first 0', 'second 0', 'second 1', 'first 2', 'second 2']);
    root.dispose();
  });

  it('a render runs outside the scope: it owns nothing', () => {
    const a = observableState(0);
    const scopes: unknown[] = [];
    const root = new Scope();
    runInScope(root, () => watch([a.asReadonly], () => scopes.push(getScope())));
    a.set(1);
    flush();
    expect(scopes).toEqual([undefined, undefined]);
    root.dispose();
  });

  it('follows any stream (its replay is not a change), mapState by its source, and writable states', () => {
    const subject = new BehaviorSubject(1);
    const a = observableState(1);
    const doubled = mapState(a.asReadonly, (v) => v * 2);
    const seen: string[] = [];
    const root = new Scope();
    runInScope(root, () => watch([subject, doubled, a], () => seen.push(`${subject.value} ${doubled.value}`)));
    expect(seen).toEqual(['1 2']);
    subject.next(2);
    flush();
    a.set(2);
    flush();
    expect(seen).toEqual(['1 2', '2 2', '2 4']);
    root.dispose();
    subject.next(3);
    flush();
    expect(subject.observed).toBe(false);
    expect(seen).toHaveLength(3);
  });

  it("a delegated event's handlers end with a flush: a click's change is rendered when it returns", () => {
    const count = observableState(0);
    const container = document.createElement('div');
    document.body.append(container);
    const text = document.createTextNode('');
    const dispose = render(() => {
      const button = (
        <button type="button" onClick={() => count.set(count.value + 1)}>
          {text}
        </button>
      ) as HTMLButtonElement;
      watch([count], () => setText(text, String(count.value)));
      return button;
    }, container);
    const button = container.querySelector('button')!;
    expect(text.data).toBe('0');
    button.click();
    expect(text.data).toBe('1');
    dispose();
    container.remove();
  });

  it("a keyed selector's key tells only the rows a change concerns", () => {
    const current = observableState(0);
    const selector = keyedSelector(current.asReadonly);
    const renders = [0, 0, 0];
    const root = new Scope();
    for (let i = 0; i < 3; i++) runInScope(root, () => watch([selector.key(i)], () => renders[i]++));
    expect(renders).toEqual([1, 1, 1]);
    current.set(2);
    flush();
    expect(renders).toEqual([2, 1, 2]);
    root.dispose();
    current.set(1);
    flush();
    expect(renders).toEqual([2, 1, 2]);
  });

  it('throws outside a scope', () => {
    expect(() => watch([], () => {})).toThrow(/outside a render scope/);
  });
});

describe('listen', () => {
  it('tells after each change, after the subscribers, and stops', () => {
    const a = observableState(0);
    const order: string[] = [];
    const sub = a.asReadonly.subscribe((v) => order.push(`sub ${v}`));
    const stop = a.asReadonly.listen!(() => order.push(`listen ${a.value}`));
    a.set(1);
    a.set(1);
    expect(order).toEqual(['sub 0', 'sub 1', 'listen 1']);
    stop();
    a.set(2);
    expect(order).toEqual(['sub 0', 'sub 1', 'listen 1', 'sub 2']);
    sub.unsubscribe();
  });

  it('a listener that sets the state sets it anew', () => {
    const a = observableState(0);
    const seen: number[] = [];
    a.asReadonly.listen!(() => {
      seen.push(a.value);
      if (a.value < 3) a.set(a.value + 1);
    });
    a.set(1);
    expect(seen).toEqual([1, 2, 3]);
    expect(a.value).toBe(3);
  });
});

describe('writes', () => {
  it('write only what differs', () => {
    const el = document.createElement('div');
    const text = document.createTextNode('a');
    el.append(text);
    const observer = new MutationObserver(() => {});
    observer.observe(el, { attributes: true, characterData: true, subtree: true });
    setText(text, 'a');
    setAttr(el, 'aria-level', 2);
    setAttr(el, 'aria-level', '2');
    setAttr(el, 'data-x', undefined);
    setClass(el, 'on', true);
    setClass(el, 'on', true);
    setStyle(el, '--depth', '3');
    setStyle(el, '--depth', '3');
    setShown(el, true);
    const records = observer.takeRecords();
    expect(records.map((r) => r.attributeName ?? r.type)).toEqual(['aria-level', 'class', 'style']);
    setText(text, 'b');
    setAttr(el, 'aria-level', false);
    setStyle(el, '--depth', '');
    setShown(el, false);
    expect(text.data).toBe('b');
    expect(el.hasAttribute('aria-level')).toBe(false);
    expect(el.style.getPropertyValue('--depth')).toBe('');
    expect(el.style.display).toBe('none');
    observer.disconnect();
  });
});
