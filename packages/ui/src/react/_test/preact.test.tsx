/** @vitest-environment jsdom */
/** @jsxImportSource preact */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Observable, of, State, Subject } from '@molstar/core/reactive';
import { createRoot } from 'preact/compat/client';
import { act } from 'preact/test-utils';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { ViewScope } from '../../dom/scope.js';
import { reactView } from '../react-view.js';
import { useReactive } from '../use-reactive.js';
import { useReactiveProp } from '../use-reactive-prop.js';

// What a Preact app's bundler does (`react` aliased to `preact/compat`): the bindings then run on Preact as written.
vi.mock('react', async () => await import('preact/compat'));

afterEach(() => {
  document.body.innerHTML = '';
});

describe('the React bindings on Preact (react aliased to preact/compat)', () => {
  it('mounts a view over a context, follows its state, and unmounts it', () => {
    const label = new State('one');
    const View = ({ ctx }: { ctx: { label: State<string> } }) => <b>{useReactive(ctx.label)}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View as any, { createRoot: createRoot as any }).mount(host, { label }, new ViewScope());
    });
    expect(host.textContent).toBe('one');
    act(() => label.set('two'));
    expect(host.textContent).toBe('two');
    act(() => dispose());
    expect(host.innerHTML).toBe('');
  });

  it('renders once on mount for a BehaviorSubject, and again only for the selected part', () => {
    const source = new BehaviorSubject({ hot: 0, cold: 'a' });
    let renders = 0;
    const View = () => {
      renders += 1;
      return <b>{useReactiveProp(source, (v) => v.cold)}</b>;
    };
    const host = document.createElement('div');
    document.body.append(host);
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View as any, { createRoot: createRoot as any }).mount(host, undefined, new ViewScope());
    });
    expect(host.textContent).toBe('a');
    expect(renders).toBe(1);
    act(() => source.next({ hot: 1, cold: 'a' }));
    expect(renders).toBe(1);
    act(() => source.next({ hot: 2, cold: 'b' }));
    expect(host.textContent).toBe('b');
    expect(renders).toBe(2);
    act(() => dispose());
  });

  it('a stream that emits while subscribing shows that value, and follows later ones (Preact checks before it subscribes)', () => {
    const later = new Subject<number>();
    const source = new Observable<number>((subscriber) => {
      subscriber.next(3);
      return later.subscribe(subscriber);
    });
    const View = ({ ctx }: { ctx: Observable<number> }) => <b>{useReactive(ctx) ?? 'none'}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View as any, { createRoot: createRoot as any }).mount(host, source, new ViewScope());
    });
    expect(host.textContent).toBe('3');
    act(() => later.next(5));
    expect(host.textContent).toBe('5');
    act(() => dispose());
    const cold = document.createElement('div');
    document.body.append(cold);
    act(() => {
      dispose = reactView(View as any, { createRoot: createRoot as any }).mount(cold, of(4), new ViewScope());
    });
    expect(cold.textContent).toBe('4');
    act(() => dispose());
  });
});
