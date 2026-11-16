/** @vitest-environment jsdom */
/** @jsxImportSource react */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  BehaviorSubject,
  derive,
  map,
  Observable,
  of,
  ReactiveConfig,
  ReactiveModel,
  type ReadonlyState,
  State,
  Subject,
} from '@molstar/core/reactive';
import { act, StrictMode, useSyncExternalStore } from 'react';
import { createRoot, type Root } from 'react-dom/client';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import type { Child } from '../../dom/index.js';
import { defineElements, ViewElement } from '../../dom/element.js';
import { ViewScope } from '../../dom/scope.js';
import { UiElementHost } from '../element-host.js';
import { reactView } from '../react-view.js';
import { useReactive } from '../use-reactive.js';
import { useReactiveModel } from '../use-reactive-model.js';
import { useReactiveProp } from '../use-reactive-prop.js';

type Equal<A, B> = (<T>() => T extends A ? 1 : 2) extends <T>() => T extends B ? 1 : 2 ? true : false;

/** Compile-time only assertion (no runtime effect); a mismatch is a type error where the test is type-checked. */
function expectTypeOf<A>(_value: A) {
  return { toEqualTypeOf<B>(..._check: Equal<A, B> extends true ? [] : [never]) {} };
}

(globalThis as typeof globalThis & { IS_REACT_ACT_ENVIRONMENT: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

let container: HTMLDivElement;
let root: Root;
beforeEach(() => {
  container = document.createElement('div');
  document.body.append(container);
  root = createRoot(container);
});
afterEach(async () => {
  await act(async () => root.unmount());
  container.remove();
});

describe('useReactive: subscribing tells React nothing', () => {
  it('a state, a derived value and a BehaviorSubject render once on mount, then once per change', async () => {
    const name = new State('one');
    const upper = derive(name, (n) => n.toUpperCase());
    const subject = new BehaviorSubject(1);
    let renders = 0;
    function Probe() {
      renders += 1;
      return (
        <span>
          {useReactive(name)}/{useReactive(upper)}/{useReactive(subject)}
        </span>
      );
    }
    await act(async () => root.render(<Probe />));
    expect(container.textContent).toBe('one/ONE/1');
    expect(renders).toBe(1);
    await act(async () => name.set('two'));
    expect(container.textContent).toBe('two/TWO/1');
    expect(renders).toBe(2);
    await act(async () => subject.next(2));
    expect(container.textContent).toBe('two/TWO/2');
    expect(renders).toBe(3);
  });

  it('the store React subscribes to is told on changes only, never by the replay on subscribing', async () => {
    // The hook's store, seen from React's side: count what React's `onStoreChange` receives.
    const told = vi.fn();
    const state = new State(1);
    const subject = new BehaviorSubject(1);
    function Spy({ source }: { source: ReadonlyState<number> }) {
      // A probe store around the hook's own: the same subscribe semantics (react §1), counted.
      useSyncExternalStore(
        (onChange) => {
          const unlisten = source.listen(() => {
            told();
            onChange();
          });
          return unlisten;
        },
        () => source.value,
      );
      return <span>{useReactive(source)}</span>;
    }
    await act(async () =>
      root.render(
        <StrictMode>
          <Spy source={state} />
          <Spy source={subject} />
        </StrictMode>,
      ),
    );
    expect(told).not.toHaveBeenCalled();
    await act(async () => state.set(2));
    await act(async () => subject.next(2));
    expect(told).toHaveBeenCalledTimes(2);
    expect(container.textContent).toBe('22');
  });

  it('under Strict Mode a state renders twice for the double render, never a third time for a replay', async () => {
    const state = new State('a');
    const subject = new BehaviorSubject('b');
    let renders = 0;
    function Probe() {
      renders += 1;
      return (
        <span>
          {useReactive(state)}
          {useReactive(subject)}
        </span>
      );
    }
    await act(async () =>
      root.render(
        <StrictMode>
          <Probe />
        </StrictMode>,
      ),
    );
    expect(container.textContent).toBe('ab');
    expect(renders).toBe(2);
  });

  it('a value changed between render and subscribe is caught when React subscribes', async () => {
    const state = new State(1);
    function Probe() {
      const value = useReactive(state);
      // Changed after the render read it, before the effect subscribes.
      if (value === 1) queueMicrotask(() => state.set(5));
      return <span>{value}</span>;
    }
    await act(async () => root.render(<Probe />));
    expect(container.textContent).toBe('5');
  });
});

describe('useReactive', () => {
  it('uses defaults only to seed a source, resets on replacement, and releases replaced subscriptions', async () => {
    const first = new Subject<number>();
    const second = new Subject<number>();
    const releases = vi.fn();
    const source = new Observable<number>((subscriber) => {
      const subscription = first.subscribe(subscriber);
      return () => {
        subscription.unsubscribe();
        releases();
      };
    });
    function Probe({ source, seed }: { source: Observable<number> | null; seed: number }) {
      return <span>{useReactive(source, { defaultValue: seed })}</span>;
    }
    await act(async () => root.render(<Probe source={source} seed={1} />));
    expect(container.textContent).toBe('1');
    await act(async () => root.render(<Probe source={source} seed={2} />));
    expect(container.textContent).toBe('1');
    await act(async () => first.next(5));
    expect(container.textContent).toBe('5');
    await act(async () => root.render(<Probe source={second} seed={2} />));
    expect(container.textContent).toBe('2');
    expect(releases).toHaveBeenCalledTimes(1);
    await act(async () => first.next(9));
    expect(container.textContent).toBe('2');
    await act(async () => root.render(<Probe source={null} seed={7} />));
    expect(container.textContent).toBe('7');
  });

  it('suppresses shallow-equal whole values and keeps the previous snapshot reference', async () => {
    const source = new Subject<{ count: number }>();
    const seen: Array<{ count: number } | undefined> = [];
    function Probe() {
      const value = useReactive(source);
      seen.push(value);
      return <span>{value?.count}</span>;
    }
    await act(async () => root.render(<Probe />));
    const value = { count: 1 };
    await act(async () => source.next(value));
    const renders = seen.length;
    await act(async () => source.next({ count: 1 }));
    expect(seen).toHaveLength(renders);
    await act(async () => root.render(<Probe />));
    expect(seen.at(-1)).toBe(value);
  });

  it('keeps the last value of a BehaviorSubject that ended in an error', async () => {
    const source = new BehaviorSubject({ count: 1 });
    const selections: Array<{ label: string }> = [];
    function Raw() {
      return <span>{useReactive(source).count}</span>;
    }
    function Selected() {
      const selection = useReactiveProp(source, (value) => ({ label: String(value.count) }));
      selections.push(selection);
      return <span>{selection.label}</span>;
    }
    const render = () =>
      root.render(
        <>
          <Raw />
          <Selected />
        </>,
      );
    const previous = ReactiveConfig.onUnhandledError;
    const handleError = vi.fn();
    ReactiveConfig.onUnhandledError = handleError;
    try {
      await act(async () => render());
      await act(async () => source.next({ count: 2 }));
      const stable = selections.at(-1);
      await act(async () => {
        source.error(new Error('source failed'));
        await new Promise((resolve) => setTimeout(resolve, 0));
      });
      await act(async () => render());
      expect(container.textContent).toBe('22');
      expect(selections.at(-1)).toBe(stable);
    } finally {
      ReactiveConfig.onUnhandledError = previous;
    }
  });

  it('receives the first emission from a stream that does not replay', async () => {
    const changes = new Subject<number>();
    function Probe() {
      return <span>{useReactive(changes, { defaultValue: 1 })}</span>;
    }
    await act(async () => root.render(<Probe />));
    await act(async () => changes.next(2));
    expect(container.textContent).toBe('2');
  });

  it('keeps what a stream emits while subscribing, under Strict Mode, without leaking subscriptions', async () => {
    const started = vi.fn();
    const released = vi.fn();
    const source = new Observable<number>((subscriber) => {
      started();
      subscriber.next(3);
      return released;
    });
    function Probe({ source }: { source: Observable<number> }) {
      return <span>{useReactive(source)}</span>;
    }
    await act(async () =>
      root.render(
        <StrictMode>
          <Probe source={source} />
        </StrictMode>,
      ),
    );
    expect(container.textContent).toBe('3');
    expect(started.mock.calls.length - released.mock.calls.length).toBe(1);
    await act(async () =>
      root.render(
        <StrictMode>
          <Probe source={of(4)} />
        </StrictMode>,
      ),
    );
    expect(container.textContent).toBe('4');
    expect(released).toHaveBeenCalledTimes(started.mock.calls.length);
  });
});

describe('useReactiveProp', () => {
  it('renders again only when the selected part changes', async () => {
    const source = new State({ hot: 0, cold: 'a' });
    let renders = 0;
    function Probe() {
      renders += 1;
      return <span>{useReactiveProp(source, (value) => value.cold)}</span>;
    }
    await act(async () => root.render(<Probe />));
    expect(renders).toBe(1);
    await act(async () => source.set({ hot: 1, cold: 'a' }));
    expect(renders).toBe(1);
    await act(async () => source.set({ hot: 2, cold: 'b' }));
    expect(renders).toBe(2);
    expect(container.textContent).toBe('b');
  });

  it('compares selected objects shallowly', async () => {
    const list = { available: [] as readonly string[], id: undefined as string | undefined };
    const source = new State({ list, view: { count: 0 } });
    let renders = 0;
    function Probe() {
      renders += 1;
      return <span>{useReactiveProp(source, (value) => value.list).id ?? 'none'}</span>;
    }
    await act(async () => root.render(<Probe />));
    await act(async () => source.set({ list, view: { count: 1 } }));
    expect(renders).toBe(1);
    await act(async () => source.set({ list: { ...list, id: 'c-1' }, view: { count: 2 } }));
    expect(renders).toBe(2);
    expect(container.textContent).toBe('c-1');
  });

  it('selects the replayed value once even when the comparer always says changed', async () => {
    const source = new BehaviorSubject(1);
    const selector = vi.fn((value: number) => ({ value }));
    let renders = 0;
    function Probe() {
      renders += 1;
      return <span>{useReactiveProp(source, selector, { compare: () => false }).value}</span>;
    }
    await act(async () => root.render(<Probe />));
    expect(container.textContent).toBe('1');
    expect(renders).toBe(1);
    expect(selector).toHaveBeenCalledTimes(1);
    await act(async () => source.next(2));
    expect(container.textContent).toBe('2');
    expect(renders).toBe(2);
    expect(selector).toHaveBeenCalledTimes(2);
  });

  it('does not select absent input, but does select defaults and emitted undefined', async () => {
    const source = new Subject<number | undefined>();
    const selector = vi.fn((value: number | undefined) => (value === undefined ? 'unavailable' : String(value)));
    function Probe() {
      return <span>{useReactiveProp(source, selector)}</span>;
    }
    await act(async () => root.render(<Probe />));
    expect(selector).not.toHaveBeenCalled();
    await act(async () => source.next(undefined));
    expect(container.textContent).toBe('unavailable');
    function MissingProbe() {
      return <span>{useReactiveProp(null, selector, { defaultValue: 2 })}</span>;
    }
    await act(async () => root.render(<MissingProbe />));
    expect(container.textContent).toBe('2');
  });

  it('keeps its selection and its one listener across inline selectors and comparers', async () => {
    const source = new State({ hot: 0, label: 'one' });
    const listen = vi.spyOn(source, 'listen');
    const selections: Array<{ label: string }> = [];
    function Probe() {
      const selection = useReactiveProp(source, (value) => ({ label: value.label }), {
        compare: (left, right) => left.label === right.label,
      });
      selections.push(selection);
      return <span>{selection.label}</span>;
    }
    await act(async () => root.render(<Probe />));
    const first = selections[0];
    await act(async () => source.set({ hot: 1, label: 'one' }));
    expect(selections).toHaveLength(1);
    await act(async () => root.render(<Probe />));
    expect(selections.at(-1)).toBe(first);
    await act(async () => source.set({ hot: 2, label: 'two' }));
    expect(container.textContent).toBe('two');
    expect(listen).toHaveBeenCalledTimes(1);
  });

  it("renders the new selection when only the selector's prop changes", async () => {
    const state = new State({ left: 'L', right: 'R' });
    function Probe({ field }: { field: 'left' | 'right' }) {
      return <span>{useReactiveProp(state, (value) => value[field])}</span>;
    }
    await act(async () => root.render(<Probe field="left" />));
    expect(container.textContent).toBe('L');
    await act(async () => root.render(<Probe field="right" />));
    expect(container.textContent).toBe('R');
  });

  it('infers current-value, defaulted, nullable and selected types', () => {
    function useCheckOverloads(
      source: ReadonlyState<number>,
      nullable: ReadonlyState<number> | undefined,
      plain: Observable<number>,
    ) {
      expectTypeOf(useReactive(source)).toEqualTypeOf<number>();
      expectTypeOf(useReactive(nullable)).toEqualTypeOf<number | undefined>();
      expectTypeOf(useReactive(plain, { defaultValue: 0 })).toEqualTypeOf<number>();
      expectTypeOf(useReactive(plain, {})).toEqualTypeOf<number | undefined>();
      expectTypeOf(useReactiveProp(source, String)).toEqualTypeOf<string>();
      expectTypeOf(useReactiveProp(nullable, String)).toEqualTypeOf<string | undefined>();
      expectTypeOf(useReactiveProp(plain, String, { defaultValue: 0 })).toEqualTypeOf<string>();
    }
    expectTypeOf<ReturnType<typeof useCheckOverloads>>(undefined).toEqualTypeOf<void>();
  });
});

describe('an unstable source', () => {
  it('a stable pipe shows the piped value and follows it', async () => {
    const base = new BehaviorSubject(5);
    const doubled = base.pipe(map((x) => x * 2));
    function Probe() {
      return <span>{useReactive(doubled, { defaultValue: 0 })}</span>;
    }
    await act(async () => root.render(<Probe />));
    expect(container.textContent).toBe('10');
    await act(async () => base.next(7));
    expect(container.textContent).toBe('14');
  });

  it('a new source on renders in a row warns once; one that changes once does not', async () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const base = new State(5);
    function Probe({ n }: { n: number }) {
      return <span>{useReactive(derive(base, (x) => x * 2 + n))}</span>;
    }
    for (const n of [0, 1, 2, 3]) await act(async () => root.render(<Probe n={n} />));
    expect(warn.mock.calls.filter(([m]) => String(m).includes('new on every render'))).toHaveLength(1);
    warn.mockClear();
    const a = new State(1);
    const b = new State(2);
    function Pick({ which }: { which: 'a' | 'b' }) {
      return <span>{useReactive(which === 'a' ? a : b)}</span>;
    }
    await act(async () => root.render(<Pick which="a" />));
    await act(async () => root.render(<Pick which="b" />));
    await act(async () => root.render(<Pick which="b" />));
    expect(container.textContent).toBe('2');
    expect(warn).not.toHaveBeenCalled();
  });
});

describe('useReactiveModel', () => {
  it('mounts and unmounts with the component, Strict Mode included, and never retires the model', async () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      constructor() {
        super();
        this.registerPermanentDispose(() => order.push('permanent cleanup'));
      }
      protected override onMount(): void {
        order.push('mount');
        this.registerDispose(() => order.push('mount cleanup'));
      }
      protected override onRetire(): void {
        order.push('retire');
      }
    }
    const model = new Example();
    function Probe() {
      useReactiveModel(model);
      return null;
    }
    await act(async () =>
      root.render(
        <StrictMode>
          <Probe />
        </StrictMode>,
      ),
    );
    expect(order).toEqual(['mount', 'mount cleanup', 'mount']);
    await act(async () => root.render(null));
    expect(order).toEqual(['mount', 'mount cleanup', 'mount', 'mount cleanup']);
    expect(order).not.toContain('retire');
    model.retire();
    expect(order.slice(-2)).toEqual(['retire', 'permanent cleanup']);
  });

  it('returns undefined for no model', async () => {
    let seen = 'pending';
    function Probe({ model }: { model: ReactiveModel | null }) {
      seen = useReactiveModel(model) === undefined ? 'undefined' : 'present';
      return null;
    }
    await act(async () => root.render(<Probe model={null} />));
    expect(seen).toBe('undefined');
  });
});

class Greeting extends ViewElement<{ name: string }> {
  static readonly tag = 'msui-test-react-greeting';
  protected render(model: { name: string }): Child {
    const span = document.createElement('span');
    span.textContent = `Hello ${model.name}`;
    return span;
  }
}
defineElements([Greeting]);

describe('views and elements', () => {
  it('mounts a React view over a context, follows its states, and unmounts it', () => {
    const label = new State('one');
    const View = ({ ctx }: { ctx: { label: State<string> } }) => <b>{useReactive(ctx.label)}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View, { createRoot }).mount(host, { label }, new ViewScope());
    });
    expect(host.textContent).toBe('one');
    act(() => label.set('two'));
    expect(host.textContent).toBe('two');
    act(() => dispose());
    expect(host.innerHTML).toBe('');
    host.remove();
  });

  it('tells the view whether it is active, following its scope', () => {
    const View = ({ active }: { ctx: undefined; active?: boolean }) => <b>{String(active)}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    const scope = new ViewScope();
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View, { createRoot }).mount(host, undefined, scope);
    });
    expect(host.textContent).toBe('true');
    act(() => scope.setActive(false));
    expect(host.textContent).toBe('false');
    act(() => scope.setActive(true));
    expect(host.textContent).toBe('true');
    act(() => dispose());
    host.remove();
  });

  it('hosts a UI element with its model set as a property', () => {
    const model = { name: 'Mol*' };
    act(() => root.render(<UiElementHost tag="msui-test-react-greeting" properties={{ model }} />));
    const el = container.querySelector('msui-test-react-greeting') as Greeting;
    expect(el.model).toBe(model);
    expect(el.textContent).toBe('Hello Mol*');
  });
});
