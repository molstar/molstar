/** @vitest-environment jsdom */
/** @jsxImportSource react */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act, StrictMode } from 'react';
import { createRoot, type Root } from 'react-dom/client';
import { BehaviorSubject, config, map, Observable, of, Subject } from 'rxjs';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { useObservable } from '../use-observable.js';
import { useObservableProp } from '../use-observable-prop.js';

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
  vi.restoreAllMocks();
});

describe('useObservable', () => {
  it('reads current values without a second mount render and supports both behavior implementations', async () => {
    let renders = 0;
    function Probe({ source }: { source: ReadonlyState<string> }) {
      renders += 1;
      const value = useObservable(source);
      return <span>{value}</span>;
    }
    const source = observableState('one');
    await act(async () => root.render(<Probe source={source.asReadonly} />));
    expect(container.textContent).toBe('one');
    expect(renders).toBe(1);
    await act(async () => source.set('two'));
    expect(container.textContent).toBe('two');
    const subject = new BehaviorSubject('four');
    await act(async () => root.render(<Probe source={subject} />));
    expect(container.textContent).toBe('four');
  });

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
      return <span>{useObservable(source, { defaultValue: seed })}</span>;
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
      const value = useObservable(source);
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

  it('retains raw and selected snapshots when an errored current-value source throws on a later render', async () => {
    const source = new BehaviorSubject({ count: 1 });
    const selections: Array<{ label: string }> = [];
    function Raw() {
      return <span>{useObservable(source).count}</span>;
    }
    function Selected() {
      const selection = useObservableProp(source, (value) => ({ label: String(value.count) }));
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
    const previousErrorHandler = config.onUnhandledError;
    const handleError = vi.fn();
    config.onUnhandledError = handleError;
    try {
      await act(async () => render());
      await act(async () => source.next({ count: 2 }));
      const stable = selections.at(-1);
      const error = new Error('source failed');
      await act(async () => {
        source.error(error);
        await new Promise((resolve) => setTimeout(resolve, 0));
      });
      expect(handleError).toHaveBeenCalledWith(error);
      await act(async () => render());
      expect(container.textContent).toBe('22');
      expect(selections.at(-1)).toBe(stable);
    } finally {
      config.onUnhandledError = previousErrorHandler;
    }
  });

  it('receives the first emission from an ordinary Observable that does not replay', async () => {
    const changes = new Subject<number>();
    const source = changes.asObservable();
    function Probe() {
      return <span>{useObservable(source, { defaultValue: 1 })}</span>;
    }
    await act(async () => root.render(<Probe />));
    await act(async () => changes.next(2));
    expect(container.textContent).toBe('2');
  });

  it('handles synchronous emissions and completion under Strict Mode without leaking subscriptions', async () => {
    const started = vi.fn();
    const released = vi.fn();
    const source = new Observable<number>((subscriber) => {
      started();
      subscriber.next(3);
      return released;
    });
    function Probe({ source }: { source: Observable<number> }) {
      return <span>{useObservable(source)}</span>;
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

describe('useObservableProp', () => {
  it('skips current-value replay before selection even when comparison always reports a change', async () => {
    const source = new BehaviorSubject(1);
    const selector = vi.fn((value: number) => ({ value }));
    let renders = 0;
    function Probe() {
      renders += 1;
      const selection = useObservableProp(source, selector, { compare: () => false });
      return <span>{selection.value}</span>;
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
    function Probe({ seed }: { seed: boolean }) {
      return <span>{useObservableProp(source, selector, seed ? { defaultValue: 1 } : undefined)}</span>;
    }
    await act(async () => root.render(<Probe seed={false} />));
    expect(selector).not.toHaveBeenCalled();
    await act(async () => source.next(undefined));
    expect(container.textContent).toBe('unavailable');
    const behavior = new BehaviorSubject<number | undefined>(undefined);
    function DefaultProbe() {
      return <span>{useObservableProp(behavior, selector, { defaultValue: 1 })}</span>;
    }
    await act(async () => root.render(<DefaultProbe />));
    expect(container.textContent).toBe('unavailable');
    function MissingProbe() {
      return <span>{useObservableProp(null, selector, { defaultValue: 2 })}</span>;
    }
    await act(async () => root.render(<MissingProbe />));
    expect(container.textContent).toBe('2');
  });

  it('retains selection and subscription across inline selector and comparator identities', async () => {
    const source = new BehaviorSubject({ hot: 0, label: 'one' });
    const subscribe = vi.spyOn(source, 'subscribe');
    const selections: Array<{ label: string }> = [];
    function Probe() {
      const selection = useObservableProp(source, (value) => ({ label: value.label }), {
        compare: (left, right) => left.label === right.label,
      });
      selections.push(selection);
      return <span>{selection.label}</span>;
    }
    await act(async () => root.render(<Probe />));
    const first = selections[0];
    await act(async () => source.next({ hot: 1, label: 'one' }));
    expect(selections).toHaveLength(1);
    await act(async () => root.render(<Probe />));
    expect(selections.at(-1)).toBe(first);
    await act(async () => source.next({ hot: 2, label: 'two' }));
    expect(container.textContent).toBe('two');
    expect(subscribe).toHaveBeenCalledTimes(1);
  });

  it('uses current selector and comparator refs on emissions without resubscribing', async () => {
    const changes = new Subject<{ left: number; right: number }>();
    const started = vi.fn();
    const source = new Observable<{ left: number; right: number }>((subscriber) => {
      started();
      return changes.subscribe(subscriber);
    });
    function Probe({ field, suppress }: { field: 'left' | 'right'; suppress: boolean }) {
      const value = useObservableProp(source, (value) => value[field], { compare: suppress ? () => true : Object.is });
      return <span>{value}</span>;
    }
    await act(async () => root.render(<Probe field="left" suppress={true} />));
    await act(async () => changes.next({ left: 1, right: 2 }));
    await act(async () => changes.next({ left: 3, right: 4 }));
    expect(container.textContent).toBe('1');
    await act(async () => root.render(<Probe field="right" suppress={false} />));
    await act(async () => changes.next({ left: 5, right: 6 }));
    expect(container.textContent).toBe('6');
    expect(started).toHaveBeenCalledTimes(1);
  });

  it('resets selection on source replacement even when the raw value is equal', async () => {
    const first = new Subject<number>();
    const second = new Subject<number>();
    function Probe({ source }: { source: Observable<number> }) {
      return <span>{useObservableProp(source, (value) => value * 2)}</span>;
    }
    await act(async () => root.render(<Probe source={first} />));
    await act(async () => first.next(2));
    expect(container.textContent).toBe('4');
    await act(async () => root.render(<Probe source={second} />));
    expect(container.textContent).toBe('');
    await act(async () => second.next(2));
    expect(container.textContent).toBe('4');
  });

  it('infers current-value, defaulted, nullable, and selected return types', () => {
    function useCheckOverloads(
      source: ReadonlyState<number>,
      nullable: ReadonlyState<number> | undefined,
      plain: Observable<number>,
    ) {
      expectTypeOf(useObservable(source)).toEqualTypeOf<number>();
      expectTypeOf(useObservable(nullable)).toEqualTypeOf<number | undefined>();
      expectTypeOf(useObservable(plain, { defaultValue: 0 })).toEqualTypeOf<number>();
      expectTypeOf(useObservable(plain, {})).toEqualTypeOf<number | undefined>();
      expectTypeOf(useObservableProp(source, String)).toEqualTypeOf<string>();
      expectTypeOf(useObservableProp(nullable, String)).toEqualTypeOf<string | undefined>();
      expectTypeOf(useObservableProp(plain, String, { defaultValue: 0 })).toEqualTypeOf<string>();
      expectTypeOf(useObservableProp(source, () => undefined)).toEqualTypeOf<undefined>();
      // @ts-expect-error a different result type requires a selector
      useObservable<number, string>(source);
      // @ts-expect-error a raw default must match the source value, not the selection
      useObservableProp<number, string>(plain, String, { defaultValue: 'zero' });
    }
    expectTypeOf<ReturnType<typeof useCheckOverloads>>(undefined).toEqualTypeOf<void>();
  });
});

// A new selector selects again
describe('a new selector selects again', () => {
  type Pair = { left: string; right: string; other: number };

  describe('useObservableProp with a prop-dependent selector', () => {
    it("renders the new selection when only the selector's prop changes", async () => {
      const state = observableState<Pair>({ left: 'L', right: 'R', other: 0 });
      function Probe({ field }: { field: 'left' | 'right' }) {
        return <span>{useObservableProp(state.asReadonly, (value) => value[field])}</span>;
      }
      await act(async () => root.render(<Probe field="left" />));
      expect(container.textContent).toBe('L');
      await act(async () => root.render(<Probe field="right" />));
      expect(container.textContent).toBe('R');
    });
  });
});

// A source's identity is its lifetime; an inline source is detected, not followed.
describe('an unstable source', () => {
  it('a stable pipe shows the piped value and follows it', async () => {
    const base = new BehaviorSubject(5);
    const doubled = base.pipe(map((x) => x * 2));
    function Probe() {
      return <span>{useObservable(doubled, { defaultValue: 0 })}</span>;
    }
    await act(async () => root.render(<Probe />));
    expect(container.textContent).toBe('10');
    await act(async () => base.next(7));
    expect(container.textContent).toBe('14');
  });

  it('a new source on renders in a row warns once', async () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const base = new BehaviorSubject(5);
    function Probe({ n }: { n: number }) {
      return <span>{useObservable(base.pipe(map((x) => x * 2 + n)), { defaultValue: 0 })}</span>;
    }
    for (const n of [0, 1, 2, 3]) await act(async () => root.render(<Probe n={n} />));
    expect(warn.mock.calls.filter(([m]) => String(m).includes('new observable on every render'))).toHaveLength(1);
  });

  it('a source that changes once does not warn', async () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const a = observableState(1);
    const b = observableState(2);
    function Probe({ which }: { which: 'a' | 'b' }) {
      return <span>{useObservable(which === 'a' ? a.asReadonly : b.asReadonly)}</span>;
    }
    await act(async () => root.render(<Probe which="a" />));
    await act(async () => root.render(<Probe which="b" />));
    await act(async () => root.render(<Probe which="b" />));
    expect(container.textContent).toBe('2');
    expect(warn).not.toHaveBeenCalled();
  });
});
