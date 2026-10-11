/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, skip } from 'rxjs';

import { observableState } from '../observable-state.js';

type Equal<A, B> = (<T>() => T extends A ? 1 : 2) extends <T>() => T extends B ? 1 : 2 ? true : false;

/** Compile-time only assertion (no runtime effect); a mismatch is a type error where the test is type-checked. */
function expectTypeOf<A>(_value: A) {
  return { toEqualTypeOf<B>(..._check: Equal<A, B> extends true ? [] : [never]) {} };
}

describe('observableState', () => {
  it('keeps one readonly projection and suppresses shallow-equal replacements for every subscriber', () => {
    const state = observableState({ count: 0, label: 'one' });
    const initial = state.value;
    const listener = jest.fn();
    const subscription = state.asReadonly.subscribe(listener);
    state.set({ count: 0, label: 'one' });
    expect(listener).toHaveBeenCalledTimes(1);
    expect(state.value).toBe(initial);
    state.set({ count: 1, label: 'one' });
    expect(listener).toHaveBeenCalledTimes(2);
    expect(state.value).toEqual({ count: 1, label: 'one' });
    expect(state.asReadonly.value).toBe(state.value);
    expect('set' in state.asReadonly).toBe(false);
    expect('next' in state.asReadonly).toBe(false);
    subscription.unsubscribe();
  });

  it('checks partial fields before merging, and retains fields not supplied by a patch', () => {
    const readLabel = jest.fn(() => 'one');
    const initial = {
      count: 0,
      get label() {
        return readLabel();
      },
    };
    const state = observableState(initial);
    state.update({ count: 0 });
    expect(readLabel).not.toHaveBeenCalled();
    expect(state.value).toBe(initial);
    state.update({ count: 1 });
    expect(state.value).toEqual({ count: 1, label: 'one' });
    const changed = state.value;
    state.update({ label: 'one' });
    expect(state.value).toBe(changed);
  });

  it('supports custom equality on replacements and changed patches', () => {
    const state = observableState(
      { locator: { id: 'a' }, count: 0 },
      {
        compare: (left, right) => left.locator.id === right.locator.id,
      },
    );
    const initial = state.value;
    state.set({ locator: { id: 'a' }, count: 1 });
    state.update({ count: 2 });
    expect(state.value).toBe(initial);
    state.update({ locator: { id: 'b' } });
    expect(state.value.locator.id).toBe('b');
  });

  it('holds scalar, array, undefined, and class values without exposing record patching', () => {
    const array = observableState<readonly number[]>([1, 2]);
    const initial = array.value;
    array.set([1, 2]);
    expect(array.value).toBe(initial);
    expect('update' in array).toBe(false);
    const scalar = observableState<number | undefined>(undefined);
    scalar.set(1);
    scalar.set(undefined);
    expect(scalar.value).toBeUndefined();
    expectTypeOf(scalar.value).toEqualTypeOf<number | undefined>();
    class Example {
      value = 1;
    }
    const model = observableState(new Example());
    const replacement = new Example();
    model.set(replacement);
    expect(model.value).toBe(replacement);
    expect('update' in model).toBe(false);
  });
});

// A set inside a subscriber reaches every subscriber in order
describe('a set inside a subscriber reaches every subscriber in order', () => {
  describe('a set() inside a subscriber reaches later subscribers out of order', () => {
    it("a later subscriber's last delivered value equals state.value", () => {
      const state = observableState(0);
      // A: reverts to the value it insists on.
      state.asReadonly.pipe(skip(1)).subscribe((v) => {
        if (v === 1) state.set(2);
      });
      const seen: number[] = [];
      state.asReadonly.subscribe((v) => seen.push(v));

      state.set(1);

      expect(state.value).toBe(2);
      expect(seen.at(-1)).toBe(state.value);
    });

    it('a later distinctUntilChanged subscriber ends on the current value', () => {
      const state = observableState('A');
      state.asReadonly.pipe(skip(1), distinctUntilChanged()).subscribe((inView) => {
        if (inView !== 'A') state.set('A');
      });
      // Like a view model's `shown` state: it keeps the last key delivered to it.
      let shown = 'A';
      state.asReadonly.pipe(distinctUntilChanged()).subscribe((key) => {
        shown = key;
      });

      state.set('B');

      expect(state.value).toBe('A');
      expect(shown).toBe('A');
    });
  });
});

// A queued value that is itself undefined is still delivered.
describe('a queued undefined', () => {
  it('is delivered', () => {
    const state = observableState<number | undefined>(0);
    state.asReadonly.subscribe((v) => {
      if (v === 1) state.set(undefined);
    });
    const seen: (number | undefined)[] = [];
    state.asReadonly.subscribe((v) => seen.push(v));
    state.set(1);
    expect(seen).toEqual([0, 1, undefined]);
    expect(state.value).toBeUndefined();
  });
});
