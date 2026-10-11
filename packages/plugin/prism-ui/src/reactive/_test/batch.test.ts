/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, observableState } from '../observable-state.js';

describe('batch', () => {
  it('emits each state once with its final value, in the order first set; reads inside see the new values', () => {
    const a = observableState(0);
    const b = observableState('x');
    const seen: string[] = [];
    a.asReadonly.subscribe((v) => seen.push(`a${v}`));
    b.asReadonly.subscribe((v) => seen.push(`b${v}`));
    seen.length = 0;
    batch(() => {
      a.set(1);
      b.set('y');
      a.set(2);
      expect(a.value).toBe(2);
      expect(a.asReadonly.value).toBe(2);
      expect(seen).toEqual([]);
    });
    expect(seen).toEqual(['a2', 'by']);
  });

  it('emits nothing for a state set back to where it started, and holds until the outermost batch ends', () => {
    const a = observableState(5);
    const seen: number[] = [];
    a.asReadonly.subscribe((v) => seen.push(v));
    seen.length = 0;
    batch(() => {
      a.set(6);
      a.set(5);
    });
    expect(seen).toEqual([]);
    batch(() => {
      batch(() => a.set(7));
      expect(seen).toEqual([]);
      a.set(8);
    });
    expect(seen).toEqual([8]);
  });

  it('delivers what a subscriber sets while the batch is delivered', () => {
    const a = observableState(0);
    const b = observableState(0);
    a.asReadonly.subscribe((v) => b.set(v * 10));
    const seen: number[] = [];
    b.asReadonly.subscribe((v) => seen.push(v));
    seen.length = 0;
    batch(() => a.set(3));
    expect(b.value).toBe(30);
    expect(seen).toEqual([30]);
  });
});
