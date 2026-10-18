/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ObjectMutator } from '../object-mutator.js';

function source() {
  return {
    kind: 'root',
    children: [
      { kind: 'a', params: { opacity: 1, point: [1, 2, 3] } },
      { kind: 'b', params: { color: 'red' }, custom: { x: { y: 1 } } },
    ],
  };
}

describe('ObjectMutator', () => {
  it('copies the paths of changed slots only and keeps the source unchanged', () => {
    const src = source();
    const plan = ObjectMutator.plan(src);
    plan.pushPrefix('children', 0, 'params');
    const opacity = plan.slot('opacity');
    const point0 = plan.slot(['point', 0]);
    plan.popPrefix();
    const m = ObjectMutator.create(plan);

    m.set(opacity, 0.5);
    m.set(point0, 10);
    const r = m.commit();
    expect(r).not.toBe(src);
    expect(r.children[0].params).toEqual({ opacity: 0.5, point: [10, 2, 3] });
    expect(r.children[1]).toBe(src.children[1]);
    expect(src).toEqual(source());
  });

  it('returns the previous root when nothing changed, and shares unchanged subtrees with it', () => {
    const src = source();
    const plan = ObjectMutator.plan(src);
    const opacity = plan.slot(['children', 0, 'params', 'opacity']);
    const y = plan.slot(['children', 1, 'custom', 'x', 'y']);
    const m = ObjectMutator.create(plan);

    m.set(opacity, 1);
    expect(m.commit()).toBe(src); // the source's value: no change

    m.set(opacity, 0.5);
    m.set(y, 2);
    const r1 = m.commit();
    m.set(opacity, 0.5);
    m.set(y, 2);
    expect(m.commit()).toBe(r1);

    m.set(opacity, 0.25);
    m.set(y, 2);
    const r2 = m.commit();
    expect(r2.children[0].params.opacity).toBe(0.25);
    expect(r2.children[1]).toBe(r1.children[1]);
    expect(r2.children[0].params.point).toBe(src.children[0].params.point);
  });

  it('compares values deeply', () => {
    const plan = ObjectMutator.plan(source());
    const point = plan.slot(['children', 0, 'params', 'point']);
    const m = ObjectMutator.create(plan);
    m.set(point, [4, 5, 6]);
    const r = m.commit();
    m.set(point, [4, 5, 6]);
    expect(m.commit()).toBe(r);
  });

  it('a slot not set takes the source value again, and a key the source does not have is removed', () => {
    const src = source();
    const plan = ObjectMutator.plan(src);
    const opacity = plan.slot(['children', 0, 'params', 'opacity']);
    const size = plan.slot(['children', 0, 'params', 'size']);
    const m = ObjectMutator.create(plan);
    expect(m.get(size)).toBeUndefined();

    m.set(opacity, 0);
    m.set(size, 3);
    const r = m.commit();
    expect(r.children[0].params).toEqual({ opacity: 0, point: [1, 2, 3], size: 3 });

    const back = m.commit();
    expect(back.children[0].params).toEqual(src.children[0].params);
    expect('size' in back.children[0].params).toBe(false);
  });

  it('get returns the value set in this frame, else the source value', () => {
    const plan = ObjectMutator.plan(source());
    const opacity = plan.slot(['children', 0, 'params', 'opacity']);
    const m = ObjectMutator.create(plan);
    expect(m.get(opacity)).toBe(1);
    m.set(opacity, 0.5);
    expect(m.get(opacity)).toBe(0.5);
    m.commit();
    expect(m.get(opacity)).toBe(1);
  });

  it('ignores slots whose path is missing in the source', () => {
    const src = source();
    const plan = ObjectMutator.plan(src);
    const missing = plan.slot(['children', 0, 'custom', 'x']);
    const m = ObjectMutator.create(plan);
    m.set(missing, 1);
    expect(m.commit()).toBe(src);
    expect(m.get(missing)).toBeUndefined();
  });

  it('applies slots inside other slots in order, from the source each frame', () => {
    const src = source();
    const plan = ObjectMutator.plan(src);
    plan.pushPrefix('children', 0, 'params');
    const point = plan.slot('point');
    const point1 = plan.slot(['point', 1]);
    plan.popPrefix();
    expect(plan.isDisjoint).toBe(false);
    const m = ObjectMutator.create(plan);

    const value = [7, 8, 9];
    m.set(point, value);
    expect(m.get(point1)).toBe(8);
    m.set(point1, 0);
    const r = m.commit();
    expect(r.children[0].params.point).toEqual([7, 0, 9]);
    expect(value).toEqual([7, 8, 9]);

    m.set(point, [7, 8, 9]);
    m.set(point1, 0);
    expect(m.commit()).toBe(r);

    m.set(point1, 5); // the outer slot is not set: the source's point
    expect(m.commit().children[0].params.point).toEqual([1, 5, 3]);
    expect(src).toEqual(source());
  });
});
