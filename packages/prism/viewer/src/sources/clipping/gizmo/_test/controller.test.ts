/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import { type CutObject, placeCut } from '../../cuts.js';
import { followCuts } from '../controller.js';

// The active cut as the cuts change: known by its number, which shifts when one before it goes (clipping.md §4).

const object = (tag?: string, x = 0): CutObject => ({
  ...placeCut('sphere', Sphere3D.create(Vec3.create(x, 0, 0), 10), Vec3.create(0, 0, 1)),
  tag,
});
const cuts = (...objects: CutObject[]) => objects.map((o, i) => ({ index: i + 1, object: o }));

describe('following the active cut', () => {
  it('follows a tagged cut to its new number when one before it goes', () => {
    const [a, b] = [object('a'), object('b')];
    const next = followCuts(new Map([[2, b]]), cuts(b), 2);
    expect([...next]).toEqual([[1, b]]);
    expect(next.get(1)).not.toBe(a);
  });

  it('follows no other cut when the active one is removed', () => {
    const [a, b] = [object('a'), object('b')];
    expect(followCuts(new Map([[1, a]]), cuts(b), 2).size).toBe(0);
  });

  it('follows an edit: by its tag, its object the new one', () => {
    const edited = object('a', 5);
    expect([...followCuts(new Map([[1, object('a')]]), cuts(edited), 1)]).toEqual([[1, edited]]);
  });

  it('follows an untagged cut by its number when it was edited, not when the count changed', () => {
    const [old, edited] = [object(undefined, 1), object(undefined, 5)];
    expect([...followCuts(new Map([[1, old]]), cuts(edited), 1)]).toEqual([[1, edited]]);
    expect(followCuts(new Map([[1, old]]), cuts(edited), 2).size).toBe(0);
  });

  it('follows an untagged cut by its values to a new number', () => {
    const [other, active] = [object(undefined, 1), object(undefined, 5)];
    expect([...followCuts(new Map([[2, active]]), cuts(active), 2)]).toEqual([[1, active]]);
    expect(followCuts(new Map([[2, active]]), cuts(other, active), 2).get(2)).toBe(active);
  });

  it('takes no cut for an edit while the old object is still among them', () => {
    const [a, b] = [object('a'), object('b')];
    // number 1 holds another cut now, and the active one is number 2
    expect([...followCuts(new Map([[1, a]]), cuts(b, a), 2)]).toEqual([[2, a]]);
  });

  it('follows nothing when nothing is active', () => {
    expect(followCuts(new Map(), cuts(object('a')), 1).size).toBe(0);
  });
});
