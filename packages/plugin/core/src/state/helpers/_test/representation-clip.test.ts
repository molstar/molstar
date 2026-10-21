/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { RepresentationClip } from '../representation-clip.js';

const object = (type: 'plane' | 'sphere', x: number): RepresentationClip.Object => ({
  type,
  invert: false,
  position: Vec3.create(x, 0, 0),
  rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
  scale: Vec3.create(1, 1, 1),
  transform: Mat4.identity(),
});
const clip = (variant: Clip.Props['variant'], ...objects: RepresentationClip.Object[]): Clip.Props => ({
  variant,
  objects,
});

describe('RepresentationClip', () => {
  it('reads and writes the clip at the type params (structure, volume, particles) or at the top (shape)', () => {
    const typed = { type: { name: 'cartoon', params: { clip: clip('pixel') } }, colorTheme: {} };
    const shape = { clip: clip('instance'), alpha: 1 };
    const next = clip('pixel', object('plane', 1));
    expect(RepresentationClip.get(RepresentationClip.set(typed, next))).toBe(next);
    expect(RepresentationClip.set(typed, next).colorTheme).toBe(typed.colorTheme);
    expect(RepresentationClip.get(RepresentationClip.set(shape, next))).toBe(next);
    expect(RepresentationClip.get({ type: { name: 'x', params: {} } })).toBeUndefined();
  });

  it("replaces an owner's objects, keeping the others and the cell's variant while it has its own", () => {
    const own = object('sphere', 5);
    const globalA = object('plane', 1);
    const globalB = object('plane', 2);
    const cell = clip('instance', own, globalA);
    const next = RepresentationClip.replaceOwned(cell, [globalA], clip('pixel', globalB));
    expect(next.objects.map((o) => o.position[0])).toEqual([5, 2]);
    expect(next.variant).toBe('instance');
    // in place: an object's index is what a clipping layer's group means
    const before = clip('instance', globalA, own);
    expect(
      RepresentationClip.replaceOwned(before, [globalA], clip('pixel', globalB)).objects.map((o) => o.position[0]),
    ).toEqual([2, 5]);
    // more after, fewer drop their places
    const more = RepresentationClip.replaceOwned(before, [globalA], clip('pixel', globalB, object('plane', 3)));
    expect(more.objects.map((o) => o.position[0])).toEqual([2, 5, 3]);
    expect(RepresentationClip.replaceOwned(before, [globalA], clip('pixel')).objects.map((o) => o.position[0])).toEqual(
      [5],
    );
    // nothing of its own: the owner's variant
    expect(RepresentationClip.replaceOwned(clip('instance', globalA), [globalA], clip('pixel')).variant).toBe('pixel');
  });
});
