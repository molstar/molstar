/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
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
    // more after; fewer: a place before the cell's own held by a `none` (its index kept), a place at the end gone
    const more = RepresentationClip.replaceOwned(before, [globalA], clip('pixel', globalB, object('plane', 3)));
    expect(more.objects.map((o) => o.position[0])).toEqual([2, 5, 3]);
    const fewer = RepresentationClip.replaceOwned(before, [globalA], clip('pixel'));
    expect(fewer.objects.map((o) => o.type)).toEqual(['none', 'sphere']);
    expect(fewer.objects[1].position[0]).toBe(5);
    expect(RepresentationClip.replaceOwned(cell, [globalA], clip('pixel')).objects.map((o) => o.position[0])).toEqual([
      5,
    ]);
    // a held place is the owner's to fill again
    const again = RepresentationClip.replaceOwned(fewer, [], clip('pixel', globalB));
    expect(again.objects.map((o) => o.position[0])).toEqual([2, 5]);
    // nothing of its own: the owner's variant
    expect(RepresentationClip.replaceOwned(clip('instance', globalA), [globalA], clip('pixel')).variant).toBe('pixel');
  });

  it('is the same cut by its tag when tagged, by its values when not', () => {
    const a = object('plane', 1);
    const moved = object('plane', 2);
    expect(RepresentationClip.sameCut(a, a)).toBe(true);
    expect(RepresentationClip.sameCut(a, moved)).toBe(false);
    // tagged: the tag, wherever it is
    expect(RepresentationClip.sameCut({ ...a, tag: 'clip-1' }, { ...moved, tag: 'clip-1' })).toBe(true);
    expect(RepresentationClip.sameCut({ ...a, tag: 'clip-1' }, { ...a, tag: 'clip-2' })).toBe(false);
    expect(RepresentationClip.sameCut({ ...a, tag: 'clip-1' }, a)).toBe(false);
    // equal in every value takes the tag too
    expect(RepresentationClip.objectsEqual({ ...a, tag: 'clip-1' }, a)).toBe(false);
  });
});
