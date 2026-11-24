/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { BaseGeometry } from '../base.js';
import { Spheres } from '../spheres/spheres.js';
import { ColorNames } from '@molstar/core/util/color/names';

function clipAt(x: number): Clip.Props {
  return {
    variant: 'pixel',
    objects: [
      {
        type: 'plane',
        invert: false,
        position: Vec3.create(x, 0, 0),
        rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
        scale: Vec3.create(1, 1, 1),
        transform: Mat4.identity(),
      },
    ],
  };
}

describe('BaseGeometry.isClipChangeOnly', () => {
  const before = { ...PD.getDefaultValues(BaseGeometry.Params), clip: clipAt(0) };

  it('is a clip moved alone, the other props the same values', () => {
    expect(BaseGeometry.isClipChangeOnly({ ...before, clip: clipAt(1) }, before)).toBe(true);
  });

  it('is not the same props, the same clip, nor another prop changed with it', () => {
    expect(BaseGeometry.isClipChangeOnly(before, before)).toBe(false);
    expect(BaseGeometry.isClipChangeOnly({ ...before }, before)).toBe(false);
    expect(BaseGeometry.isClipChangeOnly({ ...before, clip: clipAt(1), alpha: 0.5 }, before)).toBe(false);
    // a prop dropped
    const { alpha: _, ...dropped } = { ...before, clip: clipAt(1) };
    expect(BaseGeometry.isClipChangeOnly(dropped as typeof before, before)).toBe(false);
  });
});

describe('BaseGeometry.updateClip', () => {
  it('writes a clip again where values sharing its cells were written another (a merged part and its records)', () => {
    const spheres = Spheres.create(new Float32Array([0, 0, 0]), new Float32Array([0]), 1);
    const shared = Spheres.Utils.createValuesSimple(spheres, PD.getDefaultValues(Spheres.Params), ColorNames.orange, 1);
    // a record: the shared values' clip cells read through its prototype
    const record = Object.create(shared) as typeof shared;
    const a = clipAt(0);
    BaseGeometry.updateClip(record, a);
    const frameA = Array.from(shared.uClipObjectFrame.ref.value);
    BaseGeometry.updateClip(shared, clipAt(1));
    expect(Array.from(shared.uClipObjectFrame.ref.value)).not.toEqual(frameA);
    BaseGeometry.updateClip(record, a);
    expect(Array.from(shared.uClipObjectFrame.ref.value)).toEqual(frameA);
  });
});
