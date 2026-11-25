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
