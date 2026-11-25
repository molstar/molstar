/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Camera } from '../camera.js';
import { StereoCamera } from '../camera/stereo.js';
import { Viewport } from '../camera/util.js';
import { getStereoSplit, isStereoPyramidOf } from '../passes/hi-z.js';

function stereo(viewport = Viewport.create(0, 0, 400, 200)) {
  const camera = new Camera(
    { position: Vec3.create(0, 0, 100), target: Vec3.create(0, 0, 0), radius: 30, radiusMax: 50 },
    viewport,
  );
  camera.update();
  const s = new StereoCamera(camera, { eyeSeparation: 0.062, focus: 10 });
  s.update();
  return { camera, s };
}

/** a pyramid built of `s`'s eyes as they are now */
function pyramidOf(s: StereoCamera, split: number): OcclusionPyramid {
  return {
    texture: undefined as unknown as Texture,
    levelCount: 1,
    minLevel: 0,
    width: 400,
    height: 200,
    view: Mat4.clone(s.left.view),
    projection: Mat4.clone(s.left.projection),
    near: s.left.near,
    far: s.left.far,
    second: { view: Mat4.clone(s.right.view), projection: Mat4.clone(s.right.projection), split },
    version: 1,
  };
}

describe('stereo occlusion levels', () => {
  it("place the right eye's levels at the middle of the viewport, side by side", () => {
    const { camera, s } = stereo(Viewport.create(10, 20, 400, 200));
    expect(getStereoSplit(s, camera.viewport)).toBe(0.5);
  });

  it('build none where an eye is disabled or the eyes do not tile the viewport', () => {
    const { camera, s } = stereo();
    s.right.disabled = true;
    expect(getStereoSplit(s, camera.viewport)).toBeUndefined();
    s.right.disabled = false;
    expect(getStereoSplit(s, Viewport.create(0, 0, 400, 100))).toBeUndefined();
  });

  it('are of the eyes until the camera moves', () => {
    const { camera, s } = stereo();
    const pyramid = pyramidOf(s, 0.5);
    s.update();
    expect(isStereoPyramidOf(pyramid, s)).toBe(true);
    camera.setState({ position: Vec3.create(10, 0, 100) }, 0);
    camera.update();
    s.update();
    expect(isStereoPyramidOf(pyramid, s)).toBe(false);
    expect(isStereoPyramidOf({ ...pyramidOf(s, 0.5), second: undefined }, s)).toBe(false);
  });
});
