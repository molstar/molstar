/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { CullView } from '../../gl/camera.js';
import { Camera } from '../camera.js';
import { StereoCamera } from '../camera/stereo.js';

export const JitterVectors = [
  [[0, 0]],
  [
    [0, 0],
    [-4, -4],
  ],
  [
    [0, 0],
    [6, -2],
    [-6, 2],
    [2, 6],
  ],
  [
    [0, 0],
    [-1, 3],
    [5, 1],
    [-3, -5],
    [-5, 5],
    [-7, -1],
    [3, 7],
    [7, -7],
  ],
  [
    [0, 0],
    [-1, -3],
    [-3, 2],
    [4, -1],
    [-5, -2],
    [2, 5],
    [5, 3],
    [3, -5],
    [-2, 6],
    [0, -7],
    [-4, -6],
    [-6, 4],
    [-8, 0],
    [7, -4],
    [6, 7],
    [-7, -8],
  ],
  [
    [0, 0],
    [-7, -5],
    [-3, -5],
    [-5, -4],
    [-1, -4],
    [-2, -2],
    [-6, -1],
    [-4, 0],
    [-7, 1],
    [-1, 2],
    [-6, 3],
    [-3, 3],
    [-7, 6],
    [-3, 6],
    [-5, 7],
    [-1, 7],
    [5, -7],
    [1, -6],
    [6, -5],
    [4, -4],
    [2, -3],
    [7, -2],
    [1, -1],
    [4, -1],
    [2, 1],
    [6, 2],
    [0, 4],
    [4, 4],
    [2, 5],
    [7, 5],
    [5, 6],
    [3, 7],
  ],
];

JitterVectors.forEach((offsetList) => {
  offsetList.forEach((offset) => {
    // 0.0625 = 1 / 16
    offset[0] *= 0.0625;
    offset[1] *= 0.0625;
  });
});

export function getJitterOffsets(sampleLevel: number) {
  return JitterVectors[Math.max(0, Math.min(sampleLevel, 5))];
}

export function getTemporalSamplesPerFrame(sampleLevel: number) {
  return Math.pow(2, Math.max(0, sampleLevel - 2));
}

const cullViews = new WeakMap<Camera, CullView>();
const tmpInverse = Mat4();
const tmpPosition = Vec3();

/**
 * The culling of a camera's jittered samples (`RenderCamera.cullView`): the camera's unjittered, as the renderer makes
 * it (`Renderer.update`), the same to the bit. A resting camera's frame is kept over its samples (idle ticks reuse its
 * culling, `Canvas3D`), culled once unjittered: the samples draw that culling's lists (a merged object's GPU culling,
 * one dispatch a frame), not a CPU culling of each of their own views (a merged object's first builds its members'
 * grids: 1.7 s for the yeast cell). A jitter is half a pixel at most: an object no nearer the screen's edge left out
 */
function setJitterCullView(camera: Camera) {
  let view = cullViews.get(camera);
  if (!view) {
    view = { frustum: Frustum3D(), cameraPlane: Plane3D() };
    cullViews.set(camera, view);
  }
  const { frustum, cameraPlane } = view;
  Frustum3D.fromProjectionMatrix(frustum, camera.projectionView);
  Plane3D.copy(cameraPlane, frustum[Frustum3D.PlaneIndex.Near]);
  const position = Mat4.getTranslation(tmpPosition, Mat4.invert(tmpInverse, camera.view));
  cameraPlane.constant -= Plane3D.distanceToPoint(cameraPlane, position);
  camera.cullView = view;
}

export function setJitter(camera: Camera | StereoCamera, offset: number[]) {
  const { width, height } = camera.viewport;
  // the samples culled as one, from the camera unjittered (stereo's eyes cull as one of their own)
  if (camera instanceof Camera && !camera.viewOffset.enabled) setJitterCullView(camera);
  camera.viewOffset.enabled = true;
  Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height);
  updateSilently(camera);
}

/** a sample's jitter no move of the camera: its `changed` not told (a scene changed: a frame's samples, two each) */
function updateSilently(camera: Camera | StereoCamera) {
  if (camera instanceof Camera) camera.update(true);
  else camera.update(undefined, true);
}

export function clearJitter(camera: Camera | StereoCamera) {
  camera.viewOffset.enabled = false;
  if (camera instanceof Camera) camera.cullView = undefined;
  updateSilently(camera);
}
