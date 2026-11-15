/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import type { ParticleTileCamera } from '@molstar/model/model/particles/particle-tile-visibility';

/** the far plane of a tile camera: deep, as tiles not loaded are no part of the scene the canvas' far plane bounds */
const TileCameraFar = 1e8;

/**
 * A tile set's camera (`ParticleTileVisibility`) of a canvas camera's snapshot and a viewport of `width` × `height`
 * pixels: its view, its projection (perspective or orthographic, as the snapshot's mode; its near plane the snapshot's
 * `minNear`, its far plane deep: the canvas' own is the loaded scene's, which tiles not loaded are not of)
 */
export function getParticleTileCamera(
  snapshot: Pick<Camera.Snapshot, 'mode' | 'fov' | 'position' | 'target' | 'up' | 'minNear'>,
  width: number,
  height: number,
): ParticleTileCamera {
  const { position, target, up, fov } = snapshot;
  const near = Math.max(0.01, snapshot.minNear);
  const view = Mat4.lookAt(Mat4(), position, target, up);
  const aspect = width / Math.max(1, height);
  const projection = Mat4();
  const orthographic = snapshot.mode === 'orthographic';
  // the half height the projection maps to the viewport's half: at a distance of 1 (perspective), at the target
  const half = orthographic ? Math.tan(fov / 2) * Vec3.distance(position, target) : Math.tan(fov / 2);
  if (orthographic) {
    Mat4.ortho(projection, -half * aspect, half * aspect, half, -half, near, TileCameraFar);
  } else {
    const top = near * half;
    Mat4.perspective(projection, -top * aspect, top * aspect, top, -top, near, TileCameraFar);
  }
  const direction = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), target, position));
  return {
    position: Vec3.clone(position),
    direction,
    viewProjection: Mat4.mul(Mat4(), projection, view),
    orthographic,
    pixelScale: height / 2 / half,
    near,
  };
}
