/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import type { ParticleTileCamera } from '@molstar/model/model/particles/particle-tile-visibility';

/** the far plane of a tile camera: deep, as tiles not loaded are no part of the scene the canvas' far plane bounds */
const TileCameraFar = 1e8;

/** a canvas camera's snapshot as a tile set's camera reads it: its clipping's fields where known (`getClipDistances`) */
export type ParticleTileCameraSnapshot = Pick<
  Camera.Snapshot,
  'mode' | 'fov' | 'position' | 'target' | 'up' | 'minNear'
> &
  Partial<Pick<Camera.Snapshot, 'radius' | 'radiusMax' | 'clipFar' | 'minFar'>>;

/**
 * A tile set's camera (`ParticleTileVisibility`) of a canvas camera's snapshot and a viewport of `width` × `height`
 * pixels: its view, its projection (perspective or orthographic, as the snapshot's mode; its near plane the snapshot's
 * `minNear`, its far plane deep: the canvas' own is the loaded scene's, which tiles not loaded are not of). What hides
 * starts at the canvas camera's own near plane (`Camera.getClipDistances`: nearer than its focus' radius nothing is
 * drawn, the clipping slider's cut), not the projection's: a tile before it still in view (loaded, a reset framing
 * it), hiding nothing
 */
export function getParticleTileCamera(
  snapshot: ParticleTileCameraSnapshot,
  width: number,
  height: number,
): ParticleTileCamera {
  const { position, target, up, fov } = snapshot;
  const near = Math.max(0.01, snapshot.minNear);
  const { radius, radiusMax, clipFar, minFar } = snapshot;
  const clipNear =
    radius !== undefined && radiusMax !== undefined
      ? Camera.getClipDistances({
          position,
          target,
          radius,
          radiusMax,
          clipFar: clipFar ?? true,
          minNear: near,
          minFar: minFar ?? 0,
        }).near
      : near;
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
    near: Math.max(near, clipNear),
  };
}
