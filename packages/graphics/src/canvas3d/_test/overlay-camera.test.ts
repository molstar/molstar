/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Camera } from '../camera.js';
import { cameraProject } from '../camera/util.js';
import { getScreenSizeScale, OverlayCamera, type OverlayObject } from '../helper/overlay-helper.js';

/** A camera at (0, 0, `distance`) looking at the origin, its planes fitting a scene of `radius` there */
function createCamera(mode: Camera.Mode, radius = 5, distance = 100) {
  const camera = new Camera(
    { position: Vec3.create(0, 0, distance), mode, radius, radiusMax: radius, clipFar: true, minNear: 1 },
    { x: 0, y: 0, width: 200, height: 100 },
  );
  camera.update();
  return camera;
}

const cell = <T>(value: T) => ({ ref: { value } });

/** An overlay object of world bounds `sphere` */
function sphereObject(sphere: Sphere3D) {
  return { state: { visible: true }, values: { boundingSphere: cell(sphere) } } as unknown as OverlayObject;
}

/** A screen-sized overlay object (`uScreenSize`), its mesh in a unit sphere, one instance at `origin` */
function screenSizedObject(screenSize: number, origin: Vec3) {
  return {
    state: { visible: true },
    values: {
      uScreenSize: cell(screenSize),
      boundingSphere: cell(Sphere3D.create(origin, 1)),
      invariantBoundingSphere: cell(Sphere3D.create(Vec3(), 1)),
      aTransform: cell(Mat4.fromTranslation(Mat4(), origin) as unknown as Float32Array),
      instanceCount: cell(1),
    },
  } as unknown as OverlayObject;
}

/** Normalized device coordinates of `p` (world) */
function project(projectionView: ArrayLike<number>, p: Vec3) {
  const c = Vec4.transformMat4(Vec4(), Vec4.create(p[0], p[1], p[2], 1), projectionView as any);
  return Vec3.create(c[0] / c[3], c[1] / c[3], c[2] / c[3]);
}

describe('overlay camera', () => {
  for (const mode of ['perspective', 'orthographic'] as const) {
    it(`fits its planes to the overlay's bounds, the frustum otherwise the camera's (${mode})`, () => {
      const camera = createCamera(mode);
      // a gizmo well behind the scene: beyond the camera's far plane
      const sphere = Sphere3D.create(Vec3.create(0, 0, -50), 10);
      const overlay = new OverlayCamera().update(camera, [sphereObject(sphere)], 1);

      expect(camera.far).toBeLessThan(140);
      expect(overlay.near).toBeCloseTo(140 - 0.1, 1);
      expect(overlay.far).toBeCloseTo(160 + 0.1, 1);

      // in the overlay's depth range, out of the camera's
      const p = Vec3.create(3, 2, -55);
      const ndc = project(overlay.projectionView, p);
      expect(Math.abs(ndc[2])).toBeLessThan(1);
      expect(Math.abs(project(camera.projectionView, p)[2])).toBeGreaterThan(1);
      // the same place on screen
      const ndcCamera = project(camera.projectionView, p);
      expect(ndc[0]).toBeCloseTo(ndcCamera[0], 6);
      expect(ndc[1]).toBeCloseTo(ndcCamera[1], 6);
    });
  }

  it('has no fog', () => {
    const camera = createCamera('perspective');
    camera.setState({ fog: 50 });
    camera.update();
    const overlay = new OverlayCamera().update(camera, [sphereObject(Sphere3D.create(Vec3.create(0, 0, 0), 1))], 1);
    expect(overlay.state.fog).toBe(0);
  });

  it("keeps the camera's planes when it has nothing to fit (empty, all behind the camera)", () => {
    const camera = createCamera('perspective');
    const empty = new OverlayCamera().update(camera, [sphereObject(Sphere3D.create(Vec3.create(0, 0, 0), 0))], 1);
    expect([empty.near, empty.far]).toEqual([camera.near, camera.far]);
    const behind = new OverlayCamera().update(camera, [sphereObject(Sphere3D.create(Vec3.create(0, 0, 200), 10))], 1);
    expect([behind.near, behind.far]).toEqual([camera.near, camera.far]);
  });

  it('keeps the near plane of a perspective in front of the camera when inside the bounds', () => {
    const camera = createCamera('perspective');
    const overlay = new OverlayCamera().update(camera, [sphereObject(Sphere3D.create(Vec3.create(0, 0, 100), 20))], 1);
    expect(overlay.near).toBeGreaterThan(0);
    expect(overlay.far).toBeCloseTo(20.2, 1);
  });

  it('follows the camera scale (the view is in scaled units)', () => {
    const camera = createCamera('perspective');
    camera.scale = 2;
    camera.update();
    const overlay = new OverlayCamera().update(camera, [sphereObject(Sphere3D.create(Vec3.create(0, 0, 0), 10))], 1);
    // the camera at 200 in scaled units, the sphere's radius 20
    expect(overlay.near).toBeCloseTo(200 - 20.2, 1);
    expect(overlay.far).toBeCloseTo(200 + 20.2, 1);
  });
});

describe('screen size', () => {
  for (const mode of ['perspective', 'orthographic'] as const) {
    it(`is the same in CSS pixels at two zooms and pixel ratios (${mode})`, () => {
      const origin = Vec3.create(3, -2, 0);
      for (const distance of [40, 400]) {
        for (const pixelRatio of [1, 2]) {
          const camera = createCamera(mode, 5, distance);
          const s = getScreenSizeScale(camera, origin, 50, pixelRatio);
          // one unit of the mesh along the view's x axis, in window coordinates (the viewport's pixels)
          const a = cameraProject(Vec4(), origin, camera.viewport, camera.projectionView);
          const b = cameraProject(
            Vec4(),
            Vec3.create(origin[0] + s, origin[1], 0),
            camera.viewport,
            camera.projectionView,
          );
          expect((b[0] - a[0]) / pixelRatio).toBeCloseTo(50, 3);
        }
      }
    });
  }

  for (const mode of ['perspective', 'orthographic'] as const) {
    it(`keeps a screen-sized mesh's size where the viewport is cropped (${mode})`, () => {
      const camera = createCamera(mode);
      const size = getScreenSizeScale(camera, Vec3(), 50, 2);
      Camera.setCrop(camera.crop, 200, 100, 50, 25);
      Object.assign(camera.viewport, { width: 100, height: 50 });
      camera.update();
      expect(getScreenSizeScale(camera, Vec3(), 50, 2)).toBeCloseTo(size, 6);
    });
  }

  it('fits the planes to a screen-sized object at the size the camera gives it', () => {
    // its mesh in a unit sphere: in the world far larger than that when zoomed out
    for (const distance of [40, 4000]) {
      const camera = createCamera('perspective', 5, distance);
      const object = screenSizedObject(50, Vec3());
      const radius = getScreenSizeScale(camera, Vec3(), 50, 1);
      const overlay = new OverlayCamera().update(camera, [object], 1);
      expect(overlay.near).toBeCloseTo(distance - radius * 1.01, 3);
      expect(overlay.far).toBeCloseTo(distance + radius * 1.01, 3);
    }
  });
});
