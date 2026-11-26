/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Camera } from '../camera.js';
import { Viewport, cameraProject, cameraUnproject } from '../camera/util.js';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { clearJitter, setJitter } from '../passes/jitter.js';

describe('camera', () => {
  it('project/unproject', () => {
    const proj = Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 100);
    const invProj = Mat4.invert(Mat4(), proj);

    const c = Vec4();
    const po = Vec3();

    const vp = Viewport.create(0, 0, 100, 100);
    const pi = Vec3.create(0, 0, 1);
    cameraProject(c, pi, vp, proj);
    expect(Vec4.equals(c, Vec4.create(50, 50, 2.020202, -1))).toBe(true);
    cameraUnproject(po, c, vp, invProj);
    expect(Vec3.equals(po, pi)).toBe(true);

    Vec3.set(pi, 0.5, 0.5, 1);
    cameraProject(c, pi, vp, proj);
    cameraUnproject(po, c, vp, invProj);
    expect(Vec3.equals(po, pi)).toBe(true);

    Viewport.set(vp, 50, 50, 100, 100);
    Vec3.set(pi, 0.5, 0.5, 1);
    cameraProject(c, pi, vp, proj);
    cameraUnproject(po, c, vp, invProj);
    expect(Vec3.equals(po, pi)).toBe(true);
  });

  it('window coordinates have a top-left origin', () => {
    const proj = Mat4.ortho(Mat4(), -1, 1, 1, -1, 1, 100);
    const vp = Viewport.create(10, 20, 100, 100);
    const c = Vec4();
    // up in view space is up on screen, i.e. a smaller window y
    cameraProject(c, Vec3.create(0, 0.5, -10), vp, proj);
    expect(c[0]).toBeCloseTo(60);
    expect(c[1]).toBeCloseTo(45);
  });

  describe('crop', () => {
    // a 400 × 200 frame, of which the 150 × 120 part at (100, 30) is the viewport, placed at (20, 10) on the canvas
    const frame = Viewport.create(20 - 100, 10 - 30, 400, 200);
    const part = Viewport.create(20, 10, 150, 120);
    const points = [Vec3.create(0, 0, 0), Vec3.create(3, -2, 1), Vec3.create(-4, 1.5, -2), Vec3.create(1, 2.5, 3)];

    function cameras(mode: Camera.Mode, jitter?: number[]) {
      const state = { mode, position: Vec3.create(2, 3, 20), target: Vec3.create(0, 0, 0), radiusMax: 20 };
      const full = new Camera(state, Viewport.clone(frame));
      const cropped = new Camera(state, Viewport.clone(part));
      Camera.setCrop(cropped.crop, frame.width, frame.height, part.x - frame.x, part.y - frame.y);
      for (const c of [full, cropped]) {
        if (jitter) setJitter(c, jitter);
        else c.update();
      }
      return { full, cropped };
    }

    for (const mode of ['perspective', 'orthographic'] as const) {
      for (const jitter of [undefined, [0.25, -0.375]]) {
        it(`projects a point of the frame onto the same pixel (${mode}${jitter ? ', jittered' : ''})`, () => {
          const { full, cropped } = cameras(mode, jitter);
          const a = Vec4();
          const b = Vec4();
          for (const p of points) {
            full.project(a, p);
            cropped.project(b, p);
            expect(b[0]).toBeCloseTo(a[0], 6);
            expect(b[1]).toBeCloseTo(a[1], 6);
            expect(b[2]).toBeCloseTo(a[2], 6);
          }
        });
      }

      it(`measures the view by the frame (${mode})`, () => {
        const { full, cropped } = cameras(mode);
        expect(Camera.getFrame(Viewport.zero(), cropped)).toEqual(frame);
        expect(cropped.zoom).toBeCloseTo(full.zoom, 9);
        expect(cropped.getTargetDistance(5)).toBeCloseTo(full.getTargetDistance(5), 9);
        expect(cropped.getPixelSize(points[1])).toBeCloseTo(full.getPixelSize(points[1]), 9);
      });
    }
  });

  describe('clipping', () => {
    // the camera 100 from its target, the scene's radius 40 (its front at 60, its back at 140), clipped to 10 around the target
    function planes(clip: Partial<Camera.Snapshot>) {
      const camera = new Camera(
        {
          position: Vec3.create(0, 0, 100),
          target: Vec3.create(0, 0, 0),
          radius: 10,
          radiusMax: 40,
          minNear: 1,
          ...clip,
        },
        Viewport.create(0, 0, 100, 100),
      );
      camera.update();
      return [camera.near, camera.far];
    }

    it('clips both sides by the radius', () => {
      expect(planes({})).toEqual([90, 110]);
    });

    it('keeps the front with near off', () => {
      expect(planes({ clipNear: false })).toEqual([60, 110]);
    });

    it('keeps the back with far off', () => {
      expect(planes({ clipFar: false })).toEqual([90, 140]);
    });

    it('keeps both sides with near and far off', () => {
      expect(planes({ clipNear: false, clipFar: false })).toEqual([60, 140]);
    });

    it('copies and compares clipNear', () => {
      const a = Camera.createDefaultSnapshot();
      const b = Camera.copySnapshot(Camera.createDefaultSnapshot(), { clipNear: false });
      expect(b.clipNear).toBe(false);
      expect(Camera.areSnapshotsEqual(a, b)).toBe(false);
      expect(Camera.areSnapshotsEqual(a, Camera.copySnapshot(Camera.createDefaultSnapshot(), a))).toBe(true);
    });
  });

  it("tells no change for a sample's jitter, the projection jittered and restored; a move told", () => {
    const camera = new Camera(
      { position: Vec3.create(2, 3, 20), target: Vec3.create(0, 0, 0), radiusMax: 20 },
      Viewport.create(0, 0, 300, 200),
    );
    camera.update();
    const unjittered = Mat4.copy(Mat4(), camera.projection);
    let changes = 0;
    const sub = camera.changed.subscribe(() => ++changes);
    setJitter(camera, [0.25, -0.375]);
    expect(Mat4.areEqual(camera.projection, unjittered, 0)).toBe(false);
    clearJitter(camera);
    expect(Array.from(camera.projection)).toEqual(Array.from(unjittered));
    // the frame's update after the samples: nothing changed since the last told
    expect(camera.update()).toBe(false);
    expect(changes).toBe(0);
    camera.position = Vec3.create(2, 3, 25);
    expect(camera.update()).toBe(true);
    expect(changes).toBe(1);
    sub.unsubscribe();
  });

  it('culls jittered samples as the camera unjittered, to the bit, until the jitter is cleared', () => {
    const camera = new Camera(
      { position: Vec3.create(2, 3, 20), target: Vec3.create(0, 0, 0), radiusMax: 20 },
      Viewport.create(0, 0, 300, 200),
    );
    camera.update();
    // as the renderer makes them (`Renderer.update`)
    const frustum = Frustum3D.fromProjectionMatrix(Frustum3D(), camera.projectionView);
    const cameraPlane = Plane3D.copy(Plane3D(), frustum[Frustum3D.PlaneIndex.Near]);
    cameraPlane.constant -= Plane3D.distanceToPoint(
      cameraPlane,
      Mat4.getTranslation(Vec3(), Mat4.invert(Mat4(), camera.view)),
    );
    const unjittered = Mat4.copy(Mat4(), camera.projectionView);
    expect(camera.cullView).toBeUndefined();
    for (const offset of [
      [0, 0],
      [0.25, -0.375],
      [-0.5, 0.125],
    ]) {
      setJitter(camera, offset);
      const view = camera.cullView!;
      expect(view).toBeDefined();
      for (let i = 0; i < 6; ++i) {
        expect(Array.from(view.frustum[i as Frustum3D.PlaneIndex].normal)).toEqual(
          Array.from(frustum[i as Frustum3D.PlaneIndex].normal),
        );
        expect(view.frustum[i as Frustum3D.PlaneIndex].constant).toBe(frustum[i as Frustum3D.PlaneIndex].constant);
      }
      expect(Array.from(view.cameraPlane.normal)).toEqual(Array.from(cameraPlane.normal));
      expect(view.cameraPlane.constant).toBe(cameraPlane.constant);
    }
    // the samples drawn jittered all the same
    expect(Mat4.areEqual(camera.projectionView, unjittered, 1e-9)).toBe(false);
    clearJitter(camera);
    expect(camera.cullView).toBeUndefined();
  });
});
