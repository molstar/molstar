/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Camera } from '../camera.js';
import { Viewport, cameraProject, cameraUnproject } from '../camera/util.js';
import { setJitter } from '../passes/jitter.js';

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
});
