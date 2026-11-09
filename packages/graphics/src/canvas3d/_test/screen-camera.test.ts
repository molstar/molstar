/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Representation } from '@molstar/graphics/repr/representation';
import { Camera } from '../camera.js';
import { cameraProject } from '../camera/util.js';
import { ScreenCamera, ScreenHelper } from '../helper/screen-helper.js';

/** A camera of a 400 × 200 device pixel viewport, looking at the origin from an oblique direction */
function createCamera() {
  const camera = new Camera(
    { position: Vec3.create(30, 40, 100), radius: 5, radiusMax: 5 },
    { x: 0, y: 0, width: 400, height: 200 },
  );
  camera.update();
  return camera;
}

/** Window coordinates (device pixels, top-left origin) of `p` in the screen camera's view */
function project(camera: ScreenCamera, p: Vec3) {
  const c = cameraProject(Vec4(), p, camera.viewport, camera.projectionView);
  return [c[0], c[1]];
}

const placement = (p: Partial<Representation.ScreenPlacement>): Representation.ScreenPlacement => ({
  anchor: 'center',
  offsetX: 0,
  offsetY: 0,
  followRotation: true,
  ...p,
});

describe('screen camera', () => {
  // inset 10 and offset (5, 7) CSS pixels at a pixel ratio of 2: 30 and 34 device pixels from the edges anchored to
  const corners = [
    ['bottom-left', [30, 200 - 34]],
    ['bottom-right', [400 - 30, 200 - 34]],
    ['top-left', [30, 34]],
    ['top-right', [400 - 30, 34]],
    ['center', [200 + 10, 100 - 14]],
    // an edge's middle: inward from that edge, the offset as is along it (right and up)
    ['top-middle', [200 + 10, 34]],
    ['bottom-middle', [200 + 10, 200 - 34]],
    ['middle-left', [30, 100 - 14]],
    ['middle-right', [400 - 30, 100 - 14]],
  ] as const;
  for (const [anchor, expected] of corners) {
    it(`places an object's origin at its anchor (${anchor})`, () => {
      const camera = new ScreenCamera().update(
        createCamera(),
        placement({ anchor, offsetX: 5, offsetY: 7 }),
        20,
        10,
        2,
      );
      const [x, y] = project(camera, Vec3());
      expect(x).toBeCloseTo(expected[0], 3);
      expect(y).toBeCloseTo(expected[1], 3);
    });
  }

  it('has CSS pixels as its unit', () => {
    const camera = new ScreenCamera().update(createCamera(), placement({ followRotation: false }), 20, 10, 2);
    const a = project(camera, Vec3());
    const b = project(camera, Vec3.create(10, 0, 0));
    expect(b[0] - a[0]).toBeCloseTo(20, 3);
  });

  it('turns an object with the camera, or keeps it upright', () => {
    const source = createCamera();
    const rotation = Mat4.extractRotation(Mat4(), source.view);
    const turned = new ScreenCamera().update(source, placement({}), 20, 10, 1);
    const upright = new ScreenCamera().update(source, placement({ followRotation: false }), 20, 10, 1);
    const p = Vec3.create(10, 0, 0);
    // as the camera's rotation turns it
    const q = Vec3.transformMat4(Vec3(), p, rotation);
    const [x0, y0] = project(upright, Vec3());
    const [x1, y1] = project(turned, p);
    expect(x1 - x0).toBeCloseTo(q[0], 3);
    expect(y0 - y1).toBeCloseTo(q[1], 3);
    // upright: along the screen's x
    const [x2, y2] = project(upright, p);
    expect(x2 - x0).toBeCloseTo(10, 3);
    expect(y2 - y0).toBeCloseTo(0, 3);
  });

  it('stays at its place in the frame where the viewport is a cropped part of it', () => {
    const full = createCamera();
    const fullCamera = new ScreenCamera().update(full, placement({ anchor: 'bottom-left' }), 20, 10, 2);
    const [fx, fy] = project(fullCamera, Vec3());
    // the bottom-right quarter of the frame drawn: the anchor left of the viewport, as far as it is cropped
    const cropped = createCamera();
    Camera.setCrop(cropped.crop, 400, 200, 200, 100);
    Object.assign(cropped.viewport, { width: 200, height: 100 });
    const croppedCamera = new ScreenCamera().update(cropped, placement({ anchor: 'bottom-left' }), 20, 10, 2);
    const [cx, cy] = project(croppedCamera, Vec3());
    expect(cx).toBeCloseTo(fx - 200, 3);
    expect(cy).toBeCloseTo(fy - 100, 3);
  });

  it("takes the jitter's view offset as sub-pixels, in device pixels", () => {
    const source = createCamera();
    const still = project(new ScreenCamera().update(source, placement({}), 20, 10, 2), Vec3());
    source.viewOffset.enabled = true;
    Camera.setViewOffset(source.viewOffset, 400, 200, 0.5, 0.25, 400, 200);
    const [x, y] = project(new ScreenCamera().update(source, placement({}), 20, 10, 2), Vec3());
    expect(x).toBeCloseTo(still[0] - 0.5, 3);
    expect(y).toBeCloseTo(still[1] - 0.25, 3);
  });

  it('holds the object between its near and far planes', () => {
    const camera = new ScreenCamera().update(createCamera(), placement({}), 20, 10, 1);
    for (const z of [-20, 0, 20]) {
      const c = Vec4.transformMat4(Vec4(), Vec4.create(0, 0, z, 1), camera.projectionView);
      expect(Math.abs(c[2] / c[3])).toBeLessThan(1);
    }
  });
});

describe('screen layer', () => {
  it('is off in XR (no screen to anchor to)', () => {
    const xr = { session: undefined as unknown };
    const graphics = { xr } as unknown as GraphicsContext;
    const screen = new ScreenHelper(graphics);
    const owner = {};
    const scene = screen.getScene(owner, placement({}));
    Object.defineProperty(scene, 'count', { value: 1 });
    expect(screen.isEnabled).toBe(true);
    xr.session = {};
    expect(screen.isEnabled).toBe(false);
    xr.session = undefined;
    expect(screen.isEnabled).toBe(true);
    expect(screen.has(scene)).toBe(true);
  });
});

describe('ScreenHelper.scenes', () => {
  it('is the same list until a scene is added or released', () => {
    const screen = new ScreenHelper({} as unknown as GraphicsContext);
    const placement = { corner: 'bottom-left', offset: [0, 0] } as unknown as Representation.ScreenPlacement;
    expect(screen.scenes).toHaveLength(0);
    const a = {};
    const sa = screen.getScene(a, placement);
    const first = screen.scenes;
    expect(first).toEqual([sa]);
    // read again (every tick): the same list
    expect(screen.scenes).toBe(first);
    expect(screen.getScene(a, placement)).toBe(sa);
    expect(screen.scenes).toBe(first);
    const sb = screen.getScene({}, placement);
    expect(screen.scenes).toEqual([sa, sb]);
    screen.release(a);
    expect(screen.scenes).toEqual([sb]);
    screen.clear();
    expect(screen.scenes).toHaveLength(0);
  });
});
