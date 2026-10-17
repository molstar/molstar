/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';
import { type CameraSnapshot, sameView, syncCameras } from '../camera.js';

// The camera as a view to follow and set (api.md §12), on stand-in canvases with a real `Camera` (no WebGL here): the
// drawn frame is `camera.update()`, which is what announces a change.

const open: PrismViewer[] = [];
afterEach(async () => {
  for (const v of open.splice(0)) await v.dispose();
});

async function viewer() {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas'] });
  open.push(v);
  return v;
}

/** A canvas with a real camera; `frame()` is the draw loop's `camera.update()`; `draws` counts redraw requests. */
function mount(v: PrismViewer, radiusMax = 100) {
  const camera = new Camera({
    position: Vec3.create(0, 0, 100),
    target: Vec3.create(0, 0, 0),
    up: Vec3.create(0, 1, 0),
    radius: Math.min(50, radiusMax),
    radiusMax,
  });
  const canvas = { camera, draws: 0, requestDraw: () => void canvas.draws++, dispose() {} };
  Object.defineProperty(v.plugin, 'canvas3d', { value: canvas, configurable: true });
  (v.plugin as unknown as { canvas3dInit: { next(v: boolean): void } }).canvas3dInit.next(true);
  camera.update();
  return { camera, canvas, frame: () => camera.update() };
}

const at = (x: number): CameraSnapshot => ({
  mode: 'perspective',
  fov: Math.PI / 4,
  position: [x, 0, 100],
  target: [x, 0, 0],
  up: [0, 1, 0],
  radius: 50,
});

describe('the viewer camera', () => {
  it('has no view without a canvas, and cannot be set', async () => {
    const v = await viewer();
    expect(v.camera.snapshot()).toBeUndefined();
    expect(v.camera.set(at(1))).toBe(false);
  });

  it('reads the view, and announces the frames that moved it', async () => {
    const v = await viewer();
    const m = mount(v);
    expect(v.camera.snapshot()).toMatchObject({ mode: 'perspective', position: [0, 0, 100], target: [0, 0, 0] });
    const seen: CameraSnapshot[] = [];
    v.camera.changed.subscribe((s) => seen.push(s));

    // The person drags: the controls change the state, the next frame announces it once.
    Vec3.set(m.camera.state.position, 5, 0, 100);
    Vec3.set(m.camera.state.target, 5, 0, 0);
    m.frame();
    expect(seen.map((s) => s.position)).toEqual([[5, 0, 100]]);
    // Frames that redraw the same view are no change.
    m.frame();
    m.camera.changed.next();
    expect(seen).toHaveLength(1);
  });

  it('sets a view as a cut, announced unless silent', async () => {
    const v = await viewer();
    const m = mount(v);
    const seen: CameraSnapshot[] = [];
    v.camera.changed.subscribe((s) => seen.push(s));

    expect(v.camera.set(at(3))).toBe(true);
    expect(m.canvas.draws).toBe(1);
    expect(m.camera.state.target[0]).toBe(3);
    m.frame();
    expect(seen).toHaveLength(1);

    v.camera.set(at(7), { silent: true });
    m.frame();
    expect(m.camera.state.target[0]).toBe(7);
    expect(seen).toHaveLength(1);
    // After a silent set, the next real move is heard.
    v.camera.set(at(9));
    m.frame();
    expect(seen.map((s) => s.target[0])).toEqual([3, 9]);
  });

  it('ends a transition under way when set, so its frames are not announced', async () => {
    const v = await viewer();
    const m = mount(v);
    const seen: CameraSnapshot[] = [];
    v.camera.changed.subscribe((s) => seen.push(s));
    m.camera.setState({ target: Vec3.create(20, 0, 0), position: Vec3.create(20, 0, 100) }, 1000);
    m.camera.transition.tick(100);
    expect(m.camera.transition.inTransition).toBe(true);
    v.camera.set(at(5), { silent: true });
    expect(m.camera.transition.inTransition).toBe(false);
    expect(m.camera.state.target[0]).toBe(5);
    m.camera.transition.tick(200);
    m.frame();
    expect(m.camera.state.target[0]).toBe(5);
    expect(seen).toHaveLength(0);
  });

  it('compares views by what the eye sees', () => {
    expect(sameView(at(1), at(1))).toBe(true);
    expect(sameView(at(1), { ...at(1), position: [1 + 1e-9, 0, 100] })).toBe(true);
    expect(sameView(at(1), at(1.001))).toBe(false);
    expect(sameView(at(1), { ...at(1), mode: 'orthographic' })).toBe(false);
    expect(sameView(at(1), { ...at(1), radius: 40 })).toBe(false);
  });
});

describe('syncCameras', () => {
  it('moves the others with one, never echoing, and lets go on unlink', async () => {
    const [a, b, c] = [await viewer(), await viewer(), await viewer()];
    const [ma, mb, mc] = [mount(a), mount(b), mount(c)];
    const heard = { a: 0, b: 0, c: 0 };
    a.camera.changed.subscribe(() => heard.a++);
    b.camera.changed.subscribe(() => heard.b++);
    c.camera.changed.subscribe(() => heard.c++);
    // Starting together: a's view goes to the rest.
    a.camera.set(at(2), { silent: true });
    const unlink = syncCameras([a, b, c], a);
    expect([mb, mc].map((m) => m.camera.state.target[0])).toEqual([2, 2]);

    // b is dragged: a and c follow it, and nobody announces what was set on it.
    Vec3.set(mb.camera.state.position, 11, 0, 100);
    Vec3.set(mb.camera.state.target, 11, 0, 0);
    mb.frame();
    ma.frame();
    mc.frame();
    mb.frame();
    expect([ma, mb, mc].map((m) => m.camera.state.target[0])).toEqual([11, 11, 11]);
    expect(heard).toEqual({ a: 0, b: 1, c: 0 });

    // Unlinked: b moves alone.
    unlink();
    Vec3.set(mb.camera.state.target, 20, 0, 0);
    mb.frame();
    expect([ma, mb, mc].map((m) => m.camera.state.target[0])).toEqual([11, 20, 11]);
  });

  it("keeps a follower's clamp from shrinking the leader's radius", async () => {
    const [a, b] = [await viewer(), await viewer()];
    // b's scene is smaller: its camera clamps the radius to 30.
    const [ma, mb] = [mount(a), mount(b, 30)];
    const heard = { a: 0, b: 0 };
    a.camera.changed.subscribe(() => heard.a++);
    b.camera.changed.subscribe(() => heard.b++);
    const unlink = syncCameras([a, b]);
    Vec3.set(ma.camera.state.position, 4, 0, 100);
    Vec3.set(ma.camera.state.target, 4, 0, 0);
    ma.frame();
    mb.frame();
    ma.frame();
    expect(mb.camera.state.target[0]).toBe(4);
    expect(mb.camera.state.radius).toBe(30);
    expect(ma.camera.state.radius).toBe(50);
    expect(heard).toEqual({ a: 1, b: 0 });
    unlink();
  });
});
