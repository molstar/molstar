/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { degToRad } from '@molstar/core/math/misc';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import type { PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { afterEach, describe, expect, it } from 'vitest';
import { expectOk, headlessViewer } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';
import { parseVec } from '../camera.js';

// `camera.get` and `camera.set` (v5's camera info, parity.md) on a stand-in canvas (no WebGL here); `view.expand` and
// `view.fullscreen` (the navigation state the page's shell follows).

let open: PrismViewer | undefined;
afterEach(async () => {
  await open?.dispose();
  open = undefined;
});

/** A viewer with a stand-in canvas whose camera keeps what it is set to (as `requestCameraReset` would). */
async function withCamera() {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
  open = v;
  const state = {
    mode: 'perspective' as Camera.Mode,
    fov: degToRad(45),
    position: Vec3.create(0, 0, 100),
    target: Vec3.create(0, 0, 0),
    up: Vec3.create(0, 1, 0),
    radius: 50,
  };
  const snapshots: Partial<Camera.Snapshot>[] = [];
  const props: PartialCanvas3DProps[] = [];
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: {
      camera: { state, getSnapshot: () => ({ ...state }) },
      boundingSphere: { center: Vec3(), radius: 40 },
      dispose() {},
    },
    configurable: true,
  });
  v.plugin.managers.camera.setSnapshot = (s) => {
    snapshots.push(s);
    Object.assign(state, s);
  };
  v.plugin.setCanvas3DProps = (c) => {
    props.push(c);
    if (c.camera?.mode) state.mode = c.camera.mode;
    if (c.camera?.fov) state.fov = degToRad(c.camera.fov);
  };
  return { v, state, snapshots, props };
}

describe('camera.get and camera.set', () => {
  it('reads the camera as values and words', async () => {
    const { v } = await withCamera();
    const r = await v.api.camera.get();
    expectOk(r);
    expect(r.value).toMatchObject({
      mode: 'perspective',
      fov: 45,
      position: [0, 0, 100],
      target: [0, 0, 0],
      up: [0, 1, 0],
      distance: 100,
      showsAll: true,
    });
    expect(r.did).toMatch(/Perspective camera \(45°\).*at \(0, 0, 100\), looking at \(0, 0, 0\)/);
  });

  it('sets the position, target and up; the lens through the canvas props; no history step', async () => {
    const { v, state, snapshots, props } = await withCamera();
    const r = await v.api.camera.set({ position: '10, 0, 0', target: '[0, 0, 0]', up: '0 0 2', fov: 60 });
    expectOk(r);
    expect(r.did).toBe('Camera: at (10, 0, 0), up (0, 0, 1), 60°.');
    expect(Array.from(state.position)).toEqual([10, 0, 0]);
    expect(Array.from(state.up)).toEqual([0, 0, 1]);
    // The radius follows the distance (50 of 100 before: 5 of 10 now).
    expect(snapshots[0].radius).toBeCloseTo(5);
    expect(props).toEqual([{ camera: { fov: 60 } }]);
    expect(v.managers.history.steps.value).toHaveLength(0);

    const again = await v.api.camera.set({ position: '10, 0, 0', up: '0, 0, 1', fov: 60 });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(snapshots).toHaveLength(1);
  });

  it('refuses what is not a vector, a zero up, and a camera at its target', async () => {
    const { v } = await withCamera();
    expect(await v.api.camera.set({ position: '1, 2' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.camera.set({ up: '0, 0, 0' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.camera.set({ position: '0, 0, 0' })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('says so without a canvas', async () => {
    const v = await headlessViewer();
    expect(await v.api.camera.set({ fov: 30 })).toMatchObject({ ok: false, code: 'unavailable' });
    const r = await v.api.camera.get();
    expectOk(r);
    expect(r.value).toBeUndefined();
  });

  it('parses vectors loosely', () => {
    expect(parseVec('')).toBeUndefined();
    expect(Array.from(parseVec('1,2,3')!)).toEqual([1, 2, 3]);
    expect(Array.from(parseVec('[1; 2  3]')!)).toEqual([1, 2, 3]);
    expect(parseVec('1, x, 3')).toBeNull();
  });
});

describe('view.expand and view.fullscreen', () => {
  it('expands and back, already-so, and the card closes with it', async () => {
    const v = await headlessViewer();
    const nav = v.managers.navigation;
    nav.setCard({ id: 'core.overlay.export' });
    expectOk(await v.api.view.expand({ on: 'on' }));
    expect(nav.expanded.value).toBe(true);
    expect(nav.state.value.card).toBeUndefined();
    expect(await v.api.view.expand({ on: 'on' })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.view.expand({ on: 'toggle' }));
    expect(nav.expanded.value).toBe(false);
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it('opening a panel leaves the expanded view', async () => {
    const v = await headlessViewer();
    expectOk(await v.api.view.expand({}));
    v.managers.navigation.setDock('right', { active: 'home', collapsed: false });
    expect(v.managers.navigation.expanded.value).toBe(false);
  });

  it('fullscreen needs a page', async () => {
    const v = await headlessViewer();
    expect(await v.api.view.fullscreen({})).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('fullscreen asks through the navigation state in a page', async () => {
    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['page'] });
    open = v;
    expectOk(await v.api.view.fullscreen({ on: 'on' }));
    expect(v.managers.navigation.fullscreen.value).toBe(true);
    expect(await v.api.view.fullscreen({ on: 'on' })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.view.fullscreen({}));
    expect(v.managers.navigation.fullscreen.value).toBe(false);
  });
});
