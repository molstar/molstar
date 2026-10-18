/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { degToRad } from '@molstar/core/math/misc';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import { Canvas3DParams, type Canvas3DProps, type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { BehaviorSubject, firstValueFrom, skip } from 'rxjs';
import { afterEach, describe, expect, it } from 'vitest';
import { expectOk, headlessViewer } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { readViewport, type ViewportState } from '../operations.js';

// The viewport state (v5's screenshot panel, parity.md): `export.viewport`, `load.viewport` and the camera the Viewport
// section shows, on a stand-in canvas (no WebGL here).

let open: PrismViewer | undefined;
afterEach(async () => {
  await open?.dispose();
  open = undefined;
});

function snapshot(position: number[]): Camera.Snapshot {
  return {
    mode: 'perspective',
    fov: degToRad(45),
    position: Vec3.ofArray(position),
    target: Vec3.create(0, 0, 0),
    up: Vec3.create(0, 1, 0),
    radius: 50,
    radiusMax: 100,
    fog: 50,
    clipFar: true,
    minNear: 5,
    minFar: 0,
  };
}

async function withCanvas() {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
  open = v;
  const props = structuredClone(PD.getDefaultValues(Canvas3DParams)) as Canvas3DProps;
  let state = snapshot([0, 0, 100]);
  const didDraw = new BehaviorSubject(0);
  const sets: PartialCanvas3DProps[] = [];
  const cameras: Partial<Camera.Snapshot>[] = [];
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: {
      get props() {
        return structuredClone(props);
      },
      camera: {
        get state() {
          return state;
        },
        getSnapshot: () => ({ ...state }),
      },
      boundingSphere: { center: Vec3(), radius: 40 },
      didDraw,
      // A failed step puts the canvas back from its checkpoint.
      setProps() {},
      dispose() {},
    },
    configurable: true,
  });
  v.plugin.setCanvas3DProps = (c) => {
    sets.push(c);
  };
  v.plugin.managers.camera.setSnapshot = (s) => {
    cameras.push(s);
    state = { ...state, ...s };
    didDraw.next(didDraw.value + 1);
  };
  return { v, sets, cameras, didDraw };
}

describe('the viewport state', () => {
  it('exports the camera and the canvas as JSON, and loads it back', async () => {
    const { v, sets, cameras } = await withCanvas();
    const r = await v.api.export!.viewport({ name: 'view' });
    expectOk(r);
    expect(r.value!.filename).toBe('view-viewport.json');
    const saved = JSON.parse(r.value!.json) as ViewportState;
    expect(saved).toMatchObject({ kind: 'molstar-viewport', version: 1, camera: { mode: 'perspective' } });
    expect(saved.canvas.renderer.backgroundColor).toBe(Canvas3DParams.renderer.params.backgroundColor.defaultValue);
    expect(v.managers.history.steps.value).toHaveLength(0);

    const file = new File([r.value!.json], 'view-viewport.json', { type: 'application/json' });
    const loaded = await v.api.load.viewport!({ file });
    expectOk(loaded);
    expect(loaded.did).toBe('Loaded the camera and the canvas settings from view-viewport.json.');
    expect(sets).toHaveLength(1);
    expect(Array.from(cameras[0].position!)).toEqual([0, 0, 100]);
  });

  it('reads a v5 state and an MVS camera, and refuses what is neither', async () => {
    const { v, cameras } = await withCanvas();
    const v5 = { entries: [{ snapshot: { camera: { current: snapshot([5, 0, 0]) } } }] };
    expectOk(await v.api.load.viewport!({ json: JSON.stringify(v5) }));
    expect(Array.from(cameras[0].position!)).toEqual([5, 0, 0]);
    // A file of several snapshots has the camera of the one that was showing, else the first.
    const two = (current: string | undefined) => ({
      current,
      entries: [
        { snapshot: { id: 'a', camera: { current: snapshot([1, 0, 0]) } } },
        { snapshot: { id: 'b', camera: { current: snapshot([2, 0, 0]) } } },
      ],
    });
    expectOk(await v.api.load.viewport!({ json: JSON.stringify(two('b')) }));
    expect(Array.from(cameras[1].position!)).toEqual([2, 0, 0]);
    expectOk(await v.api.load.viewport!({ json: JSON.stringify(two('gone')) }));
    expect(Array.from(cameras[2].position!)).toEqual([1, 0, 0]);
    const mvs = { position: [1, 2, 3], target: [0, 0, 0], up: [0, 0, 1] };
    const r = await v.api.load.viewport!({ json: JSON.stringify(mvs) });
    expectOk(r);
    expect(r.did).toBe('Loaded the camera.');
    expect(await v.api.load.viewport!({ json: '{"a": 1}' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.viewport!({ json: 'nope' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.viewport!({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(readViewport({ kind: 'molstar-viewport', version: 2 })).toMatch(/version 2/);
    expect(readViewport({ position: [0, 0, 0], target: [0, 0, 0] })).toMatch(/at its target/);
  });

  it('the camera the section shows follows the redraws', async () => {
    const { v } = await withCanvas();
    (v.plugin.behaviors.canvas3d.initialized as BehaviorSubject<boolean>).next(true);
    const first = await firstValueFrom(v.api.export!.camera);
    expect(first?.position).toEqual([0, 0, 100]);
    const next = firstValueFrom(v.api.export!.camera.pipe(skip(1)));
    await v.api.camera.set({ position: '0, 0, 50' });
    expect((await next)?.position).toEqual([0, 0, 50]);
  });

  it('reads the camera again only when it moved', async () => {
    const { v, didDraw } = await withCanvas();
    (v.plugin.behaviors.canvas3d.initialized as BehaviorSubject<boolean>).next(true);
    const seen: unknown[] = [];
    const sub = v.api.export!.camera.subscribe((c) => seen.push(c));
    try {
      // Redraws with the camera still (accumulating samples, an animation).
      for (let i = 0; i < 3; i++) didDraw.next(didDraw.value + 1);
      await new Promise((r) => setTimeout(r, 400));
      expect(seen).toHaveLength(1);
      await v.api.camera.set({ position: '0, 0, 50' });
      await new Promise((r) => setTimeout(r, 400));
      expect(seen).toHaveLength(2);
    } finally {
      sub.unsubscribe();
    }
  });

  it('needs a canvas; its commands are there', async () => {
    const v = await headlessViewer();
    expect(await v.api.export!.viewport({})).toMatchObject({ ok: false, code: 'unavailable' });
    expect(await v.api.load.viewport!({ json: '{}' })).toMatchObject({ ok: false, code: 'unavailable' });
    expect(v.managers.commands.get('export.viewport')).toBeDefined();
    expect(v.managers.commands.get('load.viewport')).toBeDefined();
  });
});
