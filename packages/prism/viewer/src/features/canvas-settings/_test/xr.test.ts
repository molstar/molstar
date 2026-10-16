/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Subject } from 'rxjs';
import { afterEach, describe, expect, it } from 'vitest';
import { expectOk, headlessViewer } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import { PrismViewerOptions, resolveOptions } from '../../../options.js';
import { createXrSource } from '../../../sources/xr.js';
import type { PrismViewer } from '../../../viewer.js';

// `canvas.xr` (WebXR, Canvas3D's `xr`) on a stand-in canvas: there is no WebGL or headset here.

let open: PrismViewer | undefined;
afterEach(async () => {
  await open?.dispose();
  open = undefined;
});

async function withXr(supported: boolean, fail?: string) {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas'] });
  open = v;
  const xr = {
    isSupported: new BehaviorSubject(supported),
    isPresenting: new BehaviorSubject(false),
    requestFailed: new Subject<string>(),
    async request() {
      if (fail) xr.requestFailed.next(fail);
      else xr.isPresenting.next(true);
    },
    async end() {
      xr.isPresenting.next(false);
    },
  };
  Object.defineProperty(v.plugin, 'canvas3d', { value: { xr, dispose() {} }, configurable: true });
  return { v, xr };
}

describe('canvas.xr', () => {
  it('enters and leaves XR, no history step; already-so', async () => {
    const { v, xr } = await withXr(true);
    expect(v.api.describe('canvas.xr')?.available).toBe(true);
    expectOk(await v.api.canvas.xr!({ on: true }));
    expect(xr.isPresenting.value).toBe(true);
    expect(await v.api.canvas.xr!({ on: true })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.canvas.xr!({ on: false }));
    expect(xr.isPresenting.value).toBe(false);
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it('is blocked where WebXR is not, and says why a request failed', async () => {
    const none = await withXr(false);
    expect(none.v.api.describe('canvas.xr')?.available).toMatch(/not supported/);
    expect(await none.v.api.canvas.xr!({})).toMatchObject({ ok: false });
    await none.v.dispose();
    open = undefined;

    const refused = await withXr(true, 'The person declined');
    expect(await refused.v.api.canvas.xr!({ on: true })).toMatchObject({
      ok: false,
      code: 'failed',
      message: 'Could not enter XR: The person declined',
    });
  });

  it('needs a canvas; its command is there', async () => {
    const v = await headlessViewer();
    expect(await v.api.canvas.xr!({})).toMatchObject({ ok: false, code: 'unavailable' });
    expect(v.managers.commands.get('canvas.xr')).toBeDefined();
  });

  it('the XR source follows the canvas', async () => {
    const { v, xr } = await withXr(true);
    (v.plugin.behaviors.canvas3d.initialized as BehaviorSubject<boolean>).next(true);
    const source = createXrSource(v.plugin);
    expect(source.state.value).toEqual({ supported: true, presenting: false });
    xr.isPresenting.next(true);
    expect(source.state.value).toEqual({ supported: true, presenting: true });
    source.dispose();
  });

  it("the XR section sets Canvas3D's xr params through canvas.settings", async () => {
    const { v } = await withXr(true);
    const props = {
      xr: { minTargetDistance: 0.4, disablePostprocessing: true, resolutionScale: 1, sceneRadiusInMeters: 0.25 },
    };
    const sets: unknown[] = [];
    Object.defineProperty(v.plugin, 'canvas3d', {
      value: {
        ...v.plugin.canvas3d,
        get props() {
          return structuredClone(props);
        },
      },
      configurable: true,
    });
    v.plugin.setCanvas3DProps = (c: any) => {
      sets.push(c);
      Object.assign(props.xr, c.xr);
    };
    expectOk(await v.api.run('canvas.settings', { section: 'xr', values: { resolutionScale: 1.5 } }));
    expect(sets).toEqual([{ xr: { ...props.xr, resolutionScale: 1.5 } }]);
    expect(await v.api.run('canvas.settings', { section: 'xr', values: { resolutionScale: 1.5 } })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expect(await v.api.run('canvas.settings', { section: 'xr', values: { resolutionScale: 9 } })).toMatchObject({
      ok: false,
    });
  });

  it('the ui.xr option: auto by default; always and never; others refused', () => {
    expect(resolveOptions().ui.xr).toBe('auto');
    expect(resolveOptions({ ui: { xr: 'always' } }).ui.xr).toBe('always');
    expect(
      (PrismViewerOptions.ui.params.xr as unknown as { options: [string, string][] }).options.map((o) => o[0]),
    ).toEqual(['auto', 'always', 'never']);
  });
});
