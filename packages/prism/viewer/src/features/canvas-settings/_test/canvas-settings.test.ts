/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DParams, type Canvas3DProps, type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { afterEach, describe, expect, it } from 'vitest';
import { expectOk, headlessViewer } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { CanvasSettingsFeature } from '../index.js';
import { type CanvasSection, CanvasSections, sectionProps, sectionValues } from '../operations.js';

// The `canvas-settings` feature (viewer-spec.md §3.1): its operations on a stand-in canvas (there is no WebGL here).

/** Merges a partial change into props, as Canvas3D's `setProps` does for these fields. */
function merge(target: any, change: any) {
  for (const [k, v] of Object.entries(change)) {
    if (v === undefined) continue;
    if (
      v &&
      typeof v === 'object' &&
      !Array.isArray(v) &&
      !('name' in v) &&
      target[k] &&
      typeof target[k] === 'object'
    ) {
      merge(target[k], v);
    } else target[k] = v;
  }
}

let open: PrismViewer | undefined;
afterEach(async () => {
  await open?.dispose();
  open = undefined;
});

/** A viewer whose host says it has a canvas, and a stand-in canvas that keeps props. */
async function withCanvas(): Promise<{ v: PrismViewer; props: Canvas3DProps; sets: PartialCanvas3DProps[] }> {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
  open = v;
  const props = PD.getDefaultValues(Canvas3DParams) as Canvas3DProps;
  const sets: PartialCanvas3DProps[] = [];
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: { props, setProps: (c: PartialCanvas3DProps) => merge(props, c), dispose() {} },
    configurable: true,
  });
  v.plugin.setCanvas3DProps = (change) => {
    sets.push(change);
    merge(props, change);
  };
  return { v, props, sets };
}

describe('canvas settings', () => {
  it('is in the default spec: its operations in the canvas area, and its commands', async () => {
    const v = await headlessViewer();
    for (const id of [
      'canvas.settings',
      'canvas.lighting',
      'canvas.illumination',
      'canvas.camera',
      'canvas.clipping',
    ]) {
      expect(v.api.describe(id)).toBeDefined();
    }
    expect(v.api.canvas.lighting).toBeDefined();
    expect(v.managers.commands.get('canvas.lighting')).toBeDefined();
    // No canvas, nothing to change.
    expect(await v.api.canvas.lighting!({ shadow: 'on' })).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('turns lighting effects on and off, as one step; already so when they are', async () => {
    const { v, props } = await withCanvas();
    expect(props.postprocessing.shadow.name).toBe('off');
    expectOk(await v.api.canvas.lighting!({ shadow: 'on', fog: 'off' }));
    expect(props.postprocessing.shadow.name).toBe('on');
    expect(props.cameraFog.name).toBe('off');
    const again = await v.api.canvas.lighting!({ shadow: 'on', fog: 'off' });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    const none = await v.api.canvas.lighting!({});
    expect(none, JSON.stringify(none)).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('sets illumination, the camera and clipping', async () => {
    const { v, props } = await withCanvas();
    expectOk(await v.api.canvas.illumination!({ on: true }));
    expect(props.illumination.enabled).toBe(true);
    expect(await v.api.canvas.illumination!({ on: true })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.canvas.camera!({ mode: 'orthographic' }));
    expect(props.camera.mode).toBe('orthographic');
    expectOk(await v.api.canvas.clipping!({ radius: 50, far: false }));
    expect(props.cameraClipping).toMatchObject({ radius: 50, far: false });
  });

  it("applies a card section's values, and reads them back", async () => {
    const { v, props } = await withCanvas();
    expectOk(await v.api.run('canvas.settings', { section: 'background', values: { transparent: true } }));
    expect(props.transparentBackground).toBe(true);
    const same = await v.api.run('canvas.settings', { section: 'background', values: { transparent: true } });
    expect(same).toMatchObject({ ok: true, unchanged: true });
    // Every section's values map to props and back.
    const defaults = PD.getDefaultValues(Canvas3DParams) as Canvas3DProps;
    for (const section of Object.keys(CanvasSections) as CanvasSection[]) {
      const values = sectionValues(section, defaults);
      const back = { ...defaults };
      merge(back, sectionProps(section, values));
      expect(sectionValues(section, back)).toEqual(values);
    }
  });

  it('is only for built-in features', async () => {
    const { builtIn: _, ...outside } = CanvasSettingsFeature;
    await expect(headlessViewer({ extensions: [{ ...outside, id: 'outside' }] })).rejects.toThrow(/only the viewer/);
  });
});
