/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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
  // A copy: the defaults' objects are the params' own.
  const props = structuredClone(PD.getDefaultValues(Canvas3DParams)) as Canvas3DProps;
  const sets: PartialCanvas3DProps[] = [];
  // As Canvas3D: what `props` gives is a copy (a checkpoint holding it must not see later changes).
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: {
      get props() {
        return structuredClone(props);
      },
      setProps: (c: PartialCanvas3DProps) => merge(props, structuredClone(c)),
      dispose() {},
    },
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
    for (const id of ['canvas.settings', 'canvas.lighting', 'canvas.camera', 'canvas.clipping']) {
      expect(v.api.describe(id)).toBeDefined();
    }
    expect(v.api.canvas.lighting).toBeDefined();
    expect(v.managers.commands.get('canvas.lighting')).toBeDefined();
    // Illumination is the lighting's switch; its one-click command runs canvas.lighting.
    expect(v.api.describe('canvas.illumination')).toBeUndefined();
    expect(v.managers.commands.get('canvas.illumination')).toBeDefined();
    // No canvas, nothing to change.
    expect(await v.api.canvas.lighting!({ shadow: { name: 'on' } })).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('turns lighting effects on and off, as one step; already so when they are', async () => {
    const { v, props } = await withCanvas();
    expect(props.postprocessing.shadow.name).toBe('off');
    expectOk(await v.api.canvas.lighting!({ shadow: { name: 'on' }, fog: { name: 'off' } }));
    expect(props.postprocessing.shadow.name).toBe('on');
    expect(props.cameraFog.name).toBe('off');
    const again = await v.api.canvas.lighting!({ shadow: { name: 'on' }, fog: { name: 'off' } });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    const none = await v.api.canvas.lighting!({});
    expect(none, JSON.stringify(none)).toMatchObject({ ok: false, code: 'invalid' });
    // On with settings: the given ones set, the rest kept; a wrong name or setting is refused.
    expectOk(await v.api.canvas.lighting!({ shadow: { name: 'on', params: { steps: 3 } } }));
    expect((props.postprocessing.shadow.params as { steps: number }).steps).toBe(3);
    expect(await v.api.canvas.lighting!({ shadow: { name: 'maybe' } as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.canvas.lighting!({ shadow: { name: 'on', params: { nope: 1 } } as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it('sets illumination, the camera and clipping', async () => {
    const { v, props } = await withCanvas();
    expectOk(await v.api.canvas.lighting!({ illumination: { name: 'on' } }));
    expect(props.illumination.enabled).toBe(true);
    expect(await v.api.canvas.lighting!({ illumination: { name: 'on' } })).toMatchObject({ ok: true, unchanged: true });
    expectOk(await v.api.canvas.camera!({ mode: 'orthographic' }));
    expect(props.camera.mode).toBe('orthographic');
    expectOk(await v.api.canvas.clipping!({ clip: 50, far: false }));
    expect(props.cameraClipping).toMatchObject({ radius: 50, far: false });
    expect(await v.api.canvas.clipping!({ clip: 50, far: false })).toMatchObject({ ok: true, unchanged: true });
    // An absent value is kept (api.md §1): far stays off when only the clip changes.
    expectOk(await v.api.canvas.clipping!({ clip: 20 }));
    expect(props.cameraClipping).toMatchObject({ radius: 20, far: false });
    // A field of view below 10° is refused; the same again is already so.
    expect(await v.api.canvas.camera!({ fov: 5 })).toMatchObject({ ok: false, code: 'invalid' });
    expectOk(await v.api.canvas.camera!({ fov: 10 }));
    expect(props.camera.fov).toBe(10);
    expect(await v.api.canvas.camera!({ fov: 10 })).toMatchObject({ ok: true, unchanged: true });
  });

  it("applies a card section's values, and reads them back", async () => {
    const { v, props } = await withCanvas();
    expectOk(await v.api.run('canvas.settings', { section: 'background', values: { transparent: true } }));
    expect(props.transparentBackground).toBe(true);
    const same = await v.api.run('canvas.settings', { section: 'background', values: { transparent: true } });
    expect(same).toMatchObject({ ok: true, unchanged: true });
    // Every section's values map to props and back.
    const defaults = structuredClone(PD.getDefaultValues(Canvas3DParams)) as Canvas3DProps;
    for (const section of Object.keys(CanvasSections) as CanvasSection[]) {
      const values = sectionValues(section, defaults);
      const back = structuredClone(defaults);
      merge(back, sectionProps(section, values));
      expect(sectionValues(section, back)).toEqual(values);
    }
  });

  it('checks the values of a section before it sets them', async () => {
    const { v, sets } = await withCanvas();
    const run = (section: string, values: Record<string, unknown>) =>
      v.api.run('canvas.settings', { section, values } as never);
    expect(await run('lighting', { occlusion: true })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await run('lighting', { nope: 1 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await run('background', { transparent: 'yes' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(sets).toEqual([]);
    expectOk(await run('lighting', { fog: { name: 'off', params: {} } }));
  });

  it("sets the direct volumes' resolution in the advanced section, within its range", async () => {
    const { v, props } = await withCanvas();
    const run = (volumeResolutionScale: number) =>
      v.api.run('canvas.settings', { section: 'advanced', values: { volumeResolutionScale } });
    const full = await run(1);
    expectOk(full);
    expect(props.volumeResolutionScale).toBe(1);
    expect(full.did).toBe('Changed the advanced settings: Volume resolution 1.');
    expect(await run(0.1)).toMatchObject({ ok: false, code: 'invalid' });
    expect(props.volumeResolutionScale).toBe(1);
  });

  it('names the changed fields in words', async () => {
    const { v } = await withCanvas();
    const bg = await v.api.run('canvas.settings', {
      section: 'background',
      values: { color: 0xffffff, transparent: true },
    });
    expectOk(bg);
    expect(bg.did).toBe('Changed the background settings: Background white, Transparent background on.');
    const shadowOn = (Canvas3DParams.postprocessing.params.shadow as PD.Mapped<any>).map('on') as PD.Group<any>;
    const light = await v.api.run('canvas.settings', {
      section: 'lighting',
      values: { shadow: { name: 'on', params: PD.getDefaultValues(shadowOn.params) } },
    });
    expectOk(light);
    expect(light.did).toBe('Changed the lighting settings: Shadows on.');
  });

  it('is only for built-in features', async () => {
    const { builtIn: _, ...outside } = CanvasSettingsFeature;
    await expect(headlessViewer({ extensions: [{ ...outside, id: 'outside' }] })).rejects.toThrow(/only the viewer/);
  });

  it('makes each change one step, nothing when already so, and Undo restores the canvas', async () => {
    const { v, props } = await withCanvas();
    const steps = () => v.managers.history.steps.value.length;
    const before = structuredClone(props);
    expectOk(await v.api.canvas.lighting!({ shadow: { name: 'on' } }));
    expectOk(await v.api.canvas.lighting!({ illumination: { name: 'on' } }));
    expectOk(await v.api.canvas.camera!({ mode: 'orthographic' }));
    expectOk(await v.api.canvas.clipping!({ clip: 30 }));
    expectOk(await v.api.run('canvas.settings', { section: 'background', values: { transparent: true } }));
    expect(steps()).toBe(5);
    await v.api.canvas.lighting!({ illumination: { name: 'on' } });
    expect(steps()).toBe(5);
    for (let i = 0; i < 5; i++) expectOk(await v.api.history.undo());
    expect(props.postprocessing.shadow.name).toBe(before.postprocessing.shadow.name);
    expect(props.illumination.enabled).toBe(before.illumination.enabled);
    expect(props.camera.mode).toBe(before.camera.mode);
    expect(props.cameraClipping.radius).toBe(before.cameraClipping.radius);
    expect(props.transparentBackground).toBe(before.transparentBackground);
  });
});
