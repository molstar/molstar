/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DContext, Canvas3DParams, type Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import {
  canvasLeaves,
  canvasPathsByGroup,
  SettingsGroupIds,
  SettingsMapping,
  readSettings,
} from '../../sources/settings.js';
import type { PrismViewer } from '../../viewer.js';
import { expectOk, headlessViewer } from '../../_test/harness.js';

// The `settings` area (operations.md §2): the groups over the real params, `get` and `set`, Undo. Behaviors and Viewer
// work headless; the canvas groups need a canvas, which a stand-in provides here (there is no WebGL).

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

/** A viewer with a stand-in canvas and canvas context that keep their props. */
async function withCanvas() {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
  open = v;
  const props = structuredClone(PD.getDefaultValues(Canvas3DParams)) as Canvas3DProps;
  const context = structuredClone(PD.getDefaultValues(Canvas3DContext.Params));
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: {
      get props() {
        return structuredClone(props);
      },
      setProps: (c: any) => merge(props, structuredClone(c)),
      dispose() {},
    },
    configurable: true,
  });
  Object.defineProperty(v.plugin, 'canvas3dContext', {
    value: {
      get props() {
        return structuredClone(context);
      },
      setProps: (c: any) => Object.assign(context, c),
      dispose() {},
    },
    configurable: true,
  });
  v.plugin.setCanvas3DProps = (change) => {
    merge(props, structuredClone(change));
    v.plugin.events.canvas3d.settingsUpdated.next(void 0);
  };
  return { v, props, context };
}

describe('the settings mapping', () => {
  it('puts every Canvas3D param path in exactly one group', () => {
    const byGroup = canvasPathsByGroup();
    const leaves = canvasLeaves();
    expect(leaves.length).toBeGreaterThan(50);
    for (const leaf of leaves) {
      const owners = SettingsGroupIds.filter((g) => byGroup[g].includes(leaf));
      expect(owners, leaf).toHaveLength(1);
    }
    expect(Object.values(byGroup).flat().sort()).toEqual([...leaves].sort());
  });

  it('names only paths that exist, and no path in two groups', () => {
    const leaves = canvasLeaves();
    const named = SettingsMapping.flatMap((g) => g.canvas ?? []);
    expect(new Set(named).size).toBe(named.length);
    for (const entry of named) {
      expect(
        leaves.some((l) => l === entry || l.startsWith(`${entry}.`)),
        entry,
      ).toBe(true);
    }
    // No entry inside another (a node and one of its children).
    for (const a of named) for (const b of named) if (a !== b) expect(b.startsWith(`${a}.`), `${a} / ${b}`).toBe(false);
  });

  it('gives each canvas group a form over those params', async () => {
    const { v } = await withCanvas();
    const groups = readSettings(v.plugin, true);
    const ids = groups.map((g) => g.id);
    for (const id of SettingsGroupIds) if (id !== 'behaviors') expect(ids).toContain(id);
    // Every Canvas3D param is in some group's form, whole or in part.
    const keys = groups.flatMap((g) => g.canvas);
    expect([...new Set(keys)].sort()).toEqual(Object.keys(Canvas3DParams).sort());
    // The forms follow the canvas: every group's values are valid for its params.
    for (const g of groups) expect(Object.keys(g.values).sort(), g.id).toEqual(Object.keys(g.params).sort());
  });
});

describe('settings, headless', () => {
  it('lists the groups the host can show', async () => {
    const v = await headlessViewer();
    const r = await v.api.settings.get();
    expectOk(r);
    const ids = (r.value ?? []).map((g) => g.id);
    expect(ids).toContain('viewer');
    expect(ids).toContain('behaviors');
    expect(ids).not.toContain('appearance');
    expect(r.changed).toEqual({});
    const viewer = await v.api.settings.get({ group: 'viewer' });
    expectOk(viewer);
    expect(viewer.value?.[0].values).toEqual({ checkpointing: true });
    expect(await v.api.settings.get({ group: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('changes a behavior, is already so the second time, and Undo restores it', async () => {
    const v = await headlessViewer();
    const got = await v.api.settings.get({ group: 'behaviors' });
    expectOk(got);
    const values = got.value![0].values as Record<string, Record<string, unknown>>;
    const key = Object.keys(values).find((k) => typeof values[k].mark === 'boolean');
    expect(key, Object.keys(values).join()).toBeDefined();
    const before = values[key!].mark as boolean;
    const steps = () => v.managers.history.steps.value.length;
    expectOk(await v.api.settings.set({ group: 'behaviors', values: { [key!]: { mark: !before } } }));
    const read = async () =>
      ((await v.api.settings.get({ group: 'behaviors' })) as any).value[0].values[key!].mark as boolean;
    expect(await read()).toBe(!before);
    expect(steps()).toBe(1);
    const again = await v.api.settings.set({ group: 'behaviors', values: { [key!]: { mark: !before } } });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(steps()).toBe(1);
    expectOk(await v.api.history.undo());
    expect(await read()).toBe(before);
  });

  it('rejects unknown settings and bad values, changing nothing', async () => {
    const v = await headlessViewer();
    expect(await v.api.settings.set({ group: 'behaviors', values: { nope: {} } })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.settings.set({ group: 'viewer', values: { checkpointing: 'yes' } })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it('turns per-step undo off and on from the Viewer group', async () => {
    const v = await headlessViewer();
    expect(v.managers.history.checkpointing.value).toBe(true);
    expectOk(await v.api.settings.set({ group: 'viewer', values: { checkpointing: false } }));
    expect(v.managers.history.checkpointing.value).toBe(false);
    expect(await v.api.settings.set({ group: 'viewer', values: { checkpointing: false } })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expectOk(await v.api.settings.set({ group: 'viewer', values: { checkpointing: true } }));
    expect(v.managers.history.checkpointing.value).toBe(true);
  });

  it('fails unavailable for a canvas group without a canvas', async () => {
    const v = await headlessViewer();
    for (const group of ['appearance', 'camera', 'interaction', 'performance', 'advanced'] as const) {
      expect(await v.api.settings.set({ group, values: {} }), group).toMatchObject({ ok: false, code: 'unavailable' });
    }
  });
});

describe('settings, with a canvas', () => {
  it('sets canvas values, merges partial ones, and Undo restores them', async () => {
    const { v, props } = await withCanvas();
    const before = structuredClone(props);
    const r = await v.api.settings.set({
      group: 'appearance',
      values: { renderer: { backgroundColor: 0x112233 }, transparentBackground: true, cameraFog: { name: 'off' } },
    });
    expectOk(r);
    expect(r.changed.view).toBe(true);
    expect(props.renderer.backgroundColor).toBe(0x112233);
    // The rest of the renderer group is untouched.
    expect(props.renderer.exposure).toBe(before.renderer.exposure);
    expect(props.transparentBackground).toBe(true);
    expect(props.cameraFog.name).toBe('off');
    expect(v.managers.history.steps.value).toHaveLength(1);
    const again = await v.api.settings.set({
      group: 'appearance',
      values: { renderer: { backgroundColor: 0x112233 } },
    });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(v.managers.history.steps.value).toHaveLength(1);
    // A mapped value's params follow its name: 'on' takes the defaults when it was off.
    expectOk(await v.api.settings.set({ group: 'appearance', values: { cameraFog: { name: 'on' } } }));
    expect(props.cameraFog).toEqual(before.cameraFog);
    expectOk(await v.api.history.undo());
    expectOk(await v.api.history.undo());
    expect(props.renderer.backgroundColor).toBe(before.renderer.backgroundColor);
    expect(props.transparentBackground).toBe(before.transparentBackground);
    expect(props.cameraFog.name).toBe(before.cameraFog.name);
  });

  it('reads the values back, group by group', async () => {
    const { v, props } = await withCanvas();
    expectOk(await v.api.settings.set({ group: 'camera', values: { camera: { fov: 60 }, sceneRadiusFactor: 2 } }));
    const got = await v.api.settings.get({ group: 'camera' });
    expectOk(got);
    const values = got.value![0].values as any;
    expect(values.camera.fov).toBe(60);
    expect(values.sceneRadiusFactor).toBe(2);
    expect(props.camera.fov).toBe(60);
    expect(await v.api.settings.set({ group: 'camera', values: { camera: { fov: 5000 } } })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it('sets the canvas context in Performance, and says Undo does not restore it', async () => {
    const { v, context } = await withCanvas();
    const r = await v.api.settings.set({ group: 'performance', values: { pixelScale: 0.5, hiZ: { enabled: true } } });
    expectOk(r);
    expect(context.pixelScale).toBe(0.5);
    expect(r.warnings?.join(' ')).toMatch(/Undo/);
    expectOk(await v.api.history.undo());
    expect(context.pixelScale).toBe(0.5);
  });

  it('sets the picking level in Interaction', async () => {
    const { v } = await withCanvas();
    expectOk(await v.api.settings.set({ group: 'interaction', values: { granularity: 'chain' } }));
    expect(v.plugin.managers.interactivity.props.granularity).toBe('chain');
  });
});
