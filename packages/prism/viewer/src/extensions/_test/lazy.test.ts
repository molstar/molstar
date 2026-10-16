/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs, expectOk } from '../../_test/harness.js';
import { ApiArea, type ApiAreaContext, ApiResult, type OperationDefinition } from '../../api/area.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { LazyViewerExtension, PrismViewerSpec, ViewerExtension, ViewerExtensionUI } from '../../spec.js';
import type { PrismViewer, ViewerCore } from '../../viewer.js';
import type { ExtensionUIPort } from '../host.js';

// Lazy loading, all or nothing (viewer-spec.md §3–§4; gate G3/G4): stubs in core areas, a failed activation undone, a
// loaded extension held to its stub, overrides for every part, files loading what owns them.

const open: PrismViewer[] = [];
afterEach(async () => {
  for (const v of open.splice(0)) await v.dispose();
});

async function viewer(o: Partial<PrismViewerSpec> = {}): Promise<ViewerCore> {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, ...o, extensions: [...(base.extensions ?? []), ...(o.extensions ?? [])] },
  });
  open.push(v);
  return v as ViewerCore;
}

const xrOperation: OperationDefinition = {
  id: 'canvas.stereo',
  title: 'Enter XR',
  params: { on: PD.Boolean(true) },
  async run({ on }) {
    return ApiResult.ok(`XR ${on ? 'on' : 'off'}.`, {});
  },
};

/** A built-in feature, lazy, that adds `canvas.stereo` to the core's `canvas` area. */
function lazyXr(o: { ext?: Partial<ViewerExtension>; lazy?: Partial<LazyViewerExtension> } = {}) {
  const counter = { loads: 0 };
  const ext: ViewerExtension = {
    id: 'xr',
    title: 'XR',
    builtIn: true,
    contributes: () => [{ area: 'canvas', name: 'stereo', operation: xrOperation }],
    commands: (api) => [
      {
        id: 'xr.enter',
        title: 'Enter XR',
        category: 'Canvas',
        run: async () => ({ message: (await api.run('canvas.stereo', {})).ok ? 'ok' : 'no' }),
      },
    ],
    ...o.ext,
  };
  const lazy: LazyViewerExtension = {
    id: 'xr',
    title: 'XR',
    builtIn: true,
    stub: {
      parts: [],
      operations: [{ id: 'canvas.stereo', title: 'Enter XR' }],
      commands: [{ id: 'xr.enter', title: 'Enter XR' }],
    },
    load: async () => {
      counter.loads++;
      return ext;
    },
    ...o.lazy,
  };
  return { lazy, counter };
}

class ClashApi extends ApiArea {
  readonly id = 'camera';
}

describe('lazy extensions', () => {
  it("stub a built-in feature's operations in a core area, and replay the first call", async () => {
    const { lazy, counter } = lazyXr();
    const v = await viewer({ extensions: [lazy] });
    expect(v.api.describe('canvas.stereo')).toBeDefined();
    const r = await v.api.run('canvas.stereo', { on: false });
    expectOk(r);
    expect(r.did).toBe('XR off.');
    expect(counter.loads).toBe(1);
    expect(v.extensions.status('xr')?.value).toBe('loaded');
    // The real one now, its params described.
    expect(Object.keys(v.api.describe('canvas.stereo')?.params ?? {})).toEqual(['on']);
  });

  it('refuse a stub in a core area from an extension that is not built in', async () => {
    const { lazy } = lazyXr({ lazy: { builtIn: undefined } });
    await expect(viewer({ extensions: [lazy] })).rejects.toThrow(/only the viewer's features/);
  });

  it('undo a failed activation and put the stub back; a retry tries again', async () => {
    let clash = true;
    let loads = 0;
    const { lazy } = lazyXr();
    const ext = await lazy.load();
    const clashing: LazyViewerExtension = {
      ...lazy,
      stub: { ...lazy.stub, operations: [...(lazy.stub.operations ?? []), { id: 'camera.xr', title: 'XR camera' }] },
      load: async () => {
        loads++;
        // First it brings an area that exists (activation fails half way); then none for `camera.xr`.
        return clash ? { ...ext, api: { camera: (c: ApiAreaContext) => new ClashApi(c) } } : ext;
      },
    };
    const v = await viewer({ extensions: [clashing] });
    const failed: string[] = [];
    v.extensions.events.failed.subscribe((f) => failed.push(f.message));
    const r = await v.api.run('canvas.stereo', {});
    expect(r).toMatchObject({ ok: false, code: 'unavailable' });
    expect(v.extensions.status('xr')?.value).toBe('failed');
    expect(failed[0]).toMatch(/has already/);
    // The stub is back: its operations and command; nothing of the real one stayed.
    expect(v.api.describe('canvas.stereo')?.params).toEqual({});
    expect(v.api.describe('camera.xr')).toBeDefined();
    expect(v.managers.commands.get('xr.enter')).toBeDefined();
    expect(v.managers.history.steps.value).toHaveLength(0);
    expect(await v.extensions.load('xr')).toBe(false);
    expect(loads).toBe(2);
    clash = false;
    // Now it brings no `camera.xr`: still not what its stub declared, and no endless replay.
    expect(await v.api.run('camera.xr' as never, {})).toMatchObject({ ok: false, code: 'unavailable' });
    expect(failed.at(-1)).toMatch(/does not bring what its stub declared: camera\.xr/);
  });

  it('fail a load that brings a command its stub did not declare', async () => {
    const { lazy } = lazyXr({
      ext: {
        commands: () => [{ id: 'xr.other', title: 'Other', category: 'Canvas', run: async () => ({ message: '' }) }],
      },
    });
    const v = await viewer({ extensions: [lazy] });
    expect(await v.extensions.load('xr')).toBe(false);
    expect(v.managers.commands.get('xr.other')).toBeUndefined();
    expect(v.managers.commands.get('xr.enter')).toBeDefined();
  });

  it('declare the flags of their operations before the load, and fail a load whose flags differ', async () => {
    const { lazy } = lazyXr();
    const v = await viewer({
      extensions: [
        {
          ...lazy,
          stub: { ...lazy.stub, operations: [{ id: 'canvas.stereo', title: 'Enter XR', approval: 'destructive' }] },
        },
      ],
    });
    expect(v.api.describe('canvas.stereo')?.definition.approval).toBe('destructive');
    expect(await v.extensions.load('xr')).toBe(false);
  });

  it('fail a load that throws, say why, and load on a retry', async () => {
    let fail = true;
    const { lazy } = lazyXr();
    const flaky: LazyViewerExtension = {
      ...lazy,
      load: async () => {
        if (fail) throw new Error('offline');
        return lazy.load();
      },
    };
    const v = await viewer({ extensions: [flaky] });
    const r = await v.api.run('canvas.stereo', {});
    expect(r).toMatchObject({ ok: false, code: 'unavailable' });
    expect(v.managers.history.steps.value).toHaveLength(0);
    fail = false;
    expectOk(await v.api.run('canvas.stereo', {}));
    expect(v.extensions.status('xr')?.value).toBe('loaded');
  });

  it('load the extension that owns a file before opening it', async () => {
    const { lazy, counter } = lazyXr({
      lazy: {
        stub: {
          parts: [],
          operations: [{ id: 'canvas.stereo', title: 'XR' }],
          commands: [{ id: 'xr.enter', title: 'Enter XR' }],
          owns: { fileExtensions: ['.CIF'] },
        },
      },
    });
    const v = await viewer({ extensions: [lazy] });
    expectOk(await v.api.load.files({ files: [cbs()] }));
    expect(counter.loads).toBe(1);
    expect(v.extensions.status('xr')?.value).toBe('loaded');
  });

  it('open the real command after a stub command loads it', async () => {
    const { lazy } = lazyXr();
    const v = await viewer({ extensions: [lazy] });
    const opened = v.managers.commands.open('xr.enter');
    expect(opened.kind).toBe('ran');
    if (opened.kind === 'ran') expect((await opened.result)?.message).toBe('Loaded XR.');
    await new Promise((r) => setTimeout(r, 10));
    // The real command ran: its run calls canvas.xr, which records a step.
    await new Promise((r) => setTimeout(r, 10));
    expect(v.managers.history.steps.value.map((s) => s.label)).toEqual(['XR on']);
  });
});

describe('overrides', () => {
  it("remove or replace the core's commands and extensions' commands by id", async () => {
    const { lazy } = lazyXr();
    const v = await viewer({
      extensions: [lazy],
      overrides: { 'camera.reset': null, 'xr.enter': { title: 'Go immersive' } },
    });
    expect(v.managers.commands.get('camera.reset')).toBeUndefined();
    expect(v.managers.commands.get('xr.enter')?.title).toBe('Go immersive');
    await v.extensions.load('xr');
    expect(v.managers.commands.get('xr.enter')?.title).toBe('Go immersive');
  });

  it('refuse an override of a command or part nobody has', async () => {
    await expect(viewer({ overrides: { 'camera.nope': null } })).rejects.toThrow(/no command 'camera.nope'/);
    await expect(viewer({ overrides: { 'core.tab.nope': null } })).rejects.toThrow(/the core has no such part/);
    const { lazy } = lazyXr();
    await expect(viewer({ extensions: [lazy], overrides: { 'xr.tab.nope': null } })).rejects.toThrow(
      /declares no such part/,
    );
    // The shell's commands are known though a headless viewer has none.
    await viewer({ overrides: { 'help.show': null, 'core.tab.log': null } });
  });

  it("apply to an extension's views and a lazy one's placeholders; views are held to their ids", async () => {
    const views: ViewerExtensionUI = {
      sections: [{ id: 'panel.section.a', slot: 'export', title: 'A', render: { mount: () => () => {} } }],
      toolbar: [
        { id: 'panel.toolbar.a', icon: 'sun', title: 'A' },
        { id: 'panel.toolbar.b', icon: 'sun', title: 'B' },
        { id: 'other.toolbar.c', icon: 'sun', title: 'C' },
      ],
    };
    const panel: ViewerExtension = { id: 'panel', title: 'Panel', ui: () => views };
    const lazyPanel: LazyViewerExtension = {
      id: 'later',
      title: 'Later',
      stub: {
        parts: [
          { id: 'later.tab.t', kind: 'tab', title: 'T' },
          { id: 'later.tab.u', kind: 'tab', title: 'U' },
        ],
      },
      load: async () => ({ id: 'later', title: 'Later' }),
    };
    const v = await viewer({
      extensions: [panel, lazyPanel],
      overrides: { 'panel.toolbar.b': null, 'panel.section.a': { title: 'Renamed' }, 'later.tab.u': null },
    });
    const registered: ViewerExtensionUI[] = [];
    const stubs: string[][] = [];
    const failed: string[] = [];
    v.extensions.events.failed.subscribe((f) => failed.push(f.message));
    const port: ExtensionUIPort = {
      context: { api: v.api, sources: v.sources, managers: v.managers },
      register: (_, ui) => {
        registered.push(ui);
        return () => {};
      },
      registerStub: (owner, stub) => {
        // The default spec's lazy features have placeholders too: this test's own only.
        if (owner === 'later') stubs.push(stub.parts.map((p) => p.id));
        return () => {};
      },
    };
    v.attachExtensionUI(port);
    await new Promise((r) => setTimeout(r, 0));
    expect(registered[0].toolbar?.map((t) => t.id)).toEqual(['panel.toolbar.a']);
    expect(registered[0].sections?.[0].title).toBe('Renamed');
    expect(stubs).toEqual([['later.tab.t']]);
    expect(failed[0]).toMatch(/other\.toolbar\.c/);
  });
});

describe('contributed operations', () => {
  it('are undefined when removed', async () => {
    const v = await viewer();
    const remove = v.host.api.area('canvas')!.contribute('stereo', xrOperation);
    expect(v.api.describe('canvas.stereo')).toBeDefined();
    remove();
    expect(v.api.describe('canvas.stereo')).toBeUndefined();
    expect(await v.api.run('canvas.stereo' as never, {})).toMatchObject({ ok: false });
  });
});
