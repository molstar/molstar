/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { VolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/behavior';
import type { VolumeServerHeader, VolumeServerInfo } from '@molstar/plugin/behavior/dynamic/volume-streaming/model';
import {
  CreateVolumeStreamingBehavior,
  CreateVolumeStreamingInfo,
} from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer, recordRuns, until } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import { type StreamingChange, nextStreamingParams, streamingOf, streamingSummary } from '../volume.js';
import { otherOf } from '../../testing/scene-groups.js';

// `volume.streaming-view` (parity.md §2): the view, detail and channels of streamed density. The behavior needs a
// canvas (and the server), so the tests make the streaming cells by hand over a stubbed server header: the operation
// changes the behavior's params in the tree, which is what it would hand the behavior.

const sampling: VolumeServerHeader.Sampling = {
  byteOffset: 0,
  rate: 1,
  sampleCount: [10, 10, 10],
  valuesInfo: [
    { mean: 0, min: -1, max: 1, sigma: 0.1 },
    { mean: 0, min: -1, max: 1, sigma: 0.1 },
  ],
};

const header: VolumeServerHeader = {
  formatVersion: '1.0.0',
  axisOrder: [0, 1, 2],
  origin: [0, 0, 0],
  dimensions: [1, 1, 1],
  spacegroup: { number: 1, size: [10, 10, 10], angles: [90, 90, 90], isPeriodic: false },
  channels: ['2FO-FC', 'FO-FC'],
  valueType: 'float32',
  blockSize: 96,
  sampling: [sampling],
  availablePrecisions: [1, 2, 3, 4, 5, 6, 7].map((i) => ({ precision: i, maxVoxels: 2 ** (18 + i) })),
  isAvailable: true,
};

const xray: VolumeServerInfo.EntryData = {
  kind: 'x-ray',
  dataId: '1cbs',
  header,
  emDefaultContourLevel: undefined as never,
};

const defaults = (entry: VolumeServerInfo.EntryData) =>
  PD.getDefaultValues(VolumeStreaming.createParams({ data: { serverUrl: '', entries: [entry] } as never }));

const none: StreamingChange = {};

afterEach(() => vi.unstubAllGlobals());

/** Streaming on 1CBS as `volume.streaming` leaves it, minus the behavior's object (no canvas here). */
async function streamed(): Promise<PrismViewer> {
  const v = await headlessViewer({ files: [cbs()] });
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => new Response(JSON.stringify(header), { status: 200 })),
  );
  const state = v.plugin.state.data;
  const s = v.plugin.managers.structure.hierarchy.current.structures[0];
  const info = state
    .build()
    .to(s.cell)
    .apply(
      CreateVolumeStreamingInfo,
      { serverUrl: 'https://density.test', entries: [{ dataId: '1cbs', source: { name: 'x-ray', params: {} } }] },
      { tags: VolumeStreaming.RootTag },
    );
  await v.plugin.runTask(state.updateTree(info));
  const data = state.cells.get(info.ref)!.obj!.data as VolumeServerInfo.Data;
  const behavior = state
    .build()
    .to(info.ref)
    .apply(CreateVolumeStreamingBehavior, PD.getDefaultValues(VolumeStreaming.createParams({ data })) as never);
  await v.plugin.runTask(state.updateTree(behavior)).catch(() => undefined);
  return v;
}

function behaviorParams(v: PrismViewer): VolumeStreaming.Params {
  const on = streamingOf(v.plugin, v.plugin.managers.structure.hierarchy.current.structures[0]);
  if (!on) throw new Error('no streaming behavior');
  return on.params;
}

describe('nextStreamingParams', () => {
  it('changes the view, radius, detail and channels, and checks them against the entry', () => {
    const current = defaults(xray);
    expect(current.entry.params.view.name).toBe('selection-box');
    expect(nextStreamingParams(current, xray, none)).toEqual(current);

    const next = nextStreamingParams(current, xray, {
      view: 'auto',
      radius: 8,
      detail: 5,
      channels: { '2fo-fc': { isoValue: 1.2, color: '#ff0000', opacity: 0.5, wireframe: true } },
    });
    if ('ok' in next) throw new Error(next.message);
    const p = next.entry.params;
    const channels = p.channels as Record<string, VolumeStreaming.ChannelParams>;
    const was = current.entry.params.channels as Record<string, VolumeStreaming.ChannelParams>;
    expect(p.view.name).toBe('auto');
    expect((p.view.params as { radius: number }).radius).toBe(8);
    expect(p.detailLevel).toBe(4);
    expect(channels['2fo-fc']).toMatchObject({
      isoValue: { kind: 'relative', relativeValue: 1.2 },
      color: Color(0xff0000),
      opacity: 0.5,
      wireframe: true,
    });
    expect(channels['fo-fc(+ve)']).toEqual(was['fo-fc(+ve)']);
    expect(streamingSummary({ entry: xray, params: next })).toBe(
      'Auto, 8 Å, detail 5 of 7; 2Fo-Fc 1.20 σ, Fo-Fc(+ve) 3.00 σ, Fo-Fc(-ve) -3.00 σ',
    );

    // Around the camera, the detail is the dynamic one and the radius a fraction of the view.
    const camera = nextStreamingParams(current, xray, { ...none, view: 'camera-target', radius: 0.3, detail: 2 });
    if ('ok' in camera) throw new Error(camera.message);
    expect(camera.entry.params.view.params).toMatchObject({ radius: 0.3, dynamicDetailLevel: 1 });
    expect(camera.entry.params.detailLevel).toBe(current.entry.params.detailLevel);
  });

  it('refuses what the entry or the view does not have', () => {
    const current = defaults(xray);
    const fail = (change: Record<string, unknown>) =>
      nextStreamingParams(current, xray, { ...none, ...change } as never);
    expect(fail({ view: 'cell', radius: 3 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ view: 'camera-target', radius: 2 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ detail: 8 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ detail: 2.5 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ channels: { em: { isoValue: 1 } } })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ channels: { '2fo-fc': { isoValue: 'high' } } })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ channels: { '2fo-fc': { color: 'not a color' } } })).toMatchObject({ ok: false, code: 'invalid' });
    expect(fail({ channels: { '2fo-fc': { opacity: 2 } } })).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('volume.streaming-view', () => {
  it('refuses where streaming is off', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.volume.streamingView({ view: 'cell' });
    expect(r).toMatchObject({ ok: false, code: 'blocked' });
    expect(v.api.describe('volume.streaming-view')?.definition.approval).toBe('network');
    expect(await v.api.volume.streamingView({ view: 'sideways' as never })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it("updates the behavior's params as a history step; Undo restores them", async () => {
    const v = await streamed();
    const before = behaviorParams(v);
    const steps = v.api.history.steps.value.length;

    const r = await v.api.volume.streamingView({ view: 'cell', detail: 2, channels: { '2fo-fc': { isoValue: 1 } } });
    expectOk(r);
    expect(r.did).toMatch(/^Streamed density for .*: Whole Structure, detail 2 of 7; 2Fo-Fc 1\.00 σ/);
    const after = behaviorParams(v).entry.params;
    expect(after.view.name).toBe('cell');
    expect(after.detailLevel).toBe(1);
    expect(v.api.history.steps.value.length).toBe(steps + 1);

    const same = await v.api.volume.streamingView({ view: 'cell' });
    expectOk(same);
    expect(same.unchanged).toBe(true);

    const got = await v.api.volume.list();
    expectOk(got);
    expect(got.did).toMatch(/Density is streamed for .*\(Whole Structure, detail 2 of 7;/);

    expectOk(await v.api.history.undo());
    expect(behaviorParams(v).entry.params.view.name).toBe(before.entry.params.view.name);
  });
});

const channel = (v: PrismViewer, c: string) =>
  (behaviorParams(v).entry.params.channels as any)[c] as { isoValue: { relativeValue: number }; color: number };

describe('the streaming row', () => {
  it('shows under its structure with the view in its inspector, whose edits are volume.streaming-view calls', async () => {
    const v = await streamed();
    const structure = v.views.scene.roots.value[0];
    const row = await (async () => {
      await until(() => otherOf(structure).some((n) => n.label.value === 'Volume streaming'));
      return otherOf(structure).find((n) => n.label.value === 'Volume streaming')!;
    })();
    expect(row.kind).toBe('extension');
    expect(row.description.value).toMatch(/^Around Focus, 5 Å, detail 4 of 7; 2Fo-Fc 1\.50 σ/);

    const inspection = v.views.scene.inspect(row);
    const form = inspection.form!;
    expect(form.values.value).toMatchObject({ view: 'selection-box', radius: 5, detail: 4 });

    const runs = recordRuns(v);
    form.set('view', 'cell');
    form.set('detail', 2);
    await until(() => runs.some((r) => r.operation === 'volume.streaming-view' && !!r.result));
    expect(runs.filter((r) => r.operation === 'volume.streaming-view' && !!r.result)).toHaveLength(1);
    expect(behaviorParams(v).entry.params.view.name).toBe('cell');
    await until(() => row.description.value?.startsWith('Whole Structure') ?? false);

    // The channels are its levels: 2Fo-Fc, and Fo-Fc with its two signs, in σ of the whole map.
    const [twoFoFc, foFc] = inspection.levels!;
    expect(twoFoFc.title).toBe('2Fo-Fc');
    expect(foFc.levels.value.map((l) => l.id)).toEqual(['fo-fc(+ve)', 'fo-fc(-ve)']);
    const { mean, sigma } = twoFoFc.stats;
    expect((twoFoFc.levels.value[0].value - mean) / sigma).toBeCloseTo(1.5, 5);
    twoFoFc.change('2fo-fc', { value: mean + 1.1 * sigma, color: 0xff0000 });
    twoFoFc.settle();
    await until(() => Math.abs(channel(v, '2fo-fc').isoValue.relativeValue - 1.1) < 1e-6);
    expect(channel(v, '2fo-fc').color).toBe(0xff0000);

    // The eye turns the view off; remove turns streaming off (no network: off needs none).
    row.toggleVisibility();
    await until(() => row.isHidden.value);
    expect(behaviorParams(v).entry.params.view.name).toBe('off');

    // And back on it shows the view the person chose, not the default one.
    row.toggleVisibility();
    await until(() => !row.isHidden.value);
    expect(behaviorParams(v).entry.params.view.name).toBe('cell');
  });
});
