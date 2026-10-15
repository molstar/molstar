/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import type { VolumeInfo } from '../volume.js';
import { mapFile } from './volume-fixtures.js';

// `volume` (operations.md §2, plan 2.5) on a headless viewer, over small CCP4 maps made in memory.

const infos = async (v: PrismViewer) => {
  const r = await v.api.volume.get();
  expectOk(r);
  return r.value as VolumeInfo[];
};

const reps = (v: PrismViewer) =>
  v.plugin.managers.volume.hierarchy.current.volumes.flatMap((x) =>
    x.representations.map((r) => (r.cell.transform.params as any).type.name as string),
  );

describe('volume', () => {
  it('describes the volumes: label, grid facts, representations with their iso values', async () => {
    const v = await headlessViewer();
    const none = await v.api.volume.get();
    expectOk(none);
    expect(none.did).toMatch(/No volumes/);

    await v.api.load.files({ files: [mapFile()] });
    const [map] = await infos(v);
    expect(map.label).toBe('test.ccp4');
    expect(map.stats.grid).toEqual([12, 12, 12]);
    expect(map.stats.voxelSize[0]).toBeCloseTo(1.5, 3);
    expect(map.stats.sigma).toBeGreaterThan(0);
    expect(map.representations).toHaveLength(1);
    expect(map.representations[0].type).toBe('isosurface');
    expect(map.representations[0].isoValue?.kind).toBe('relative');
    const r = await v.api.volume.get();
    expectOk(r);
    expect(r.did).toMatch(/test\.ccp4 .*12 × 12 × 12 grid.*isosurface/);
  });

  it('sets the iso value, relative or absolute, and says when it is already so', async () => {
    const v = await headlessViewer();
    expect((await v.api.volume.isovalue({ value: 1 })).ok).toBe(false);
    await v.api.load.files({ files: [mapFile()] });

    expectOk(await v.api.volume.isovalue({ value: 1, kind: 'relative' }));
    let iso = (await infos(v))[0].representations[0].isoValue!;
    expect(iso).toMatchObject({ kind: 'relative', relative: 1 });

    const again = await v.api.volume.isovalue({ value: 1, kind: 'relative' });
    expectOk(again);
    expect(again.unchanged).toBe(true);

    expectOk(await v.api.volume.isovalue({ value: 0.5, kind: 'absolute' }));
    iso = (await infos(v))[0].representations[0].isoValue!;
    expect(iso).toMatchObject({ kind: 'absolute', absolute: 0.5 });
    // The relative field of a kind that was replaced does not stay behind in the state.
    const params = v.plugin.managers.volume.hierarchy.current.volumes[0].representations[0].cell.transform
      .params as any;
    expect(params.type.params.isoValue).toEqual({ kind: 'absolute', absoluteValue: 0.5 });

    expect(((await v.api.volume.isovalue({ value: 0.5, kind: 'absolute' })) as any).unchanged).toBe(true);
  });

  it('refuses an iso value outside the map and a representation that has none, naming what is valid', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    const out = await v.api.volume.isovalue({ value: 5, kind: 'absolute' });
    expect(out).toMatchObject({ ok: false, code: 'invalid' });
    expect((out as any).hint).toMatch(/range from 0\.0+\d* to 0\.9/);

    const missing = await v.api.volume.isovalue({ value: 1, representation: 'slice' });
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
    expect((missing as any).hint).toMatch(/isosurface/);
  });

  it('adds a representation, updates it when it is there, and leaves it when nothing changes', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    expect(reps(v)).toEqual(['isosurface']);

    expectOk(await v.api.volume.representation({ type: 'slice' }));
    expect(reps(v).sort()).toEqual(['isosurface', 'slice']);

    const same = await v.api.volume.representation({ type: 'slice' });
    expectOk(same);
    expect(same.unchanged).toBe(true);
    expect(reps(v)).toHaveLength(2);

    expectOk(await v.api.volume.representation({ type: 'isosurface', params: { alpha: 0.5 } }));
    const iso = v.plugin.managers.volume.hierarchy.current.volumes[0].representations.find(
      (r) => (r.cell.transform.params as any).type.name === 'isosurface',
    )!;
    expect((iso.cell.transform.params as any).type.params.alpha).toBe(0.5);
    // The rest of its params stay.
    expect((iso.cell.transform.params as any).type.params.isoValue.kind).toBe('relative');
    expect(((await v.api.volume.representation({ type: 'isosurface', params: { alpha: 0.5 } })) as any).unchanged).toBe(
      true,
    );

    // An iso value given as a whole replaces the old one.
    expectOk(
      await v.api.volume.representation({
        type: 'isosurface',
        params: { isoValue: { kind: 'absolute', absoluteValue: 0.4 } },
      }),
    );
    expect((await infos(v))[0].representations.find((r) => r.type === 'isosurface')?.isoValue).toMatchObject({
      kind: 'absolute',
      absolute: 0.4,
    });

    const bad = await v.api.volume.representation({ type: 'blob' });
    expect(bad).toMatchObject({ ok: false });
    expect((bad as any).hint).toMatch(/isosurface/);
    const unknown = await v.api.volume.representation({ type: 'isosurface', params: { nope: 1 } });
    expect(unknown).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('updates a representation from the inspector form', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    const rep = v.plugin.managers.volume.hierarchy.current.volumes[0].representations[0];
    const ref = rep.cell.transform.ref;
    const current = rep.cell.transform.params as any;
    // The form gives all of the values.
    const params = { type: { name: 'isosurface', params: { ...current.type.params, alpha: 0.3 } } };
    expectOk(await v.api.volume.update({ representation: ref, params }));
    expect((rep.cell.transform.params as any).type.params.alpha).toBe(0.3);
    const again = await v.api.volume.update({ representation: ref, params });
    expectOk(again);
    expect(again.unchanged).toBe(true);
    expect(await v.api.volume.update({ representation: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('fails cleanly where a representation cannot be made (direct volume needs a canvas), leaving the scene as it was', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    const r = await v.api.volume.representation({ type: 'direct-volume' });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(reps(v)).toEqual(['isosurface']);
    expect([...v.plugin.state.data.cells.values()].filter((c) => c.status === 'error')).toHaveLength(0);
  });

  it('shows and hides a volume or one of its representations', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    const hidden = () => v.plugin.managers.volume.hierarchy.current.volumes[0].cell.state.isHidden;

    expectOk(await v.api.volume.visibility({ visible: false }));
    expect(hidden()).toBe(true);
    expect((await infos(v))[0].hidden).toBe(true);
    expect(((await v.api.volume.visibility({ visible: false })) as any).unchanged).toBe(true);
    expectOk(await v.api.volume.visibility({ visible: true }));
    expect(hidden()).toBeFalsy();

    expectOk(await v.api.volume.visibility({ visible: false, representation: 'isosurface' }));
    expect((await infos(v))[0].representations[0].hidden).toBe(true);
    expect(hidden()).toBeFalsy();
    expect(await v.api.volume.visibility({ visible: false, representation: 'slice' })).toMatchObject({
      ok: false,
      code: 'not-found',
    });
  });

  it('removes a representation or the volume, and undoes', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile()] });
    expectOk(await v.api.volume.representation({ type: 'slice' }));
    expect(reps(v)).toHaveLength(2);

    expectOk(await v.api.volume.remove({ representation: 'slice' }));
    expect(reps(v)).toEqual(['isosurface']);
    expectOk(await v.api.history.undo());
    expect(reps(v)).toHaveLength(2);

    const removed = await v.api.volume.remove({});
    expectOk(removed);
    expect(removed.did).toMatch(/Removed test\.ccp4/);
    expect(await infos(v)).toHaveLength(0);
    // The data and format it was read from go too.
    expect(v.plugin.state.data.cells.size).toBe(1);
    expect(await v.api.volume.remove({})).toMatchObject({ ok: false, code: 'not-found' });

    expectOk(await v.api.history.undo());
    expect((await infos(v)).map((x) => x.label)).toEqual(['test.ccp4']);
    expect(reps(v)).toHaveLength(2);
  });

  it('asks which volume when there are several, and takes a label or a state ref', async () => {
    const v = await headlessViewer();
    await v.api.load.files({ files: [mapFile('a.ccp4')] });
    await v.api.load.files({ files: [mapFile('b.ccp4')] });

    const ambiguous = await v.api.volume.isovalue({ value: 1 });
    expect(ambiguous).toMatchObject({ ok: false, code: 'ambiguous' });
    expect((ambiguous as any).hint).toMatch(/a\.ccp4, b\.ccp4/);

    expectOk(await v.api.volume.isovalue({ value: 1, volume: 'B.CCP4' }));
    const [a, b] = await infos(v);
    expect(a.representations[0].isoValue?.relative).toBeCloseTo(2);
    expect(b.representations[0].isoValue?.relative).toBeCloseTo(1);

    expectOk(await v.api.volume.isovalue({ value: 3, volume: a.ref }));
    expect((await infos(v))[0].representations[0].isoValue?.relative).toBeCloseTo(3);
    expect(await v.api.volume.isovalue({ value: 1, volume: 'c.ccp4' })).toMatchObject({
      ok: false,
      code: 'not-found',
    });
    // A representation named by its state ref needs no volume.
    expectOk(await v.api.volume.isovalue({ value: 1.25, representation: b.representations[0].ref }));
    expect((await infos(v))[1].representations[0].isoValue?.relative).toBeCloseTo(1.25);
  });

  it('streams density only for a structure, and is described as a network operation', async () => {
    const v = await headlessViewer();
    expect(await v.api.volume.streaming({ enabled: true })).toMatchObject({ ok: false, code: 'not-found' });
    const d = v.api.describe('volume.streaming');
    expect(d?.definition.approval).toBe('network');
    expect(d?.definition.modules).toEqual(['network']);
    expect(d?.params.enabled).toBeDefined();

    // Off where it was never on is already so, without touching the network.
    await v.api.load.files({ files: [cbs()] });
    const off = await v.api.volume.streaming({ enabled: false });
    expectOk(off);
    expect(off.unchanged).toBe(true);
  });
});
