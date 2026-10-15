/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { exampleFile, expectOk, headlessViewer } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { LoadersFeature, loadsFromUrl } from '../index.js';

// parity.md §1 (phase 6): `load.data` (a file's contents as text, base64 or bytes, with a format) and the structure
// options of `load.url` and `load.data` (every model, an assembly; v5's `loadAllModelsOrAssemblyFromUrl`), with the
// URL's `all-models` and `assembly`. No network: `load.url` reads a stubbed `fetch`.

afterEach(() => {
  vi.unstubAllGlobals();
});

async function viewer() {
  const v = await headlessViewer({ extensions: [LoadersFeature] });
  v.plugin.state.setSnapshotParams({ camera: false });
  return v;
}
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures;
const cells = (v: PrismViewer) => v.plugin.state.data.cells.size;
const text = async (name: string) => exampleFile(name).text();
const bytes = async (name: string) => new Uint8Array(await exampleFile(name).arrayBuffer());
const typeOf = (v: PrismViewer) =>
  (structures(v)[0].cell.transform.params as { type: { name: string; params?: { id?: string } } }).type;

describe('load.data', () => {
  it('loads mmCIF text as one step, names it, and undoes', async () => {
    const v = await viewer();
    const empty = cells(v);
    const r = await v.api.load.data!({ data: await text('1cbs_updated.cif'), format: 'mmcif' });
    expectOk(r);
    expect(r.did).toBe('Loaded 1CBS.');
    expect(r.value?.loaded).toEqual(['1CBS']);
    expect(structures(v)).toHaveLength(1);
    expect(steps(v)).toBe(1);
    expectOk(await v.api.history.undo());
    expect(cells(v)).toBe(empty);
  });

  it('loads binary data from base64 and from bytes, and reads the format from the label', async () => {
    const v = await viewer();
    const bcif = await bytes('1cbs_full.bcif');
    const base64 = `base64,${Buffer.from(bcif).toString('base64')}`;
    expectOk(await v.api.load.data!({ data: base64, format: 'mmcif' }));
    expectOk(await v.api.load.data!({ bytes: bcif, format: 'mmcif' }));
    expectOk(await v.api.load.data!({ data: await text('adrenalin.mol2'), format: 'auto', label: 'adrenalin.mol2' }));
    expect(structures(v)).toHaveLength(3);
    expect(steps(v)).toBe(3);
  });

  it('refuses what it cannot load, and a failed load leaves nothing', async () => {
    const v = await viewer();
    const empty = cells(v);
    expect(await v.api.load.data!({ data: '' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.data!({ data: 'x', format: 'nonsense' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.data!({ data: 'x', format: 'auto' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.data!({ data: 'x', bytes: new Uint8Array([1]) })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    const bad = await v.api.load.data!({ data: 'not a structure', format: 'mmcif' });
    expect(bad).toMatchObject({ ok: false, code: 'failed' });
    expect(cells(v)).toBe(empty);
    expect(steps(v)).toBe(0);
  });

  it('loads every model, or one, of an ensemble', async () => {
    const v = await viewer();
    const data = await text('1grm_updated.cif');
    expectOk(await v.api.load.data!({ data, format: 'mmcif' }));
    expect(structures(v)).toHaveLength(1);
    await v.api.history.undo();
    expectOk(await v.api.load.data!({ data, format: 'mmcif', models: 'all' }));
    expect(structures(v).length).toBeGreaterThan(1);
  });

  it('builds an assembly, and refuses structure options for a volume', async () => {
    const v = await viewer();
    expectOk(await v.api.load.data!({ data: await text('1cbs_updated.cif'), assembly: '1' }));
    expect(typeOf(v)).toMatchObject({ name: 'assembly', params: { id: '1' } });
    const volume = await v.api.load.data!({ data: 'x', format: 'ccp4', models: 'all' });
    expect(volume).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('has a command', async () => {
    const v = await viewer();
    expect(v.managers.commands.get('load.data')).toBeDefined();
    expect(v.api.load.data?.definition.approval).toBeUndefined();
    expect(v.api.load.data?.definition.modules).toBeUndefined();
  });
});

describe('load.url: every model, an assembly', () => {
  it('loads every model of a file at an address', async () => {
    const data = await bytes('1grm_updated.cif');
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response(data as BodyInit, { status: 200 })),
    );
    const v = await viewer();
    expectOk(await v.api.load.url!({ url: 'https://x.org/1grm.cif', format: 'mmcif', models: 'all' }));
    expect(structures(v).length).toBeGreaterThan(1);
    expect(steps(v)).toBe(1);
  });
});

describe('the URL: all-models and assembly', () => {
  it('passes them to the structure loads only', () => {
    expect(loadsFromUrl('?structure-url=https://x.org/a.cif&all-models=1&assembly=2')[0].values).toMatchObject({
      models: 'all',
      assembly: '2',
    });
    expect(loadsFromUrl('?url=https://x.org/a.cif&all-models=true')[0].values).toMatchObject({ models: 'all' });
    expect(loadsFromUrl('?pdb=1cbs&assembly=1')).toEqual([
      { operation: 'load.entry', values: { source: 'pdb', id: '1cbs', assembly: '1' } },
    ]);
    expect(loadsFromUrl('?emdb=EMD-1832&assembly=1')[0].values).toEqual({ source: 'emdb', id: 'EMD-1832' });
    // Alone they load nothing.
    expect(loadsFromUrl('?assembly=1&all-models=1')).toEqual([]);
  });
});
