/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { mapFile } from '../../../api/_test/volume-fixtures.js';
import { exampleFile, expectOk, headlessViewer, recordRuns, until, withExtensions } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { LoadersFeature } from '../index.js';
import { sessionTypeOf } from '../../../api/session.js';

// The `loaders` feature (viewer-spec.md §3.1): its operations over a stubbed `fetch` that serves the repository's
// example files by address (nothing reaches the network). The plugin downloads through `fetch` in Node.

type Served = Uint8Array | string;

/** Answers each address from `routes` (by substring), else 404. Returns the stub, to count and read the requests. */
function serve(routes: Record<string, Served | Promise<Served>>) {
  const stub = vi.fn(async (input: string | URL | Request) => {
    const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
    const hit = Object.keys(routes).find((k) => url.includes(k));
    if (!hit) return new Response('not found', { status: 404 });
    return new Response((await routes[hit]) as BodyInit, { status: 200 });
  });
  vi.stubGlobal('fetch', stub);
  return stub;
}

// A MolViewSpec of one structure, a cartoon (the examples' labels draw text on a canvas, which Node has not).
const simpleMvs = JSON.stringify({
  metadata: { version: '1', timestamp: '2023-11-24T10:38:17.483Z' },
  root: {
    kind: 'root',
    children: [
      {
        kind: 'download',
        params: { url: 'https://www.ebi.ac.uk/pdbe/entry-files/1cbs.bcif' },
        children: [
          {
            kind: 'parse',
            params: { format: 'bcif' },
            children: [
              {
                kind: 'structure',
                params: { type: 'model' },
                children: [
                  {
                    kind: 'component',
                    params: { selector: 'polymer' },
                    children: [{ kind: 'representation', params: { type: 'cartoon' } }],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
  },
});

const bytes = async (name: string) => new Uint8Array(await exampleFile(name).arrayBuffer());

afterEach(() => {
  vi.unstubAllGlobals();
});

// No canvas here: snapshots leave the camera out, as a viewer without one must.
async function viewer(spec = {}) {
  const v = await headlessViewer({ extensions: [LoadersFeature], spec });
  v.plugin.state.setSnapshotParams({ camera: false });
  return v;
}
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures;
const cells = (v: PrismViewer) => v.plugin.state.data.cells.size;

describe('loaders: the feature', () => {
  it('adds the operations to the load area, and the commands', async () => {
    const v = await viewer();
    for (const name of ['entry', 'url', 'data', 'mvs', 'svs', 'snapshot', 'trajectory', 'volume'] as const) {
      expect(v.api.load[name], name).toBeTypeOf('function');
      // Replacing the scene is destructive, which wins over network (api/area.ts).
      const replaces = name === 'mvs' || name === 'snapshot';
      expect(v.api.load[name]?.definition.approval, name).toBe(replaces ? 'destructive' : 'network');
      // The loaders that take url, data or a file check the input's module per call, not as the operation's.
      const perCall = ['data', 'mvs', 'svs', 'snapshot', 'volume'].includes(name);
      expect(v.api.load[name]?.definition.modules, name).toEqual(perCall ? undefined : ['network']);
    }
    for (const id of ['load.entry', 'load.url', 'load.mvs', 'load.svs', 'load.snapshot'])
      expect(v.managers.commands.get(id)).toBeDefined();
  });

  it('brings no views: Home’s Download Structure is the form for an id', async () => {
    expect(LoadersFeature.ui).toBeUndefined();
  });

  it('keeps the core load operations', async () => {
    const v = await viewer();
    expect(v.api.load.files).toBeTypeOf('function');
    expect(v.api.load.action).toBeTypeOf('function');
  });
});

describe('load.entry', () => {
  it('loads a PDB entry as one step, and names it', async () => {
    const stub = serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    const r = await v.api.load.entry!({ source: 'pdb', id: '1cbs' });
    expectOk(r);
    expect(r.did).toBe('Loaded 1CBS.');
    expect(r.value?.loaded).toEqual(['1CBS']);
    expect(structures(v)).toHaveLength(1);
    expect(steps(v)).toBe(1);
    expect(stub).toHaveBeenCalledTimes(1);
    expect(String(stub.mock.calls[0][0])).toContain('pdbe');
  });

  it('is already so the second time: nothing changes, no step, no request', async () => {
    const stub = serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const r = await v.api.load.entry!({ source: 'pdb', id: '1CBS' });
    expect(r).toMatchObject({ ok: true, unchanged: true, did: '1CBS is already loaded.' });
    expect(structures(v)).toHaveLength(1);
    expect(steps(v)).toBe(1);
    expect(stub).toHaveBeenCalledTimes(1);
  });

  it('is already so when asked twice at once: the second, queued behind the first, finds it loaded', async () => {
    const stub = serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    const [a, b] = await Promise.all([
      v.api.load.entry!({ source: 'pdb', id: '1cbs' }),
      v.api.load.entry!({ source: 'pdb', id: '1cbs' }),
    ]);
    expectOk(a);
    expect(b).toMatchObject({ ok: true, unchanged: true, did: '1CBS is already loaded.' });
    expect(structures(v)).toHaveLength(1);
    expect(steps(v)).toBe(1);
    expect(stub).toHaveBeenCalledTimes(1);
  });

  it('loads another with copy', async () => {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const r = await v.api.load.entry!({ source: 'pdb', id: '1cbs' }, { copy: true });
    expectOk(r);
    expect(r.unchanged).toBeUndefined();
    expect(structures(v)).toHaveLength(2);
    expect(steps(v)).toBe(2);
  });

  it('loads again after the structure was removed', async () => {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    await v.api.history.undo();
    expect(structures(v)).toHaveLength(0);
    const r = await v.api.load.entry!({ source: 'pdb', id: '1cbs' });
    expectOk(r);
    expect(r.unchanged).toBeUndefined();
    expect(structures(v)).toHaveLength(1);
  });

  it('undoes as one step', async () => {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    const empty = cells(v);
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    expectOk(await v.api.history.undo());
    expect(cells(v)).toBe(empty);
  });

  it('fails on a 404 and leaves nothing', async () => {
    serve({});
    const v = await viewer();
    const empty = cells(v);
    const r = await v.api.load.entry!({ source: 'pdb', id: '1xyz' });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    if (!r.ok) expect(r.message).toContain('404');
    expect(cells(v)).toBe(empty);
    expect(structures(v)).toHaveLength(0);
    expect(steps(v)).toBe(0);
  });

  it('reads the server from the viewer options: PDBj serves mmCIF text', async () => {
    const stub = serve({ '1cbs.cif': await bytes('1cbs_updated.cif') });
    const v = await viewer({ options: { data: { pdbProvider: 'pdbj' } } });
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    expect(String(stub.mock.calls[0][0])).toContain('pdbj');
    expect(structures(v)).toHaveLength(1);
  });

  it('takes several ids, and loads only the ones that are not there', async () => {
    const stub = serve({ '1cbs.bcif': await bytes('1cbs_full.bcif'), '1crn.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const r = await v.api.load.entry!({ source: 'pdb', id: '1cbs, 1crn' });
    expectOk(r);
    expect(stub).toHaveBeenCalledTimes(2);
    expect(String(stub.mock.calls[1][0])).toContain('1crn');
    expect(structures(v)).toHaveLength(2);
  });

  it('accepts pdb-dev as pdb-ihm', async () => {
    const stub = serve({ 'pdb-ihm.org/bcif/': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb-dev', id: '8zzc' }));
    expect(String(stub.mock.calls[0][0])).toContain('pdb-ihm.org/bcif/8zzc.bcif');
    const again = await v.api.load.entry!({ source: 'pdb-ihm', id: '8ZZC' });
    expect(again).toMatchObject({ ok: true, unchanged: true });
  });

  it('loads AlphaFold DB through its API, with the pLDDT preset', async () => {
    const stub = serve({
      'alphafold.ebi.ac.uk/api/prediction/Q8W3K0': JSON.stringify([{ bcifUrl: 'https://af.test/AF-Q8W3K0-F1.bcif' }]),
      'af.test/': await bytes('1cbs_full.bcif'),
    });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'afdb', id: 'q8w3k0' }));
    expect(stub).toHaveBeenCalledTimes(2);
    expect(structures(v)).toHaveLength(1);
    expect(await v.api.load.entry!({ source: 'afdb', id: 'q8w3k0' })).toMatchObject({ unchanged: true });
  });

  it('loads Model Archive', async () => {
    const stub = serve({ 'modelarchive.org/doi/10.5452/ma-bak-cepc-0003.cif': await bytes('1cbs_updated.cif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'model-archive', id: 'ma-bak-cepc-0003' }));
    expect(stub).toHaveBeenCalledTimes(1);
    expect(structures(v)).toHaveLength(1);
  });

  it('checks its values: no id, an unknown preset', async () => {
    serve({});
    const v = await viewer();
    expect(await v.api.load.entry!({ source: 'pdb', id: '' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.entry!({ source: 'pdb', id: '1cbs', representation: 'nonsense' })).toMatchObject({
      ok: false,
      code: 'not-found',
    });
  });

  it('builds the assembly asked for', async () => {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs', assembly: '1' }));
    const params = [...v.plugin.state.data.cells.values()].find((c) =>
      c.transform.transformer.id.includes('structure-from-model'),
    )?.transform.params as { type?: { name: string } } | undefined;
    expect(params?.type?.name).toBe('assembly');
  });
});

describe('load.url', () => {
  it('loads a structure from an address, and is already so the second time', async () => {
    serve({ 'files.test/1cbs.cif': await bytes('1cbs_updated.cif') });
    const v = await viewer();
    const r = await v.api.load.url!({ url: 'https://files.test/1cbs.cif', format: 'mmcif' });
    expectOk(r);
    expect(structures(v)).toHaveLength(1);
    expect(r.value?.loaded).toHaveLength(1);
    expect(steps(v)).toBe(1);
    const again = await v.api.load.url!({ url: 'https://files.test/1cbs.cif', format: 'mmcif' });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(steps(v)).toBe(1);
    expectOk(await v.api.load.url!({ url: 'https://files.test/1cbs.cif', format: 'mmcif' }, { copy: true }));
    expect(structures(v)).toHaveLength(2);
  });

  it('is already so when asked twice at once', async () => {
    serve({ 'files.test/1cbs.cif': await bytes('1cbs_updated.cif') });
    const v = await viewer();
    const url = 'https://files.test/1cbs.cif';
    const [a, b] = await Promise.all([
      v.api.load.url!({ url, format: 'mmcif' }),
      v.api.load.url!({ url, format: 'mmcif' }),
    ]);
    expectOk(a);
    expect(b).toMatchObject({ ok: true, unchanged: true, did: `${url} is already loaded.` });
    expect(structures(v)).toHaveLength(1);
  });

  it('reads the format and whether it is binary from the address with auto', async () => {
    serve({ 'files.test/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.url!({ url: 'https://files.test/1cbs.bcif', format: 'auto' }));
    expect(structures(v)).toHaveLength(1);
  });

  it('takes bcif as the name of the binary mmCIF format', async () => {
    serve({ 'files.test/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.url!({ url: 'https://files.test/1cbs.bcif', format: 'bcif' }));
    expect(structures(v)).toHaveLength(1);
  });

  it('fails on a 404 and on a file that does not parse, leaving nothing', async () => {
    serve({ 'files.test/bad.cif': 'this is not mmCIF' });
    const v = await viewer();
    const empty = cells(v);
    expect(await v.api.load.url!({ url: 'https://files.test/missing.cif', format: 'mmcif' })).toMatchObject({
      ok: false,
      code: 'failed',
    });
    expect(await v.api.load.url!({ url: 'https://files.test/bad.cif', format: 'mmcif' })).toMatchObject({
      ok: false,
      code: 'failed',
    });
    expect(cells(v)).toBe(empty);
    expect(steps(v)).toBe(0);
  });

  it('checks its values: an empty address, an unknown format', async () => {
    serve({});
    const v = await viewer();
    expect(await v.api.load.url!({ url: ' ' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.load.url!({ url: 'https://files.test/a.nope', format: 'nope' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });
});

describe('load.mvs', () => {
  it('loads a MolViewSpec from an address', async () => {
    const stub = serve({
      'files.test/1cbs.mvsj': simpleMvs,
      'entry-files/1cbs.bcif': await bytes('1cbs_full.bcif'),
    });
    const v = await viewer();
    const r = await v.api.load.mvs!({ url: 'https://files.test/1cbs.mvsj', keepCamera: true });
    expectOk(r);
    expect(structures(v).length).toBeGreaterThan(0);
    expect(stub.mock.calls.length).toBeGreaterThanOrEqual(2);
  });

  it('loads a MolViewSpec given as data', async () => {
    serve({ 'entry-files/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.mvs!({ data: simpleMvs, keepCamera: true }));
    expect(structures(v).length).toBeGreaterThan(0);
  });

  it('loads data without the network module, and needs it for an address', async () => {
    serve({ 'entry-files/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const spec = DefaultPrismViewerSpec();
    const v = await createPrismViewerHeadless({
      spec: { ...spec, extensions: withExtensions(spec.extensions ?? [], [LoadersFeature]) },
      modules: ['files'],
    });
    try {
      expect(v.api.describe('load.mvs')?.available).toBe(true);
      expect(await v.api.load.mvs!({ url: 'https://files.test/1cbs.mvsj' })).toMatchObject({
        ok: false,
        code: 'unavailable',
      });
      expectOk(await v.api.load.mvs!({ data: simpleMvs, keepCamera: true }));
      expect(structures(v).length).toBeGreaterThan(0);
    } finally {
      await v.dispose();
    }
  });

  it('fails on a 404, on invalid data, and on both a url and data', async () => {
    serve({});
    const v = await viewer();
    const empty = cells(v);
    expect(await v.api.load.mvs!({ url: 'https://files.test/none.mvsj' })).toMatchObject({ ok: false, code: 'failed' });
    expect(await v.api.load.mvs!({ data: '{"not": "mvs"}' })).toMatchObject({ ok: false, code: 'failed' });
    expect(await v.api.load.mvs!({ url: 'https://a.test/x.mvsj', data: '{}' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(await v.api.load.mvs!({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(cells(v)).toBe(empty);
  });
});

describe('load.snapshot', () => {
  async function savedState(): Promise<string> {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const a = await viewer();
    expectOk(await a.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    return (await a.plugin.managers.snapshot.serialize({ type: 'molj' })).text();
  }

  it('opens a state from an address', async () => {
    const molj = await savedState();
    serve({ 'files.test/state.molj': molj, '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    const r = await v.api.load.snapshot!({ url: 'https://files.test/state.molj', type: 'molj' });
    expectOk(r);
    await until(() => structures(v).length === 1);
  });

  it('follows the address’s extension when no type is given, and obeys one that is', async () => {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const a = await viewer();
    expectOk(await a.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const molx = new Uint8Array(await (await a.plugin.managers.snapshot.serialize({ type: 'molx' })).arrayBuffer());
    const molj = await (await a.plugin.managers.snapshot.serialize({ type: 'molj' })).text();

    serve({
      'files.test/state.molx': molx,
      'files.test/state.molj': molj,
      'files.test/download': molx,
      '1cbs.bcif': await bytes('1cbs_full.bcif'),
    });
    // A .molx address, no type: read as a session (as a state it would fail to parse).
    const v = await viewer();
    expectOk(await v.api.load.snapshot!({ url: 'https://files.test/state.molx?download=1' }));
    await until(() => structures(v).length === 1);
    // A .molj address, no type: a state.
    const w = await viewer();
    expectOk(await w.api.load.snapshot!({ url: 'https://files.test/state.molj' }));
    await until(() => structures(w).length === 1);
    // An address without an extension: the type says.
    const x = await viewer();
    expectOk(await x.api.load.snapshot!({ url: 'https://files.test/download', type: 'molx' }));
    await until(() => structures(x).length === 1);
  });

  it('reads the kind from an address’s path', () => {
    expect(sessionTypeOf('https://a.test/s.molx')).toBe('molx');
    expect(sessionTypeOf('https://a.test/S.MOLJ?x=1#y')).toBe('molj');
    expect(sessionTypeOf('https://a.test/download?name=s.molx')).toBeUndefined();
  });

  it('opens a remote state by id from the state server', async () => {
    const molj = JSON.parse(await savedState());
    const stub = serve({
      'state.test/get/abc': JSON.stringify({ data: molj }),
      '1cbs.bcif': await bytes('1cbs_full.bcif'),
    });
    const v = await viewer({ options: { data: { stateServer: 'https://state.test' } } });
    expectOk(await v.api.load.snapshot!({ id: 'abc' }));
    expect(String(stub.mock.calls[0][0])).toBe('https://state.test/get/abc');
    await until(() => structures(v).length === 1);
  });

  it('loads the lazy extensions that own the state’s transformers before applying it, by id too', async () => {
    const molj = JSON.parse(await savedState());
    serve({ 'state.test/get/abc': JSON.stringify({ data: molj }), '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer({ options: { data: { stateServer: 'https://state.test' } } });
    const ensure = vi.spyOn(v.extensions, 'ensureFor');
    expectOk(await v.api.load.snapshot!({ id: 'abc' }));
    expect(ensure).toHaveBeenCalledTimes(1);
    expect(ensure.mock.calls[0][0].transformers).toContain('ms-plugin.download');
    await until(() => structures(v).length === 1);
  });

  it('fails on a state that is not one, and on a 404, by id', async () => {
    serve({ 'state.test/get/bad': JSON.stringify({ data: { nothing: true } }) });
    const v = await viewer({ options: { data: { stateServer: 'https://state.test' } } });
    expect(await v.api.load.snapshot!({ id: 'bad' })).toMatchObject({ ok: false, code: 'failed' });
    expect(await v.api.load.snapshot!({ id: 'missing' })).toMatchObject({ ok: false, code: 'failed' });
  });

  it('fails on a 404 and on a file that is not a state', async () => {
    serve({ 'files.test/bad.molj': 'not json at all' });
    const v = await viewer();
    expect(await v.api.load.snapshot!({ url: 'https://files.test/none.molj' })).toMatchObject({
      ok: false,
      code: 'failed',
    });
    expect(await v.api.load.snapshot!({ url: 'https://files.test/bad.molj' })).toMatchObject({
      ok: false,
      code: 'failed',
    });
    expect(await v.api.load.snapshot!({})).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('load.trajectory', () => {
  it('loads a model with coordinates', async () => {
    serve({
      'files.test/protein.pdb': await bytes('trajectory/protein.pdb'),
      'files.test/protein.dcd': await bytes('trajectory/protein.dcd'),
    });
    const v = await viewer();
    const r = await v.api.load.trajectory!({
      model: 'https://files.test/protein.pdb',
      modelFormat: 'pdb',
      coordinates: 'https://files.test/protein.dcd',
      coordinatesFormat: 'dcd',
    });
    expectOk(r);
    expect(structures(v)).toHaveLength(1);
    expect(steps(v)).toBe(1);
    expectOk(await v.api.history.undo());
    expect(structures(v)).toHaveLength(0);
  });

  it('fails when the coordinates are missing, leaving nothing', async () => {
    serve({ 'files.test/protein.pdb': await bytes('trajectory/protein.pdb') });
    const v = await viewer();
    const empty = cells(v);
    const r = await v.api.load.trajectory!({
      model: 'https://files.test/protein.pdb',
      modelFormat: 'pdb',
      coordinates: 'https://files.test/none.dcd',
      coordinatesFormat: 'dcd',
    });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(cells(v)).toBe(empty);
    expect(await v.api.load.trajectory!({ model: 'https://files.test/protein.pdb' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });
});

describe('load.volume', () => {
  it('loads a map with an isosurface per iso value', async () => {
    serve({ 'files.test/map.ccp4': new Uint8Array(await mapFile().arrayBuffer()) });
    const v = await viewer();
    const r = await v.api.load.volume!({ url: 'https://files.test/map.ccp4', format: 'ccp4', binary: true });
    expectOk(r);
    expect(v.plugin.managers.volume.hierarchy.current.volumes).toHaveLength(1);
    expect(v.plugin.managers.volume.hierarchy.current.volumes[0].representations).toHaveLength(1);
    expect(steps(v)).toBe(1);
    expect(
      await v.api.load.volume!({ url: 'https://files.test/map.ccp4', format: 'ccp4', binary: true }),
    ).toMatchObject({
      unchanged: true,
    });
  });

  it('is already so when asked twice at once', async () => {
    serve({ 'files.test/map.ccp4': new Uint8Array(await mapFile().arrayBuffer()) });
    const v = await viewer();
    const params = { url: 'https://files.test/map.ccp4', format: 'ccp4', binary: true };
    const [a, b] = await Promise.all([v.api.load.volume!(params), v.api.load.volume!(params)]);
    expectOk(a);
    expect(b).toMatchObject({ ok: true, unchanged: true });
    expect(v.plugin.managers.volume.hierarchy.current.volumes).toHaveLength(1);
  });

  it('fails on a 404 and on a format that is unknown', async () => {
    serve({});
    const v = await viewer();
    const empty = cells(v);
    expect(await v.api.load.volume!({ url: 'https://files.test/none.ccp4', format: 'ccp4' })).toMatchObject({
      ok: false,
      code: 'failed',
    });
    expect(await v.api.load.volume!({ url: 'https://files.test/a.map', format: 'nope' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(cells(v)).toBe(empty);
  });
});

describe('the runs', () => {
  it('record the loads with their origin', async () => {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    const runs = recordRuns(v);
    await v.api.as({ kind: 'url' }).run('load.entry', { source: 'pdb', id: '1cbs' });
    await until(() => runs.length > 0);
    expect(runs[0]).toMatchObject({ operation: 'load.entry', origin: { kind: 'url' } });
  });
});
