/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { zipFiles } from '@molstar/seq-core/format';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { type ApiContext, ApiResult } from '../../../api/area.js';
import { mapFile } from '../../../api/_test/volume-fixtures.js';
import { cbs, exampleFile, expectOk, headlessViewer } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { chooseFormat, type FormatScope, formatNameOf } from '../../../api/formats.js';
import { LoadersFeature } from '../index.js';

// One loading vocabulary (operations.md §1, "Inputs"; the branch review's A6): `format` is 'auto' by default, told by
// Mol*'s format registry from the name's extension (else the data), binary told from the extension or the format; each
// loader takes exactly one of url, data (text or bytes) and file; `session.open` takes them or a state server's id, and
// `load.snapshot` is the same operation. Nothing reaches the network: `fetch` is stubbed.

type Served = Uint8Array | string;

/** Answers each address from `routes` (by substring), else 404. */
function serve(routes: Record<string, Served>) {
  const stub = vi.fn(async (input: string | URL | Request) => {
    const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
    const hit = Object.keys(routes).find((k) => url.includes(k));
    if (!hit) return new Response('not found', { status: 404 });
    return new Response(routes[hit] as BodyInit, { status: 200 });
  });
  vi.stubGlobal('fetch', stub);
  return stub;
}

afterEach(() => {
  vi.unstubAllGlobals();
});

async function viewer(spec = {}) {
  const v = await headlessViewer({ extensions: [LoadersFeature], spec });
  v.plugin.state.setSnapshotParams({ camera: false });
  return v;
}
const bytes = async (name: string) => new Uint8Array(await exampleFile(name).arrayBuffer());
const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures;
const volumes = (v: PrismViewer) => v.plugin.managers.volume.hierarchy.current.volumes;
const mapBytes = async () => new Uint8Array(await mapFile().arrayBuffer());

// A MolViewSpec of one structure, a cartoon (no labels: they draw text on a canvas, which Node has not).
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
const mvsx = () => zipFiles({ 'index.mvsj': new TextEncoder().encode(simpleMvs) });

describe('the format, told', () => {
  // Scoped, as each loader is: a `.cif` may also be a density server's or a particle list's until the data says.
  const Structures: FormatScope = { categories: ['Trajectory'], what: 'a structure' };
  const Volumes: FormatScope = { categories: ['Volume'], what: 'a volume' };
  const Coordinates: FormatScope = { categories: ['Coordinates'], what: 'a coordinates' };

  async function picker() {
    const v = await viewer();
    const ctx = { plugin: v.plugin, extensions: v.extensions } as unknown as ApiContext;
    /** The format and binary chosen, or the failure. */
    type Picked = { ok?: false; code?: string; message?: string; hint?: string; format?: string; binary?: boolean };
    return async (name: string, format?: string, scope?: FormatScope, binary?: boolean): Promise<Picked> => {
      const c = await chooseFormat(ctx, { format, name, binary, ...(scope ? { scope } : {}) });
      return ApiResult.isFailure(c) ? c : { format: formatNameOf(c), binary: c.binary };
    };
  }

  it("reads the format and whether it is binary from the extension, by Mol*'s registry", async () => {
    const pick = await picker();
    expect(await pick('1abc.cif', 'auto', Structures)).toEqual({ format: 'mmcif', binary: false });
    expect(await pick('1abc.cif.gz', 'auto', Structures)).toEqual({ format: 'mmcif', binary: false });
    expect(await pick('1abc.bcif', 'auto', Structures)).toEqual({ format: 'mmcif', binary: true });
    expect(await pick('a.pdb', 'auto', Structures)).toEqual({ format: 'pdb', binary: false });
    expect(await pick('a.gro', 'auto', Structures)).toEqual({ format: 'gro', binary: false });
    for (const ext of ['map', 'mrc', 'ccp4']) {
      expect(await pick(`emd.${ext}`, 'auto', Volumes)).toEqual({ format: 'ccp4', binary: true });
    }
    expect(await pick('x.dsn6', 'auto', Volumes)).toEqual({ format: 'dsn6', binary: true });
    expect(await pick('x.cube', 'auto', Volumes)).toEqual({ format: 'cube', binary: false });
    for (const ext of ['dcd', 'xtc', 'trr']) {
      expect(await pick(`md.${ext}`, 'auto', Coordinates)).toEqual({ format: ext, binary: true });
    }
  });

  it('infers binary from a format given by name, and obeys binary when given', async () => {
    const pick = await picker();
    expect(await pick('download', 'ccp4')).toEqual({ format: 'ccp4', binary: true });
    expect(await pick('download', 'mmcif')).toEqual({ format: 'mmcif', binary: false });
    expect(await pick('download', 'bcif')).toEqual({ format: 'mmcif', binary: true });
    expect(await pick('1abc.bcif', 'mmcif')).toEqual({ format: 'mmcif', binary: true });
    expect(await pick('download', 'mmcif', undefined, true)).toEqual({ format: 'mmcif', binary: true });
    expect(await pick('emd.map', 'ccp4', undefined, false)).toEqual({ format: 'ccp4', binary: false });
    // A density server's cell at an address without an extension: BinaryCIF, as v5 read it.
    expect(await pick('https://ds.example/em/emd-1234/cell', 'dscif')).toEqual({ format: 'dscif', binary: true });
    expect(await pick('cell.cif', 'dscif')).toEqual({ format: 'dscif', binary: false });
  });

  it('fails invalid where it cannot tell, naming the formats', async () => {
    const pick = await picker();
    const none = await pick('download');
    expect(none).toMatchObject({ ok: false, code: 'invalid', message: 'Which format is download?' });
    if (!none.ok) {
      expect(none.hint).toMatch(/^Give its format: auto, /);
      expect(none.hint).toContain('mmcif');
      expect(none.hint).toContain('ccp4');
    }
    expect(await pick('a.nope')).toMatchObject({ ok: false, code: 'invalid' });
    // A structure is not a volume: the volume formats only.
    const pdb = await pick('a.pdb', 'auto', Volumes);
    expect(pdb).toMatchObject({ ok: false, code: 'invalid' });
    if (!pdb.ok) {
      expect(pdb.hint).toContain('ccp4');
      expect(pdb.hint).not.toContain('mmcif');
    }
    expect(await pick('a', 'pdb', Volumes)).toMatchObject({ ok: false, code: 'invalid' });
    expect(await pick('a', 'nope')).toMatchObject({ ok: false, code: 'invalid' });
  });
});

describe('load.url and load.data: auto by default', () => {
  it('loads BinaryCIF from a .bcif address with neither format nor binary given', async () => {
    serve({ 'files.test/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.url!({ url: 'https://files.test/1cbs.bcif' }));
    expect(structures(v)).toHaveLength(1);
  });

  it('refuses an address that does not say its format before fetching it', async () => {
    const stub = serve({ 'files.test/download': await bytes('1cbs_updated.cif') });
    const v = await viewer();
    const r = await v.api.load.url!({ url: 'https://files.test/download' });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    expect(stub).not.toHaveBeenCalled();
    expectOk(await v.api.load.url!({ url: 'https://files.test/download', format: 'mmcif' }));
  });

  it('loads data without a name by what it is: mmCIF text, PDB records, BinaryCIF bytes', async () => {
    const v = await viewer();
    expectOk(await v.api.load.data!({ data: await exampleFile('1cbs_updated.cif').text() }));
    expectOk(await v.api.load.data!({ data: await exampleFile('trajectory/protein.pdb').text() }));
    expectOk(await v.api.load.data!({ data: await bytes('1cbs_full.bcif') }));
    expect(structures(v)).toHaveLength(3);
    expect(await v.api.load.data!({ data: 'neither' })).toMatchObject({ ok: false, code: 'invalid' });
  });

  it('loads a file, its format from its name', async () => {
    const v = await viewer();
    const r = await v.api.load.data!({ file: cbs() });
    expectOk(r);
    expect(r.value?.loaded).toEqual(['1CBS']);
    expectOk(await v.api.load.data!({ file: exampleFile('1cbs_full.bcif', 'x.bcif') }));
    expect(structures(v)).toHaveLength(2);
  });

  it('loads from a url too, with load.url’s identity', async () => {
    serve({ 'files.test/1cbs.cif': await bytes('1cbs_updated.cif') });
    const v = await viewer();
    expectOk(await v.api.load.data!({ url: 'https://files.test/1cbs.cif' }));
    expect(await v.api.load.url!({ url: 'https://files.test/1cbs.cif' })).toMatchObject({ ok: true, unchanged: true });
  });
});

describe('one input: url, data or file', () => {
  it('fails invalid with none, and with several', async () => {
    const v = await viewer();
    const file = cbs();
    for (const op of ['load.data', 'load.volume', 'load.mvs', 'load.svs'] as const) {
      expect(await v.api.run(op, {}), op).toMatchObject({
        ok: false,
        code: 'invalid',
        message: 'Load what? Give url, data or file.',
      });
      expect(await v.api.run(op, { url: 'https://x.test/a.cif', file }), op).toMatchObject({
        ok: false,
        code: 'invalid',
        message: 'Give one of url, data or file.',
      });
      expect(await v.api.run(op, { data: 'x', file }), op).toMatchObject({ ok: false, code: 'invalid' });
    }
    expect(await v.api.session.open({})).toMatchObject({
      ok: false,
      code: 'invalid',
      message: 'Open what? Give url, data, file or id.',
    });
    expect(await v.api.session.open({ data: '{}', id: 'abc' })).toMatchObject({
      ok: false,
      code: 'invalid',
      message: 'Give one of url, data, file or id.',
    });
  });

  it('takes data as text or bytes, nothing else', async () => {
    const v = await viewer();
    expect(await v.api.run('load.data', { data: 42 })).toMatchObject({
      ok: false,
      code: 'invalid',
      message: "'data' takes text or bytes.",
    });
  });

  it('checks the module per input: a url needs the network, data does not', async () => {
    const spec = DefaultPrismViewerSpec();
    const v = await createPrismViewerHeadless({
      spec: { ...spec, extensions: [...(spec.extensions ?? []).filter((e) => e.id !== 'loaders'), LoadersFeature] },
      modules: ['files'],
    });
    try {
      for (const op of ['load.data', 'load.volume', 'session.open']) {
        expect(await v.api.run(op, { url: 'https://x.test/a.ccp4' }), op).toMatchObject({
          ok: false,
          code: 'unavailable',
        });
      }
      expectOk(await v.api.load.volume!({ data: await mapBytes() }));
    } finally {
      await v.dispose();
    }
  });
});

describe('load.volume: auto by default', () => {
  it('loads a map from an address by its extension, binary told', async () => {
    serve({ 'files.test/emd.map': await mapBytes() });
    const v = await viewer();
    expectOk(await v.api.load.volume!({ url: 'https://files.test/emd.map' }));
    expect(volumes(v)).toHaveLength(1);
  });

  it('loads a map from bytes (what it is: a CCP4 header) and from a file', async () => {
    const v = await viewer();
    expectOk(await v.api.load.volume!({ data: await mapBytes() }));
    expectOk(await v.api.load.volume!({ file: mapFile('emd.mrc') }));
    expect(volumes(v)).toHaveLength(2);
  });

  it('reads lazily from an address only', async () => {
    const v = await viewer();
    expect(await v.api.load.volume!({ data: await mapBytes(), lazy: true })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });
});

describe('load.mvs: mvsj or mvsx, told', () => {
  it('loads an archive from an .mvsx address, from bytes and from a file, with no format given', async () => {
    const archive = await mvsx();
    serve({ 'files.test/story.mvsx': archive, 'entry-files/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.mvs!({ url: 'https://files.test/story.mvsx', keepCamera: true }));
    expect(structures(v)).toHaveLength(1);
    expectOk(await v.api.load.mvs!({ data: archive, keepCamera: true }));
    expect(structures(v)).toHaveLength(1);
    expectOk(
      await v.api.load.mvs!({ file: new File([archive as Uint8Array<ArrayBuffer>], 'story.mvsx'), keepCamera: true }),
    );
    expect(structures(v)).toHaveLength(1);
  });

  it('loads JSON from a file', async () => {
    serve({ 'entry-files/1cbs.bcif': await bytes('1cbs_full.bcif') });
    const v = await viewer();
    expectOk(await v.api.load.mvs!({ file: new File([simpleMvs], 'story.mvsj'), keepCamera: true }));
    expect(structures(v)).toHaveLength(1);
  });
});

describe('session.open: url, data, file or id; load.snapshot is the same', () => {
  /** A session (.molx) of 1CBS opened from a file, and a state (.molj) of 1CBS downloaded. */
  async function saved() {
    serve({ '1cbs.bcif': await bytes('1cbs_full.bcif') });
    const a = await headlessViewer({ extensions: [LoadersFeature], files: [cbs()] });
    a.plugin.state.setSnapshotParams({ camera: false });
    const molx = new Uint8Array(await (await a.plugin.managers.snapshot.serialize({ type: 'molx' })).arrayBuffer());
    const b = await viewer();
    expectOk(await b.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const molj = await (await b.plugin.managers.snapshot.serialize({ type: 'molj' })).text();
    return { molx, molj };
  }

  it('opens a session from bytes and from a file, its kind told by the content or the name', async () => {
    const { molx } = await saved();
    const v = await viewer();
    const r = await v.api.session.open({ data: molx });
    expectOk(r);
    expect(r.did).toBe('Opened the session.');
    expect(r.value?.loaded).toHaveLength(1);
    expect(structures(v)).toHaveLength(1);
    const w = await viewer();
    expectOk(await w.api.session.open({ file: new File([molx as Uint8Array<ArrayBuffer>], 'mine.molx') }));
    expect(structures(w)).toHaveLength(1);
  });

  it('opens a state from its JSON as data', async () => {
    const { molj } = await saved();
    const v = await viewer();
    expectOk(await v.api.session.open({ data: molj }));
    expect(structures(v)).toHaveLength(1);
  });

  it('opens from an address without an extension (a zip is a session), and asks as destructive', async () => {
    const { molx } = await saved();
    serve({ 'files.test/download': molx });
    const v = await viewer();
    const r = await v.api.session.open({ url: 'https://files.test/download' });
    expectOk(r);
    expect(r.did).toBe('Opened download.');
    expect(structures(v)).toHaveLength(1);
    const def = v.api.session.open.definition;
    expect(def.approval).toBe('destructive');
    expect(def.modules).toBeUndefined();
    expect(def.addresses?.({ url: 'https://files.test/download' } as never)).toEqual(['https://files.test/download']);
  });

  it('opens a state of the state server by id', async () => {
    const { molj } = await saved();
    serve({
      'state.test/get/abc': JSON.stringify({ data: JSON.parse(molj) }),
      '1cbs.bcif': await bytes('1cbs_full.bcif'),
    });
    const v = await viewer({ options: { data: { stateServer: 'https://state.test' } } });
    const r = await v.api.session.open({ id: 'abc' });
    expectOk(r);
    expect(r.did).toBe('Opened state abc.');
    expect(structures(v)).toHaveLength(1);
  });

  it('keeps load.snapshot as the same operation, its params and its run', async () => {
    const { molx } = await saved();
    const v = await viewer();
    const alias = v.api.load.snapshot!.definition;
    expect(alias.params).toBe(v.api.session.open.definition.params);
    expect(alias.run).toBe(v.api.session.open.definition.run);
    expect(alias.approval).toBe('destructive');
    expectOk(await v.api.load.snapshot!({ data: molx }));
    expect(structures(v)).toHaveLength(1);
  });

  it('fails on data that is no session, leaving what was there', async () => {
    const v = await viewer();
    expect(await v.api.session.open({ data: 'not a state' })).toMatchObject({ ok: false, code: 'failed' });
    expect(await v.api.session.open({ data: new Uint8Array([0x50, 0x4b, 1, 2]) })).toMatchObject({
      ok: false,
      code: 'failed',
    });
  });
});
