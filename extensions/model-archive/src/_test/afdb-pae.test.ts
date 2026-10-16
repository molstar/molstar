/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { afdbEntryOf, paeValuesFromAfdb, resolveAfdbPaeUrl } from '../pae/afdb.js';
import { ModelArchiveExtension } from '../prism-viewer.js';

// AlphaFold DB's PAE (parity.md §6, "PAE plot"): which structure is an entry, where its PAE JSON is, and the fetch into
// the state tree on the recorded P04637 model. Every request is a stub: no network.

const Entry = { id: 'AF-P04637-F1', accession: 'P04637', fragment: 1 };
const Api = 'https://alphafold.ebi.ac.uk/api/prediction/P04637';
const Json = 'https://alphafold.ebi.ac.uk/files/AF-P04637-F1-predicted_aligned_error_v6.json';

describe('AlphaFold DB entry detection', () => {
  it('reads the id from any of the words a structure carries, first match wins, with the file version', () => {
    expect(afdbEntryOf([undefined, 'AF-P04637-F1'])).toEqual(Entry);
    expect(afdbEntryOf(['1cbs', 'https://alphafold.ebi.ac.uk/files/af-p04637-f1-model_v4.cif'])).toEqual({
      ...Entry,
      version: 4,
    });
    expect(afdbEntryOf(['AF-A0A0B4J2F0-F12.bcif'])).toEqual({
      id: 'AF-A0A0B4J2F0-F12',
      accession: 'A0A0B4J2F0',
      fragment: 12,
    });
    expect(afdbEntryOf(['1CBS', 'model.cif', undefined, ''])).toBeUndefined();
  });
});

const response = (ok: boolean, body: unknown = null, status = ok ? 200 : 404) =>
  Promise.resolve({ ok, status, json: async () => body }) as never;

describe('resolving the PAE address', () => {
  it("takes the API entry's paeDocUrl for this id and fragment", async () => {
    const fetcher = vi.fn(async (url: string) =>
      url === Api
        ? response(true, [
            { entryId: 'AF-P04637-F2', paeDocUrl: 'https://alphafold.ebi.ac.uk/files/other.json' },
            { entryId: 'AF-P04637-F1', paeDocUrl: Json },
          ])
        : response(false),
    );
    expect(await resolveAfdbPaeUrl(Entry, fetcher)).toBe(Json);
  });

  it('never follows an address outside AlphaFold DB; falls back to the conventional file, probed', async () => {
    const fetcher = vi.fn(async (url: string, init?: { method?: string }) => {
      if (url === Api) return response(true, [{ entryId: Entry.id, paeDocUrl: 'https://evil.example/pae.json' }]);
      return init?.method === 'HEAD' && url === Json ? response(true) : response(false);
    });
    expect(await resolveAfdbPaeUrl(Entry, fetcher)).toBe(Json);
  });

  it('says so when the entry has none, or the service cannot be reached', async () => {
    await expect(resolveAfdbPaeUrl(Entry, async () => response(false))).rejects.toThrow(/has no PAE for AF-P04637-F1/);
    await expect(
      resolveAfdbPaeUrl(Entry, async () => {
        throw new TypeError('fetch failed');
      }),
    ).rejects.toThrow(/Could not reach AlphaFold DB \(fetch failed\)/);
  });

  it('reads the JSON as a row-major matrix (row: aligned, column: scored), or says why not', () => {
    const m = paeValuesFromAfdb([
      {
        predicted_aligned_error: [
          [0, 1],
          [2, 0],
        ],
        max_predicted_aligned_error: 31.75,
      },
    ]);
    expect(m).toMatchObject({ size: 2, max: 2 });
    if (typeof m === 'string') throw new Error(m);
    expect([...m.values]).toEqual([0, 1, 2, 0]);
    expect(typeof paeValuesFromAfdb({ nope: 1 })).toBe('string');
  });
});

async function until(check: () => boolean, ms = 5000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

const fixture = () =>
  new File(
    [
      readFileSync(
        new URL(
          '../../../../packages/prism/viewer/src/sequence/molstar/_test/fixtures/structures/AF-P04637-F1-model_v6.bcif',
          import.meta.url,
        ),
      ),
    ],
    'AF-P04637-F1.bcif',
  );

const viewers: PrismViewer[] = [];
afterEach(async () => {
  vi.unstubAllGlobals();
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer() {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ModelArchiveExtension] },
  });
  viewers.push(v);
  const r = await v.api.load.files({ files: [fixture()] });
  if (!r.ok) throw new Error(r.message);
  return v;
}

const structureRef = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0].cell.transform.ref;
const residuesOf = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data.model.atomicHierarchy.residues._rowCount;

/** A stub of the two requests, the API and the JSON: a PAE of `n` residues whose cell (i, j) is `(7i + j) % 30`. */
function stubAfdb(n: number, status = 200) {
  const rows = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (7 * i + j) % 30));
  const seen: string[] = [];
  vi.stubGlobal(
    'fetch',
    vi.fn(async (input: unknown, init?: { method?: string }) => {
      const url = String(input);
      seen.push(`${init?.method ?? 'GET'} ${url}`);
      if (url === Api) {
        return new Response(JSON.stringify([{ entryId: Entry.id, paeDocUrl: Json }]), { status: 200 });
      }
      if (url === Json)
        return new Response(status === 200 ? JSON.stringify([{ predicted_aligned_error: rows }]) : 'not found', {
          status,
        });
      return new Response('', { status: 404 });
    }),
  );
  return { rows, seen };
}

describe('fetching the PAE of an AlphaFold DB model', () => {
  it('attaches it as state cells; the area reads it as the embedded PAE is read, and undo removes it', async () => {
    const v = await viewer();
    const area = v.api['model-archive']!;
    const ref = structureRef(v);
    expect(area.afdbEntry(ref)).toMatchObject({ id: Entry.id });
    expect(area.hasPae(ref)).toBe(false);
    expect(area.canFetchPae(ref)).toBe(true);
    const before = await area.pae({});
    expect(before.ok && before.value?.hasPae).toBe(false);

    const n = residuesOf(v);
    const { rows, seen } = stubAfdb(n);
    const r = await area.fetchPae({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(seen).toEqual([`GET ${Api}`, `GET ${Json}`]);

    expect(area.hasPae(ref)).toBe(true);
    expect(area.canFetchPae(ref)).toBe(false);
    expect(area.paeSize(ref)).toEqual({ residues: n, chains: 1 });
    const m = area.paeMatrix(ref)!;
    expect(m.size).toBe(n);
    expect(m.values[3 * n + 5]).toBe(rows[3][5]);
    expect(m.values[5 * n + 3]).toBe(rows[5][3]);
    expect(m.residues[0]).toMatchObject({ chain: 'A', seq: 1 });
    expect(m.plddt?.length).toBe(n);
    const info = await area.pae({});
    expect(info.ok && info.value).toMatchObject({ hasPae: true, size: n });
    // Fetching again is already so.
    const again = await area.fetchPae({});
    expect(again.ok && again.unchanged).toBe(true);

    await v.managers.history.undo();
    expect(area.hasPae(ref)).toBe(false);
  });

  it("offers the fetch in the structure row's menu and as a command, until the PAE is there", async () => {
    const v = await viewer();
    await until(() => v.views.scene.roots.value.length > 0);
    const entry = () => v.views.scene.roots.value[0].menu?.find((m) => m.id === 'model-archive.fetch-pae');
    expect(entry()?.label).toBe('Fetch PAE from AlphaFold DB');
    expect(v.managers.commands.get('model-archive.fetch-pae')).toBeDefined();
    stubAfdb(residuesOf(v));
    entry()!.run!();
    await until(() => v.api['model-archive']!.hasPae(structureRef(v)));
    await until(() => !entry());
  });

  it('fails with a message, and leaves nothing in the tree, when the file is missing or does not fit', async () => {
    const v = await viewer();
    const area = v.api['model-archive']!;
    const ref = structureRef(v);
    const cells = v.plugin.state.data.cells.size;

    stubAfdb(residuesOf(v), 404);
    const missing = await area.fetchPae({});
    expect(!missing.ok && missing.message).toMatch(/Could not get the PAE of AF-P04637-F1/);
    expect(v.plugin.state.data.cells.size).toBe(cells);

    stubAfdb(residuesOf(v) + 1);
    const wrong = await area.fetchPae({});
    expect(!wrong.ok && wrong.message).toMatch(/does not fit/);
    expect(v.plugin.state.data.cells.size).toBe(cells);
    expect(area.hasPae(ref)).toBe(false);
  });

  it('fails without a request for a structure that is no AlphaFold DB entry', async () => {
    const base = DefaultPrismViewerSpec();
    const v = await createPrismViewerHeadless({
      spec: { ...base, extensions: [...(base.extensions ?? []), ModelArchiveExtension] },
    });
    viewers.push(v);
    const cif = new File(
      [readFileSync(new URL('../../../../data/examples/1cbs_updated.cif', import.meta.url))],
      '1cbs.cif',
    );
    await v.api.load.files({ files: [cif] });
    const fetchStub = vi.fn();
    vi.stubGlobal('fetch', fetchStub);
    const r = await v.api['model-archive']!.fetchPae({});
    expect(!r.ok && r.code).toBe('blocked');
    expect(fetchStub).not.toHaveBeenCalled();
  });
});
