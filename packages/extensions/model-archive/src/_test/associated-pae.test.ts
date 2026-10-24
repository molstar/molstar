/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { associatedPaeOf } from '../pae/associated.js';
import { ModelArchiveExtension } from '../prism-viewer.js';
import { archive, frameOf, modelCif, Path, pae, paeCif, Zip } from './associated-fixture.js';

// Model Archive's PAE in an associated file (`pae/associated.ts`, parity.md §6 "PAE plot") over the generated model of
// `associated-fixture.ts`; every request is a stub: no network.

describe('the associated file a model names for its PAE', () => {
  it('finds an entry of an archive and a file of its own, by its content', async () => {
    expect(associatedPaeOf(await frameOf(modelCif({ kind: 'zip' })))).toEqual({ url: Zip, path: Path });
    const url = 'https://modelarchive.org/files/pae.cif';
    expect(associatedPaeOf(await frameOf(modelCif({ kind: 'cif', url })))).toEqual({ url });
    expect(associatedPaeOf(await frameOf(modelCif()))).toBeUndefined();
  });

  it('resolves a relative address against where the model came from, and drops it without one', async () => {
    const frame = await frameOf(modelCif({ kind: 'cif', url: 'pae.cif' }));
    expect(associatedPaeOf(frame, 'https://example.org/models/m.cif')).toEqual({
      url: 'https://example.org/models/pae.cif',
    });
    expect(associatedPaeOf(frame)).toBeUndefined();
  });
});

const viewers: PrismViewer[] = [];
afterEach(async () => {
  vi.unstubAllGlobals();
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(model: string) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ModelArchiveExtension] },
  });
  viewers.push(v);
  const r = await v.api.load.files({ files: [new File([model], 'ma-test-0001.cif')] });
  if (!r.ok) throw new Error(r.message);
  return v;
}

/** Answers `url` with `body`, anything else with a 404; the requests made. */
function stub(url: string, body: Uint8Array | string) {
  const seen: string[] = [];
  vi.stubGlobal(
    'fetch',
    vi.fn(async (input: unknown) => {
      seen.push(String(input));
      return String(input) === url
        ? new Response(body as BodyInit, { status: 200 })
        : new Response('', { status: 404 });
    }),
  );
  return seen;
}

const structureRef = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures[0].cell.transform.ref;
const paeCells = (v: PrismViewer) =>
  [...v.plugin.state.data.cells.values()].filter((c) => c.transform.tags?.includes('model-archive.associated-pae'));

describe('fetching the PAE a Model Archive model names', () => {
  it('unpacks it from the zip into a cell under the model; the area reads it, and undo removes it', async () => {
    const v = await viewer(modelCif({ kind: 'zip' }));
    const area = v.api['model-archive']!;
    const ref = structureRef(v);
    // The model lists the metric but carries no value: no PAE yet, one to fetch.
    expect(area.hasPae(ref)).toBe(false);
    expect(area.canFetchPae(ref)).toBe(true);
    expect(area.afdbEntry(ref)).toBeUndefined();
    // The agent's policy sees the address before the run (Model Archive is a known service).
    expect(area.fetchPae.definition.addresses?.({ structure: ref })).toEqual([Zip]);
    // Named by its label (which the run resolves), or not at all: the address is listed all the same.
    expect(area.fetchPae.definition.addresses?.({ structure: 'ma-test-0001' })).toEqual([Zip]);
    expect(area.fetchPae.definition.addresses?.({ structure: '' })).toEqual([Zip]);

    const seen = stub(Zip, await archive({ [Path]: paeCif(), 'A.a3m': '>A\nAGS\n' }));
    const r = await area.fetchPae({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(r.ok && r.did).toMatch(/from Model Archive: 3 × 3 residues/);
    expect(seen).toEqual([Zip]);

    expect(area.hasPae(ref)).toBe(true);
    expect(area.canFetchPae(ref)).toBe(false);
    expect(area.paeSize(ref)).toEqual({ residues: 3, chains: 1 });
    const m = area.paeMatrix(ref)!;
    expect(m.size).toBe(3);
    expect(m.values[0 * 3 + 2]).toBe(pae(0, 2));
    expect(m.values[2 * 3 + 1]).toBe(pae(2, 1));
    expect(m.residues.map((x) => x.seq)).toEqual([1, 2, 3]);

    // A saved session keeps the address and the entry's path, under the model.
    const cells = paeCells(v);
    expect(cells).toHaveLength(1);
    expect(cells[0].transform.params).toMatchObject({ url: Zip, path: Path });

    await v.managers.history.undo();
    expect(area.hasPae(ref)).toBe(false);
    expect(paeCells(v)).toHaveLength(0);
  });

  it('reads a PAE file of its own, offered in the structure row menu', async () => {
    const url = 'https://modelarchive.org/files/ma-test-0001_pae.cif';
    const v = await viewer(modelCif({ kind: 'cif', url }));
    await until(() => v.views.scene.roots.value.length > 0);
    const entry = () => v.views.scene.roots.value[0].menu?.find((m) => m.id === 'model-archive.fetch-pae');
    expect(entry()?.label).toBe('Fetch PAE');
    stub(url, paeCif());
    entry()!.run!();
    await until(() => v.api['model-archive']!.hasPae(structureRef(v)));
    expect(v.api['model-archive']!.paeMatrix(structureRef(v))?.size).toBe(3);
    await until(() => !entry());
  });

  it('fails with a message, and leaves nothing in the tree, when the file is missing or does not fit', async () => {
    const v = await viewer(modelCif({ kind: 'zip' }));
    const area = v.api['model-archive']!;
    const cells = v.plugin.state.data.cells.size;

    stub(Zip, await archive({ 'other.cif': paeCif() }));
    const missing = await area.fetchPae({});
    expect(!missing.ok && missing.message).toMatch(/Could not get the PAE of ma-test-0001 .*has no/);
    expect(v.plugin.state.data.cells.size).toBe(cells);

    // An empty file fails as well, every time, and leaves nothing behind.
    stub(Zip, await archive({ [Path]: '# nothing\n' }));
    for (let k = 0; k < 2; k++) {
      const empty = await area.fetchPae({});
      expect(!empty.ok && empty.message).toMatch(/no data block/);
      expect(v.plugin.state.data.cells.size).toBe(cells);
    }

    stub(Zip, await archive({ [Path]: paeCif('Z') }));
    const wrong = await area.fetchPae({});
    expect(!wrong.ok && wrong.message).toMatch(/does not fit/);
    expect(v.plugin.state.data.cells.size).toBe(cells);
    expect(area.hasPae(structureRef(v))).toBe(false);
  });
});

async function until(ok: () => boolean, ms = 5000) {
  const end = Date.now() + ms;
  while (!ok()) {
    if (Date.now() > end) throw new Error('timed out');
    await new Promise((r) => setTimeout(r, 5));
  }
}
