/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { zip } from '@molstar/core/util/zip/zip';
import { CIF } from '@molstar/io/reader/cif';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { associatedPaeOf } from '../pae/associated.js';
import { ModelArchiveExtension } from '../prism-viewer.js';

// Model Archive's PAE in an associated file (`pae/associated.ts`, parity.md §6 "PAE plot"): the model names it, in a zip
// with the MSA (as ma-bak-cepc-0003) or as a CIF of its own. A tiny ModelCIF and its PAE file are written here (chain A,
// 3 residues; the PAE is 2·|i − j| + 1); every request is a stub: no network.

const Residues = [
  [1, 'ALA'],
  [2, 'GLY'],
  [3, 'SER'],
] as const;
const pae = (i: number, j: number) => 2 * Math.abs(i - j) + 1;
const Zip = 'https://modelarchive.org/api/projects/ma-test-0001?type=materials_procedures__accompanying_data_file_name';
const Path = 'ma-test-0001_predicted_aligned_error_v1.cif';

const metrics = [
  'loop_',
  '_ma_qa_metric.id',
  '_ma_qa_metric.name',
  '_ma_qa_metric.type',
  '_ma_qa_metric.mode',
  '_ma_qa_metric.software_group_id',
  '1 PAE PAE local-pairwise .',
  '#',
];

/** The model, naming its PAE file: `zip` (an entry of an archive), `cif` (a file of its own, at `url`), or none. */
function modelCif(associated?: { kind: 'zip' } | { kind: 'cif'; url: string }): string {
  const atoms = Residues.map(
    ([seq, comp], k) =>
      `ATOM ${k + 1} C CA . ${comp} A 1 ${seq} ? ${(3.8 * k).toFixed(3)} 0.000 0.000 1 90.0 ${seq} A 1`,
  );
  const files =
    associated?.kind === 'zip'
      ? [
          'loop_',
          ...['id', 'entry_id', 'file_url', 'file_type', 'file_format', 'file_content', 'details'].map(
            (f) => `_ma_entry_associated_files.${f}`,
          ),
          `1 ma-test-0001 '${Zip}' archive zip 'archive with multiple files' .`,
          '#',
          'loop_',
          ...['id', 'archive_file_id', 'file_path', 'file_format', 'file_content', 'description'].map(
            (f) => `_ma_associated_archive_file_details.${f}`,
          ),
          `1 1 ${Path} cif 'local pairwise QA scores' ?`,
          "2 1 A.a3m a3m 'multiple sequence alignments' ?",
          '#',
        ]
      : associated?.kind === 'cif'
        ? [
            'loop_',
            ...['id', 'entry_id', 'file_url', 'file_type', 'file_format', 'file_content'].map(
              (f) => `_ma_entry_associated_files.${f}`,
            ),
            `1 ma-test-0001 '${associated.url}' file cif 'local pairwise QA scores'`,
            '#',
          ]
        : [];
  return [
    'data_ma-test-0001',
    '_entry.id ma-test-0001',
    '#',
    'loop_',
    '_entity.id',
    '_entity.type',
    '1 polymer',
    '#',
    'loop_',
    ...[
      'group_PDB',
      'id',
      'type_symbol',
      'label_atom_id',
      'label_alt_id',
      'label_comp_id',
      'label_asym_id',
      'label_entity_id',
      'label_seq_id',
      'pdbx_PDB_ins_code',
      'Cartn_x',
      'Cartn_y',
      'Cartn_z',
      'occupancy',
      'B_iso_or_equiv',
      'auth_seq_id',
      'auth_asym_id',
      'pdbx_PDB_model_num',
    ].map((f) => `_atom_site.${f}`),
    ...atoms,
    '#',
    // The model lists the metric; its values are in the associated file.
    ...metrics,
    ...files,
    '',
  ].join('\n');
}

/** The associated file: the pairwise metric over the model's residues (`chain`: another model's chain, fitting none). */
function paeCif(chain = 'A'): string {
  const pairs: string[] = [];
  let n = 0;
  Residues.forEach(([s1, r1], i) => {
    Residues.forEach(([s2, r2], j) => {
      pairs.push(`${++n} 1 ${chain} ${s1} ${r1} ${chain} ${s2} ${r2} 1 ${pae(i, j).toFixed(2)}`);
    });
  });
  return [
    'data_ma-test-0001_predicted_aligned_error_v1',
    ...metrics,
    'loop_',
    ...[
      'ordinal_id',
      'model_id',
      'label_asym_id_1',
      'label_seq_id_1',
      'label_comp_id_1',
      'label_asym_id_2',
      'label_seq_id_2',
      'label_comp_id_2',
      'metric_id',
      'metric_value',
    ].map((f) => `_ma_qa_metric_local_pairwise.${f}`),
    ...pairs,
    '#',
    '',
  ].join('\n');
}

const bytes = (s: string) => new TextEncoder().encode(s);
const archive = async (entries: Record<string, string>) =>
  zip(SyncRuntimeContext, Object.fromEntries(Object.entries(entries).map(([k, v]) => [k, bytes(v)])));

async function frameOf(cif: string) {
  const parsed = await CIF.parse(cif).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result.blocks[0];
}

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
