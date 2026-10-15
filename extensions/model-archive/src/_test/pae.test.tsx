/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { firstValueFrom } from 'rxjs';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { type PaeMatrix, paeCellAt, paeCellLabel, paeChains, paeImage, paePairTarget } from '../pae/data.js';
import { paeUI } from '../pae/ui.js';
import { ModelArchiveExtension } from '../prism-viewer.js';

// The PAE (parity.md §6, "PAE plot"): no recorded fixture carries one (AlphaFold DB ships PAE apart from the model), so
// a tiny ModelCIF is written here: chains A (3 residues) and B (1), and a full `ma_qa_metric_local_pairwise` whose
// value is 2·|i − j| + 1. The drawing is checked from plain arrays; the card in jsdom with a stand-in 2D context.

const Residues = [
  ['A', 1, 'ALA'],
  ['A', 2, 'GLY'],
  ['A', 3, 'SER'],
  ['B', 1, 'LYS'],
] as const;
const pae = (i: number, j: number) => 2 * Math.abs(i - j) + 1;

function modelCif(): string {
  const atoms = Residues.map(
    ([chain, seq, comp], k) =>
      `ATOM ${k + 1} C CA . ${comp} ${chain} 1 ${seq} ? ${(3.8 * k).toFixed(3)} 0.000 0.000 1 90.0 ${seq} ${chain} 1`,
  );
  const pairs: string[] = [];
  let n = 0;
  Residues.forEach(([c1, s1, r1], i) => {
    Residues.forEach(([c2, s2, r2], j) => {
      pairs.push(`${++n} 1 ${c1} ${s1} ${r1} ${c2} ${s2} ${r2} 2 ${pae(i, j).toFixed(2)}`);
    });
  });
  return [
    'data_PAETEST',
    '_entry.id PAETEST',
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
    'loop_',
    '_ma_qa_metric.id',
    '_ma_qa_metric.name',
    '_ma_qa_metric.description',
    '_ma_qa_metric.type',
    '_ma_qa_metric.mode',
    '_ma_qa_metric.software_group_id',
    "1 pLDDT 'Predicted lDDT' pLDDT local .",
    "2 PAE 'Predicted aligned error' PAE local-pairwise .",
    '#',
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

const modelFile = () => new File([modelCif()], 'pae-test.cif');
// jsdom's `import.meta.url` is not a file URL: the path from the folder.
const cbs = () =>
  new File([readFileSync(resolve(import.meta.dirname, '../../../../data/examples/1cbs_updated.cif'))], '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

async function viewer(files: File[]) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ModelArchiveExtension] },
  });
  viewers.push(v);
  for (const f of files) {
    const r = await v.api.load.files({ files: [f] });
    if (!r.ok) throw new Error(r.message);
  }
  return v;
}

async function until(check: () => boolean, ms = 5000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

const refOf = (v: PrismViewer) => v.sources.structures.value[0].ref;

/** A 2×2 matrix by hand: one triangle given, the other mirrored in the picture. */
const tiny: PaeMatrix = {
  name: 'PAE',
  size: 2,
  values: new Float32Array([0, 10, Number.NaN, Number.NaN]),
  max: 10,
  residues: [
    { chain: 'A', seq: 5, insCode: '', name: 'ALA' },
    { chain: 'B', seq: 7, insCode: 'A', name: 'GLY' },
  ],
};

describe('the PAE drawing', () => {
  it('colors cells from dark (low) to light (high), mirrors a missing triangle, and leaves no value transparent', () => {
    const rgba = paeImage(tiny);
    const px = (x: number, y: number) => [...rgba.slice((y * 2 + x) * 4, (y * 2 + x) * 4 + 4)];
    expect(px(0, 0)).toEqual([0x00, 0x44, 0x1b, 255]);
    expect(px(0, 1)).toEqual([0xf7, 0xfc, 0xf5, 255]);
    // (1, 0) is the mirror of (0, 1); (1, 1) has none.
    expect(px(1, 0)).toEqual(px(0, 1));
    expect(px(1, 1)[3]).toBe(0);
    // The buffer is reused when it fits.
    expect(paeImage(tiny, rgba)).toBe(rgba);
  });

  it('names cells, their residues as a target, and the chains along the axis', () => {
    expect(paeCellAt(tiny, 0.75, 0.25)).toEqual({ i: 1, j: 0 });
    expect(paeCellAt(tiny, 1, 0.5)).toBeUndefined();
    expect(paePairTarget(tiny, 0, 1)).toBe('A:5, B:7A');
    expect(paePairTarget(tiny, 1, 1)).toBe('B:7A');
    expect(paeCellLabel(tiny, 1, 0)).toBe('B 7A GLY · A 5 ALA · 10.0 Å');
    expect(paeChains(tiny)).toEqual([
      { label: 'A', from: 0, to: 1 },
      { label: 'B', from: 1, to: 2 },
    ]);
  });
});

describe('model-archive.pae', () => {
  it("reads a ModelCIF's pairwise metric: whether it has one, its size and range", async () => {
    const v = await viewer([modelFile()]);
    expect(v.api.describe('model-archive.pae')?.definition).toMatchObject({ read: true });
    const steps = v.managers.history.steps.value.length;
    const r = await v.api['model-archive']!.pae({});
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(r.ok && r.value).toMatchObject({ hasPae: true, metric: 'PAE', size: 4, max: 7, chains: ['A', 'B'] });
    // A read: no step.
    expect(v.managers.history.steps.value.length).toBe(steps);

    const m = v.api['model-archive']!.paeMatrix(refOf(v))!;
    expect(m.residues.map((x) => `${x.chain}${x.seq}${x.name}`)).toEqual(['A1ALA', 'A2GLY', 'A3SER', 'B1LYS']);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) expect(m.values[i * 4 + j]).toBe(pae(i, j));
    expect(v.api['model-archive']!.hasPae(refOf(v))).toBe(true);
  });

  it('says so when a structure has none, and fails when nothing is loaded', async () => {
    const v = await viewer([cbs()]);
    const r = await v.api['model-archive']!.pae({});
    expect(r.ok && r.value).toMatchObject({ hasPae: false });
    expect(v.api['model-archive']!.hasPae(refOf(v))).toBe(false);
    const empty = await viewer([]);
    expect((await empty.api['model-archive']!.pae({})).ok).toBe(false);
  });

  it('brings the commands', async () => {
    const v = await viewer([]);
    expect(v.managers.commands.get('model-archive.pae')).toBeDefined();
    expect(v.managers.commands.get('model-archive.show-pae')?.title).toBe('Show the PAE plot');
  });
});

describe('the PAE card', () => {
  it('draws the matrix, names the cell under the pointer, and a click selects its two residues', async () => {
    const put = vi.fn();
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(
      () =>
        ({
          createImageData: (w: number, h: number) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
          putImageData: put,
        }) as never,
    );
    vi.spyOn(HTMLCanvasElement.prototype, 'getBoundingClientRect').mockReturnValue({
      left: 0,
      top: 0,
      width: 400,
      height: 400,
    } as DOMRect);
    const v = await viewer([modelFile()]);
    const ui = paeUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    expect(ui.overlays?.map((o) => o.id)).toEqual(['model-archive.overlay.pae']);
    expect(ui.toolbar?.[0]).toMatchObject({ id: 'model-archive.toolbar.pae', overlay: 'model-archive.overlay.pae' });
    expect(await firstValueFrom(ui.toolbar![0].visible!)).toBe(true);

    const host = document.createElement('div');
    document.body.append(host);
    ui.overlays![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    const canvas = host.querySelector('canvas')!;
    expect(canvas.width).toBe(4);
    expect(put).toHaveBeenCalledTimes(1);
    expect(host.textContent).toContain('B');

    // (300, 50): residue 3 (B 1) against residue 0 (A 1).
    canvas.dispatchEvent(new MouseEvent('pointermove', { clientX: 300, clientY: 50 }));
    await until(() => host.textContent?.includes('B 1 LYS · A 1 ALA · 7.0 Å') ?? false);
    canvas.dispatchEvent(new MouseEvent('click', { clientX: 300, clientY: 50 }));
    await until(() => v.sources.interaction.selection.value !== undefined);
    const selected = await v.api.select.get({});
    expect(selected.ok && JSON.stringify(selected.value)).toMatch(/"residues":2|2 residues/);
  });

  it('has no toolbar button while no loaded model has a PAE', async () => {
    const v = await viewer([cbs()]);
    const ui = paeUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    expect(await firstValueFrom(ui.toolbar![0].visible!)).toBe(false);
  });
});
