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
import { StructureElement } from '@molstar/model/model/structure';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { firstValueFrom } from 'rxjs';
import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  type PaeMatrix,
  paeBoxMean,
  paeCellAt,
  paeCellLabel,
  paeChains,
  paeImage,
  paePairTarget,
  paeRangeLabel,
  paeRangeTarget,
  paeTicks,
} from '../pae/data.js';
import { paeUI } from '../pae/ui.js';
import { PlddtBands, PlddtNoScore, plddtColor } from '../quality-assessment/plddt-bands.js';
import { QualityAssessment } from '../quality-assessment/prop.js';
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
/** Their pLDDT: one in each band. */
const Plddt = [95, 80, 60, 30];

/** `partialTo`: the metric covers only the residues of that chain (a domain's PAE of a bigger model). */
function modelCif(partialTo?: string): string {
  const atoms = Residues.map(
    ([chain, seq, comp], k) =>
      `ATOM ${k + 1} C CA . ${comp} ${chain} 1 ${seq} ? ${(3.8 * k).toFixed(3)} 0.000 0.000 1 90.0 ${seq} ${chain} 1`,
  );
  const pairs: string[] = [];
  let n = 0;
  Residues.forEach(([c1, s1, r1], i) => {
    Residues.forEach(([c2, s2, r2], j) => {
      if (partialTo && (c1 !== partialTo || c2 !== partialTo)) return;
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
    ...['ordinal_id', 'model_id', 'label_asym_id', 'label_seq_id', 'label_comp_id', 'metric_id', 'metric_value'].map(
      (f) => `_ma_qa_metric_local.${f}`,
    ),
    ...Residues.map(([chain, seq, comp], k) => `${k + 1} 1 ${chain} ${seq} ${comp} 1 ${Plddt[k].toFixed(2)}`),
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

const modelFile = (partialTo?: string) => new File([modelCif(partialTo)], 'pae-test.cif');
// jsdom's `import.meta.url` is not a file URL: the path from the folder.
const cbs = () =>
  new File([readFileSync(resolve(import.meta.dirname, '../../../../../data/examples/1cbs_updated.cif'))], '1cbs.cif');

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
    // x is the scored residue (j), y the aligned (i), as AlphaFold DB draws it.
    expect(paeCellAt(tiny, 0.75, 0.25)).toEqual({ i: 0, j: 1 });
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
    // The residues' pLDDT, for the plot's strip; a band each.
    expect([...m.plddt!]).toEqual(Plddt);
    expect([...m.plddt!].map(plddtColor)).toEqual(PlddtBands.map((b) => b.color));
    expect(plddtColor(Number.NaN)).toBe(PlddtNoScore);
    expect(v.api['model-archive']!.hasPae(refOf(v))).toBe(true);
  });

  it("says a PAE's size without reading its matrix, then the matrix's once read", async () => {
    const v = await viewer([modelFile()]);
    const area = v.api['model-archive']!;
    const read = vi.spyOn(QualityAssessment, 'pairwiseMetricFromModelArchiveCIF');
    // Before the card: the model's polymer residues and chains, no pairwise row read.
    expect(area.paeSize(refOf(v))).toEqual({ residues: 4, chains: 2 });
    expect(read).not.toHaveBeenCalled();
    area.paeMatrix(refOf(v));
    expect(read).toHaveBeenCalledTimes(1);
    expect(area.paeSize(refOf(v))).toEqual({ residues: 4, chains: 2 });
    const cbsViewer = await viewer([cbs()]);
    expect(cbsViewer.api['model-archive']!.paeSize(refOf(cbsViewer))).toBeUndefined();
  });

  it("sizes a partial metric by the residues it covers, not the model's polymer", async () => {
    const v = await viewer([modelFile('A')]);
    const area = v.api['model-archive']!;
    expect(area.paeSize(refOf(v))).toEqual({ residues: 3, chains: 1 });
    expect(area.paeMatrix(refOf(v))?.size).toBe(3);
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
  it('draws the matrix; the pointer lights its cell, a press selects it, a box marks its two ranges, moves and resizes until cleared', async () => {
    const put = vi.fn();
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(
      () =>
        ({
          createImageData: (w: number, h: number) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
          putImageData: put,
        }) as never,
    );
    vi.spyOn(HTMLElement.prototype, 'getBoundingClientRect').mockReturnValue({
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
    const unmount = ui.overlays![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    const canvas = host.querySelector('canvas')!;
    expect(canvas.width).toBe(4);
    // The matrix, and the two pLDDT strips.
    expect(put).toHaveBeenCalledTimes(3);
    const plot = host.querySelector<HTMLElement>('[role="img"]')!;

    // The pLDDT strips: along both axes with the bands' legend, until the toggle hides them.
    const toggle = [...host.querySelectorAll('button')].find((b) => b.textContent === 'pLDDT')!;
    const legend = host.querySelector<HTMLElement>('[aria-label="pLDDT bands"]')!;
    expect(toggle.getAttribute('aria-pressed')).toBe('true');
    expect(legend.style.display).toBe('flex');
    expect(legend.textContent).toContain('Very high > 90');
    toggle.click();
    expect(toggle.getAttribute('aria-pressed')).toBe('false');
    expect(legend.style.display).toBe('none');
    toggle.click();
    expect(legend.style.display).toBe('flex');
    const pointer = (type: string, x: number, y: number) =>
      // The main button held but on the release (a drag ends on a move without it, as a lost release).
      plot.dispatchEvent(
        new MouseEvent(type, {
          clientX: x,
          clientY: y,
          button: 0,
          buttons: type === 'pointerup' ? 0 : 1,
          bubbles: true,
        }),
      );

    // (300, 50): aligned on residue 0 (A 1), residue 3 (B 1) scored; lit in the sequence too (the marks' hover).
    pointer('pointermove', 300, 50);
    await until(() => host.textContent?.includes('A 1 ALA · B 1 LYS · 7.0 Å') ?? false);
    await until(() => v.sources.marks.hover.value !== undefined);
    expect(v.sources.marks.hover.value?.owner).toBe('highlight');

    // A press without a drag picks the cell: its two residues selected.
    pointer('pointerdown', 300, 50);
    pointer('pointerup', 300, 50);
    await until(() => v.sources.interaction.selection.value !== undefined);
    const selected = await v.api.select.get({});
    expect(selected.ok && JSON.stringify(selected.value)).toMatch(/"residues":2|2 residues/);
    // The same pair again: unselected; and once more, selected.
    pointer('pointerdown', 300, 50);
    pointer('pointerup', 300, 50);
    await until(() => v.sources.interaction.selection.value === undefined);
    pointer('pointerdown', 300, 50);
    pointer('pointerup', 300, 50);
    await until(() => v.sources.interaction.selection.value !== undefined);

    // A drag from (0, 0) to (3, 2): aligned A 1–3 and B 1, scored A 1–3; marked, summed up, its mean error.
    const steps = v.managers.history.steps.value.length;
    pointer('pointerdown', 50, 50);
    pointer('pointermove', 250, 350);
    pointer('pointerup', 250, 350);
    await until(() => v.sources.marks.marks.value.length === 2);
    const marks = v.sources.marks.marks.value;
    expect(marks.map((m) => [m.owner, m.id, m.tone, m.label])).toEqual([
      ['model-archive.pae', 'aligned', 'highlight', 'Aligned A 1–3, B 1'],
      ['model-archive.pae', 'scored', 'accent', 'Scored A 1–3'],
    ]);
    expect(marks.map((m) => StructureElement.Loci.size(m.loci))).toEqual([4, 3]);
    // The box took the pair the plot selected off again: a pair or a box, one at a time.
    await until(() => v.sources.interaction.selection.value === undefined);
    // And while the box is shown, a click selects nothing.
    pointer('pointermove', 300, 50);
    pointer('pointerdown', 300, 50);
    pointer('pointerup', 300, 50);
    await new Promise((r) => setTimeout(r, 50));
    expect(v.sources.interaction.selection.value).toBeUndefined();
    expect(host.textContent).toContain('Aligned A 1–3, B 1');
    expect(host.textContent).toContain('Mean error 3.3 Å');
    // Marks are view state: no step.
    expect(v.managers.history.steps.value.length).toBe(steps);

    // Dragged from inside, the box moves (one column right); by its left edge, it resizes (back to column 0).
    const scored = () => v.sources.marks.marks.value.find((m) => m.id === 'scored')?.label;
    pointer('pointermove', 150, 200);
    await new Promise((r) => requestAnimationFrame(r));
    expect(plot.style.cursor).toBe('move');
    pointer('pointerdown', 150, 200);
    pointer('pointermove', 250, 200);
    pointer('pointerup', 250, 200);
    await until(() => scored() === 'Scored A 2–3, B 1');
    pointer('pointerdown', 102, 200);
    pointer('pointermove', 10, 200);
    pointer('pointerup', 10, 200);
    await until(() => scored() === 'Scored A 1–3, B 1');

    // The ✕ on the box's corner clears the box and its marks.
    const corner = host.querySelector<HTMLButtonElement>('button[aria-label="Clear the box"]')!;
    expect(corner.style.display).toBe('grid');
    corner.click();
    await until(() => v.sources.marks.marks.value.length === 0);
    expect(host.textContent).not.toContain('Mean error');

    // A drag whose capture Chrome drops just before its release, then sends the release elsewhere in the page: the
    // release ends it, and the box drawn stays (after the ✕ there is none to go back to).
    pointer('pointerdown', 50, 50);
    pointer('pointermove', 250, 350);
    plot.dispatchEvent(new MouseEvent('lostpointercapture', { clientX: 250, clientY: 350, bubbles: true }));
    document.body.dispatchEvent(new MouseEvent('pointerup', { clientX: 250, clientY: 350, button: 0, bubbles: true }));
    await until(() => v.sources.marks.marks.value.length === 2);
    expect(host.textContent).toContain('Mean error 3.3 Å');
    expect(corner.style.display).toBe('grid');

    // A pair picked, then a selection of the person's own: a box leaves that selection alone.
    corner.click();
    await until(() => v.sources.marks.marks.value.length === 0);
    pointer('pointerdown', 300, 50);
    pointer('pointerup', 300, 50);
    await until(() => v.sources.interaction.selection.value !== undefined);
    expect((await v.api.select.all({})).ok).toBe(true);
    pointer('pointerdown', 50, 50);
    pointer('pointermove', 250, 350);
    pointer('pointerup', 250, 350);
    await until(() => v.sources.marks.marks.value.length === 2);
    const kept = await v.api.select.get({});
    expect(kept.ok && JSON.stringify(kept.value)).toMatch(/"residues":4/);
    unmount();
  });

  it('brings the crosshair back after a box drag ends on the cell it started on', async () => {
    vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(
      () =>
        ({
          createImageData: (w: number, h: number) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
          putImageData: vi.fn(),
        }) as never,
    );
    vi.spyOn(HTMLElement.prototype, 'getBoundingClientRect').mockReturnValue({
      left: 0,
      top: 0,
      width: 400,
      height: 400,
    } as DOMRect);
    const v = await viewer([modelFile()]);
    const ui = paeUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    const host = document.createElement('div');
    document.body.append(host);
    const unmount = ui.overlays![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    const plot = host.querySelector<HTMLElement>('[role="img"]')!;
    const pointer = (type: string, x: number, y: number) =>
      // The main button held but on the release (a drag ends on a move without it, as a lost release).
      plot.dispatchEvent(
        new MouseEvent(type, {
          clientX: x,
          clientY: y,
          button: 0,
          buttons: type === 'pointerup' ? 0 : 1,
          bubbles: true,
        }),
      );
    const frame = () => new Promise((r) => requestAnimationFrame(r));
    const hair = () => plot.querySelectorAll('line[stroke="rgba(255, 255, 255, 0.7)"][visibility="visible"]').length;

    pointer('pointermove', 300, 50);
    await frame();
    expect(hair()).toBe(2);
    // Down on the cell, away (the crosshair goes), back and up on the same cell: it comes back.
    pointer('pointerdown', 300, 50);
    pointer('pointermove', 100, 250);
    await frame();
    expect(hair()).toBe(0);
    pointer('pointermove', 300, 50);
    pointer('pointerup', 300, 50);
    await frame();
    expect(hair()).toBe(2);
    unmount();
  });

  it('has no toolbar button while no loaded model has a PAE', async () => {
    const v = await viewer([cbs()]);
    const ui = paeUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    expect(await firstValueFrom(ui.toolbar![0].visible!)).toBe(false);
  });

  it('has a toolbar button for an AlphaFold DB entry that has no PAE yet, to fetch it from the card', async () => {
    const af = new File(
      [
        readFileSync(
          resolve(
            import.meta.dirname,
            '../../../../../packages/prism/viewer/src/sequence/molstar/_test/fixtures/structures/AF-P04637-F1-model_v6.bcif',
          ),
        ),
      ],
      'AF-P04637-F1.bcif',
    );
    const v = await viewer([af]);
    const ui = paeUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    expect(await firstValueFrom(ui.toolbar![0].visible!)).toBe(true);
  });

  it('ticks residue numbers at a round step per chain, names runs as targets and words, and averages a box', () => {
    // Chain A 1–120, chain B 1–80.
    const residues = [
      ...Array.from({ length: 120 }, (_, k) => ({ chain: 'A', seq: k + 1, insCode: '', name: 'ALA' })),
      ...Array.from({ length: 80 }, (_, k) => ({ chain: 'B', seq: k + 1, insCode: '', name: 'GLY' })),
    ];
    const size = residues.length;
    const values = new Float32Array(size * size).map((_, k) => (Math.floor(k / size) < 120 ? 2 : 20));
    const big: PaeMatrix = { name: 'PAE', size, values, max: 20, residues };
    expect(paeTicks(big, 6).map((t) => `${residues[t.at].chain}${t.label}`)).toEqual(['A50', 'A100', 'B50']);
    expect(paeRangeTarget(big, 54, 99)).toBe('A:55-100');
    expect(paeRangeTarget(big, 110, 124)).toBe('A:111-120, B:1-5');
    expect(paeRangeLabel(big, 124, 110)).toBe('A 111–120, B 1–5');
    // Rows 0–1 are 2, rows 120 on are 20.
    expect(paeBoxMean(big, 0, 1, 0, 9)).toBe(2);
    expect(paeBoxMean(big, 119, 120, 0, 0)).toBe(11);
  });
});
