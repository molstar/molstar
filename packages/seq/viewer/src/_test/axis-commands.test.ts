/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The axis in the viewer (_spec/README.md "Commands", "ViewModel" History, "Bindings"): `section.reverse` and
// `section.collapse` as patches with their keys in `view.keys`; back and forward across them landing on the same
// positions; a pan whose window ends on a break showing nothing hidden; hits on a reversed collapsed axis with residue
// rows through each reader (the cursor, a drag, a click's selection, the focus); a break's hit target.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { hitView } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { hitTarget } from '../models/hover.js';
import { positionAt } from '../models/axis.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Protein = '/views/main/protein';
const Gene = '/views/main/gene';

// p53-like and small: 4 residues on a 30-nucleotide gene, exons [2, 7) and [12, 19), minus strand
function spec(gene: Record<string, unknown> = {}) {
  const b = createSeqSpec({ title: 'axis' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGTACGGATCCGTAGCATGCATCGATCGA' });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [2, 12], end: [7, 19] }, strand: '-' });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p });
    protein.track('letters').layer('l', { type: 'letters' });
    // `lead: true`: the gene leads (no follow)
    const { lead, ...rest } = gene;
    const s = v.section('gene', { of: g, ...(lead ? {} : { follow: protein }), ...(rest as object) });
    s.track('residues').layer('r', { type: 'letters', data: p });
    s.track('bases').layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make(gene: Record<string, unknown> = {}) {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    history: uiHistory(),
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewers.push(v);
  await v.load(spec(gene));
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

const param = (v: SeqViewer, address: string, key: string) => {
  let n = v.spec.root.value;
  for (const id of address.split('/').filter(Boolean)) n = n.children!.find((c) => c.id === id)!;
  return n.params?.[key];
};
const axisOf = (v: SeqViewer, address: string) =>
  v.view.resolved.value!.sections.find((s) => s.address === address)!.axis;

describe('section.reverse and section.collapse', () => {
  it('each is one patch and one undo step; follow becomes the explicit opposite of what it resolved to', async () => {
    const v = await make();
    expect(axisOf(v, Gene).dir).toBe(-1);
    await v.execute('section.reverse', { address: Gene });
    expect(param(v, Gene, 'direction')).toBe('genomic');
    expect(axisOf(v, Gene).dir).toBe(1);
    await v.execute('section.reverse', { address: Gene });
    expect(param(v, Gene, 'direction')).toBe('reversed');
    await v.execute('section.collapse', { address: Gene });
    expect(param(v, Gene, 'pieces')).toBe('coding');
    expect(axisOf(v, Gene).ranges).toEqual([12, 19, 2, 7]);
    await v.execute('section.collapse', { address: Gene });
    expect(param(v, Gene, 'pieces')).toBeUndefined();
    // undo returns each change, one step each
    expect((await (v.spec.history as UiHistory).stack.undo()).ok).toBe(true);
    expect(param(v, Gene, 'pieces')).toBe('coding');
  });

  it('acts on the keyboard cursor’s section, else the first leading one; its keys are in `view.keys`', async () => {
    const v = await make();
    expect(v.config.view.keys['section.reverse']).toEqual(['Shift+R']);
    expect(v.config.view.keys['section.collapse']).toEqual(['Shift+C']);
    await v.execute('section.reverse', {});
    expect(param(v, Protein, 'direction')).toBe('reversed');
    v.selection.setCursor({ section: Gene, system: '/data/g', position: 5 });
    await v.execute('section.collapse', {});
    expect(param(v, Gene, 'pieces')).toBe('coding');
  });
});

describe('the view’s history over a changed axis', () => {
  const windowOf = (v: SeqViewer, address: string) => {
    const g = v.view.get(address)!;
    const a = positionAt(axisOf(v, address), g.start, 'start')!;
    const b = positionAt(axisOf(v, address), g.end, 'end')!;
    return [Math.min(a, b), Math.max(a, b)];
  };

  it('a reverse, then back: the exact window it left; a collapse, then forward: the exact window', async () => {
    const v = await make({ lead: true });
    v.view.jump(() => v.view.setRange(Gene, 13, 18));
    v.view.jump(() => v.view.setRange(Gene, 3, 6));
    await v.execute('section.reverse', { address: Gene });
    expect(axisOf(v, Gene).dir).toBe(-1);
    expect(v.view.back()).toBe(true);
    expect(windowOf(v, Gene)).toEqual([13, 18]);
    await v.execute('section.collapse', { address: Gene });
    expect(axisOf(v, Gene).pieces).toBe('coding');
    expect(v.view.forward()).toBe(true);
    expect(windowOf(v, Gene)).toEqual([3, 6]);
    expect(v.view.back()).toBe(true);
    expect(windowOf(v, Gene)).toEqual([13, 18]);
  });

  it('back and forward across a reverse patch and a collapse patch land on the same positions', async () => {
    const v = await make();
    // a jump on the gene's leader: its window on the gene
    v.view.jump(() => v.view.setRange(Protein, 1, 3));
    const window = () => {
      const g = v.view.get(Gene)!;
      const a = positionAt(axisOf(v, Gene), g.start, 'start')!;
      const b = positionAt(axisOf(v, Gene), g.end, 'end')!;
      return [Math.min(a, b), Math.max(a, b)];
    };
    v.view.jump(() => v.view.setRange(Protein, 0, 2));
    const before = window();
    await v.execute('section.collapse', { address: Gene });
    await v.execute('section.reverse', { address: Protein });
    // back: to the protein's [1, 3), on the axes as they are now
    expect(v.view.back()).toBe(true);
    const p = v.view.get(Protein)!;
    const pa = positionAt(axisOf(v, Protein), p.start, 'start')!;
    const pb = positionAt(axisOf(v, Protein), p.end, 'end')!;
    expect([Math.min(pa, pb), Math.max(pa, pb)]).toEqual([1, 3]);
    expect(v.view.forward()).toBe(true);
    // the gene shows what it showed, its positions the same though its axis turned and collapsed
    const after = window();
    expect(after[0]).toBeGreaterThanOrEqual(before[0] - 1e-9);
    expect(after[1]).toBeLessThanOrEqual(before[1] + 1e-9);
  });

  it('a pan whose window ends on a break shows nothing hidden', async () => {
    const v = await make({ pieces: '/data/m', direction: 'genomic', lead: true });
    // the gene leads: its coding axis [2, 7) then [12, 19); a window ending at offset 5 ends on the break
    v.view.setRange(Gene, 1, 5);
    const g = v.view.get(Gene)!;
    expect(positionAt(axisOf(v, Gene), g.end, 'end')).toBe(7);
    expect(positionAt(axisOf(v, Gene), g.start, 'start')).toBe(3);
  });
});

describe('hits on a reversed collapsed axis with a residue row', () => {
  it('through each reader: the hit’s position, a click’s selection of the residue, the cursor, the focus', async () => {
    const v = await make({ pieces: 'coding' });
    v.draw();
    const layout = v.frame!.layout;
    const sec = layout.sections[1];
    const axis = sec.axis!;
    // the residue row: residue 2 (nucleotides 12, 6, 5) is split by the intron
    const box = sec.tracks[0].box;
    axis.cell(12);
    const h = hitView(layout, (axis.sx0 + axis.sx1) / 2, box.y + 2)!;
    expect(h.position).toBe(12);
    expect(h.extent?.parts.map((p) => [p.system, [...p.ranges]])).toEqual([['/data/p', [2, 3]]]);
    expect(hitTarget(h)).toBe('cell');
    // a click selects the residue on the protein; its codon is what the gene shows of it
    v.selection.select(h, axisOf(v, Gene));
    expect(v.selection.value.value.location.parts.map((p) => [p.system, [...p.ranges]])).toEqual([['/data/p', [2, 3]]]);
    // the cursor walks the gene in display order across the break, and says so
    v.selection.setCursor({ section: Gene, system: '/data/g', position: 12, keyboard: true });
    const r = await v.execute('cursor.next', { section: Gene });
    expect(v.selection.cursor.value?.position).toBe(6);
    expect(r).toContain('crossed a break, 5 positions hidden');
  });

  it('a drag through the input on the reversed collapsed bases selects what it crossed, in display order', async () => {
    const v = await make({ pieces: 'coding' });
    v.draw();
    const sec = v.frame!.layout.sections[1];
    const axis = sec.axis!;
    const box = sec.tracks[1].box;
    const at = (p: number) => {
      axis.cell(p);
      return (axis.sx0 + axis.sx1) / 2;
    };
    // from 13 to 5 on screen: 13, 12, then across the break 6, 5; each end snapped to its whole codon (EC-149), so the
    // codons (13, 14, 15) and (12, 6, 5) are selected whole
    expect(v.input.down(at(13), box.y + 2)).toBe(true);
    v.input.move(at(5), box.y + 2);
    v.input.up();
    const parts = v.selection.value.value.location.parts.map((p) => [p.system, [...p.ranges]]);
    expect(parts).toEqual([['/data/g', [5, 7, 12, 16]]]);
  });

  it('a break is its own hit target near its notch; a focus into a hidden stretch moves to its break and says so', async () => {
    const v = await make({ pieces: '/data/m', direction: 'reversed', lead: true });
    v.draw();
    const sec = v.frame!.layout.sections[1];
    const [x] = sec.breaks;
    expect(sec.breaks.length).toBe(3);
    const h = hitView(v.frame!.layout, x + 1, sec.band!.y - 2)!;
    expect(hitTarget(h)).toBe('break');
    expect(h.extent?.parts.map((p) => [...p.ranges])).toEqual([[7, 12]]);
    // the gene leads here: a focus on an intronic position
    const out = v.addresses.focus('/data/g/(range :from 8 :to 10)', { section: Gene });
    expect(out.message).toContain('hidden by collapsed introns');
  });

  it('a focus outside the section’s range is not said to be hidden by introns', async () => {
    const v = await make({ range: [10, 25], lead: true });
    const out = v.addresses.focus('/data/g/(range :from 2 :to 4)', { section: Gene });
    expect(out.message).not.toContain('hidden by collapsed introns');
    expect(v.addresses.movedToBreak).toBe(false);
  });
});

describe('a web font arriving', () => {
  it('measures the label column again and gives the sections their data width with it', async () => {
    const v = await make();
    const set = vi.spyOn(v.view, 'setWidth');
    v.invalidateText();
    expect(set).toHaveBeenCalled();
  });
});
