/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover, highlights, selection and the channels' marks (IN3; seqspec:interaction#hover, seqspec:interaction#channels,
// seqspec:interaction#selection, seqspec:interaction#focus): the hover's owners (ported from the managers test) and its
// external target, each channel's closure default, the selection's close and widen, the channels lit per section and
// rebuilt only when their own channel or the view changes (counters), rows sections' rectangles (a gap cell's too), the
// shown set's own systems and a host's Mol*-like external ones (a complex across chains), and held targets carried to
// every set (the host highlight across a lineage).

import { SelectionGranularity as G } from '@molstar/seq-core/format';
import {
  CoordinateSystemKind as K,
  type Loci,
  type MoleculeInput,
  OrderedSet,
  type SeqIndex,
  type SeqLoci,
  type SeqSet,
  SeqSets,
} from '@molstar/seq-core/model';
import { type SystemDescriptor, TranslateMaps } from '@molstar/seq-core/translate';
import { describe, expect, it } from 'vitest';
import { case7 } from '../../core/_test/format/cases.js';
import { gene, rows } from '../../core/_test/format/project-fixtures.js';
import { Highlight } from '../src/highlight.js';
import {
  ChannelMarks,
  createSeqViewer,
  type HighlightEntry,
  Hit,
  HitRegion,
  HoverBehavior,
  HoverKind,
  MarkBit,
  SelectBehavior,
  type SeqViewerContext,
} from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/spec.js';
import { shown } from './projection-fixtures.js';

type Shown = Awaited<ReturnType<typeof shown>>;

/** A hit on a residue of a unit (a letter of section `section`, block `block`). */
function seqHit(set: SeqSet, unit: number, index: number, section = 0, block = 0) {
  const hit = Hit.create();
  Object.assign(hit, { section, region: HitRegion.data, kind: HoverKind.seq, block, index, shown: index });
  Object.assign(hit.seq, { set, unit, index });
  Object.assign(hit.system, { kind: K.unit, id: unit, row: -1 });
  return hit;
}
const seqLoci = (set: SeqSet, unit: number, indices: number[]): SeqLoci =>
  ({
    kind: 'seq-loci',
    set,
    elements: [{ unit, indices: OrderedSet.ofSortedArray(Int32Array.from(indices)) }],
  }) as unknown as SeqLoci;
/** Per block label, the runs a channel lights (every section), after `sync`. */
function lit(t: Pick<Shown, 'set' | 'view' | 'label'>, ch: ChannelMarks, bit: number) {
  const out: Record<string, number[][]> = {};
  for (let s = 0; s < t.view.sections.count; s++) {
    const m = ch.marks(s)?.runs(bit);
    for (let b = 0; m && b < t.view.sections.blockCount[s]; b++) {
      const runs: number[][] = [];
      for (let r = m.from(b); r < m.to(b); r++) runs.push([m.start[r], m.stop[r]]);
      out[t.label(t.view.sections.blockStart[s] + b)] = runs;
    }
  }
  return out;
}
const keyed = (l: Loci | null) => {
  if (l?.kind !== 'seq-loci') return l;
  const values = (s: OrderedSet) => Array.from({ length: OrderedSet.size(s) }, (_, i) => OrderedSet.getAt(s, i));
  return l.elements.map((e) => [l.set.strings.get(l.set.units.key[e.unit]), ...values(e.indices)]);
};
/** Case 14's gene placed on a chromosome: gene 77 closes to residue 11, which lights 45, 76, 77 (VW5). */
async function geneOnChr() {
  const t = await shown(gene({ chromosome: true }));
  return { ...t, ch: new ChannelMarks(t.ctx), gu: t.unit('gene/unit'), pu: t.unit('protein/unit') };
}
const HIS12 = {
  gene: [
    [45, 46],
    [76, 78],
  ],
  protein: [[11, 12]],
  chr: [
    [145, 146],
    [176, 178],
  ],
};

describe('hover', () => {
  it("one owner at a time; a stale owner's clear does nothing; the version bumps per element only", async () => {
    const t = await geneOnChr();
    const h = t.ctx.hover.live;
    let bumps = 0;
    h.version.subscribe(() => bumps++);
    bumps = 0;
    const hit = seqHit(t.set, t.gu, 5);
    t.ctx.hover.setFromHit(hit, 'pointer');
    hit.shown = 5; // a move inside the element: the same element
    t.ctx.hover.setFromHit(hit, 'pointer');
    expect([bumps, h.current.owner, h.current.seq.index]).toEqual([1, 'pointer', 5]);
    hit.seq.index = hit.index = 6 as SeqIndex;
    t.ctx.hover.setFromHit(hit, 'pointer');
    expect(bumps).toBe(2);
    expect(keyed(t.ctx.hover.toLoci())).toEqual([['gene/unit', 6]]); // in the UTR: no residue contains it

    const loci = t.ctx.refs.toLoci('protein', { start: 3, stop: 4 }) as Loci;
    t.ctx.hover.setLoci(loci, '3d'); // a later owner replaces an earlier
    t.ctx.hover.clear('pointer'); // a late pointer-leave keeps the 3D hover
    expect([h.current.owner, h.current.kind, h.current.external, bumps]).toEqual(['3d', HoverKind.none, loci, 3]);
    expect(t.ctx.hover.toLoci()).toBe(loci); // external at `element`: as given
    t.ctx.hover.setLoci(null, 'pointer');
    expect(h.current.owner).toBe('3d');
    t.ctx.hover.clear('3d');
    expect([h.current.owner, h.current.external, bumps]).toEqual(['', null, 4]);
  });

  it("toLoci and size resolve once a version, whatever the hover's kind (VS-11)", async () => {
    const t = await geneOnChr();
    const hit = seqHit(t.set, t.gu, 77);
    t.ctx.hover.setFromHit(hit, 'pointer');
    const a = t.ctx.hover.toLoci();
    expect(t.ctx.hover.toLoci()).toBe(a);
    expect(t.ctx.hover.size()).toBe(t.ctx.hover.size());
    hit.seq.index = hit.index = 78 as SeqIndex;
    t.ctx.hover.setFromHit(hit, 'pointer');
    expect(t.ctx.hover.toLoci()).not.toBe(a);
  });

  it('toLoci is what is lit: a nucleotide closed to its residue, in their own units', async () => {
    const t = await geneOnChr();
    t.ctx.hover.setFromHit(seqHit(t.set, t.gu, 77), 'pointer');
    expect(keyed(t.ctx.hover.toLoci())).toEqual([
      ['gene/unit', 77],
      ['protein/unit', 11],
    ]);
    t.ctx.hover.setFromHit(seqHit(t.set, t.gu, 77), 'pointer', 'element');
    expect(keyed(t.ctx.hover.toLoci())).toEqual([['gene/unit', 77]]);
  });
});

describe('copies', () => {
  it('a molecule residue hovered lights and gives that residue in every alive unit of its molecule', async () => {
    const t = await shown(case7(), 'complex-each');
    const ch = new ChannelMarks(t.ctx);
    const a = t.unit('1brs/A');
    const hit = seqHit(t.set, a, 30);
    Object.assign(hit.system, { kind: K.molecule, id: t.set.units.molecule[a] });
    t.ctx.hover.setFromHit(hit, 'pointer');
    expect(keyed(t.ctx.hover.toLoci())).toEqual([
      ['1brs/A', 30],
      ['1brs/B', 30],
      ['1brs/C', 30],
    ]);
    expect(lit(t, ch, MarkBit.hover)).toEqual({ A: [[30, 31]], B: [[30, 31]], C: [[30, 31]], D: [], E: [], F: [] });
  });
});

describe("each channel's closure default", () => {
  it('pointer hover containing; external hover, highlights element; the select behavior containing', async () => {
    const t = await geneOnChr();
    t.ctx.hover.setFromHit(seqHit(t.set, t.gu, 77), 'pointer');
    expect(t.ctx.hover.live.current.closure).toBe('containing');
    expect(lit(t, t.ch, MarkBit.hover)).toEqual(HIS12);

    const nucleotide = seqLoci(t.set, t.gu, [77]);
    t.ctx.hover.setLoci(nucleotide, '3d');
    expect(t.ctx.hover.live.current.closure).toBe('element');
    // Element: no closure, but a nucleotide still projects onto the residue that contains it (through a path).
    const element = { gene: [[77, 78]], protein: [[11, 12]], chr: [[177, 178]] };
    expect(lit(t, t.ch, MarkBit.hover)).toEqual(element);

    t.ctx.highlight.set('host', nucleotide);
    expect(t.ctx.highlight.channels.value.get('host')?.closure).toBe('element');
    expect(lit(t, t.ch, MarkBit.highlight)).toEqual(element);
    t.ctx.highlight.set('link', nucleotide, { closure: 'containing', style: 'accent' });
    expect(lit(t, t.ch, MarkBit.accent)).toEqual(HIS12);

    // The hover and select behaviors' default closure (seqspec:viewer#behaviors): a click selects what the hover lit.
    expect([HoverBehavior.params.closure, SelectBehavior.params.closure]).toEqual(['containing', 'containing']);
  });

  it('protein, gene and chromosome light one another, both ways, from any of the three', async () => {
    const t = await geneOnChr();
    const chr = t.unit('chr/unit');
    for (const [unit, index, section] of [
      [t.gu, 77, 0],
      [t.pu, 11, 1],
      [chr, 177, 2],
    ]) {
      t.ctx.hover.setFromHit(seqHit(t.set, unit, index, section), 'pointer');
      expect(lit(t, t.ch, MarkBit.hover)).toEqual(HIS12);
    }
  });

  it("channel parity: hover, toLoci passed back, and a click's selection light the same in every section", async () => {
    const t = await geneOnChr();
    t.ctx.hover.setFromHit(seqHit(t.set, t.gu, 45), 'pointer');
    const hovered = lit(t, t.ch, MarkBit.hover);
    const loci = t.ctx.hover.toLoci() as Loci;
    t.ctx.highlight.set('host', loci); // closed already: element
    t.ctx.selection.set(t.ctx.selection.close(seqLoci(t.set, t.gu, [45]), 'containing'), 'pointer');
    expect(lit(t, t.ch, MarkBit.highlight)).toEqual(hovered);
    expect(lit(t, t.ch, MarkBit.selection)).toEqual(hovered);
    expect(hovered).toEqual(HIS12);
  });
});

describe('highlight', () => {
  it('channels are cleared independently', async () => {
    const t = await geneOnChr();
    const [a, b] = [t.ctx.refs.toLoci('gene') as Loci, t.ctx.refs.toLoci('protein') as Loci];
    t.ctx.highlight.set('host', a);
    t.ctx.highlight.set('links:l1', b, { style: 'accent' });
    expect([...t.ctx.highlight.channels.value.keys()]).toEqual(['host', 'links:l1']);
    const version = t.ctx.highlight.version.value;
    t.ctx.highlight.clear('host');
    expect([...t.ctx.highlight.channels.value]).toEqual([
      ['links:l1', { target: b, style: 'accent', closure: 'element', section: null }],
    ]);
    t.ctx.highlight.clear('host'); // nothing to clear
    t.ctx.highlight.set('links:l1', null);
    expect([t.ctx.highlight.channels.value.size, t.ctx.highlight.version.value]).toEqual([0, version + 2]);
  });

  it('a list of loci is refused at `set`; a channel whose carry throws is dropped alone, the others carried (VS-02)', async () => {
    const t = await geneOnChr();
    const [a, b] = [t.ctx.refs.toLoci('gene') as Loci, t.ctx.refs.toLoci('protein') as Loci];
    expect(() => t.ctx.highlight.set('mine', [a, b] as never)).toThrow(/union/);
    expect(t.ctx.highlight.channels.value.has('mine')).toBe(false);
    const failed: string[] = [];
    const h = new Highlight((m) => failed.push(m));
    h.set('good', a);
    h.channels.set(
      new Map(h.channels.value).set('bad', {
        ...(h.channels.value.get('good') as HighlightEntry),
        target: { kind: 'seq-loci', set: {} } as never,
      }),
    );
    h.set('after', b, { style: 'accent' });
    const next = { ...t.set, lineage: t.set.lineage } as SeqSet;
    h.update(next);
    expect([...h.channels.value.keys()]).toEqual(['good', 'after']);
    expect(failed).toEqual([expect.stringMatching(/"bad"/)]);
    h.dispose();
  });

  it('an item hover lights its rows; a pair from any position of either end (marks.any)', async () => {
    const t = await shown(case7(), 'complex-each');
    const ch = new ChannelMarks(t.ctx);
    const ann = t.set.annotations.find((a) => a && t.set.strings.get(a.key) === 'contacts');
    if (!ann) throw new Error('no contacts');
    const hit = Hit.create();
    Object.assign(hit, { section: 0, region: HitRegion.data, kind: HoverKind.item, block: 0, index: 0 });
    Object.assign(hit.itemLoc, { set: t.set, annotation: ann.id, item: 0 });
    Object.assign(hit.system, { kind: K.item, id: ann.id, row: 0 });
    t.ctx.hover.setFromHit(hit, 'pointer');
    const marks = ch.marks(0);
    const runs = lit(t, ch, MarkBit.hover);
    const ends = Object.entries(runs).filter(([, r]) => r.length > 0);
    expect(ends.length).toBe(2); // a contact's two ends: a barnase chain and a barstar chain
    const [label, [[s, e]]] = ends[0];
    const b = Object.keys(runs).indexOf(label);
    expect(marks?.any(b, s, e, MarkBit.hover)).toBe(true);
  });
});

describe('selection: close, then widen', () => {
  it("close crosses coordinate systems by containment; widen stays in each element's own", async () => {
    const t = await geneOnChr();
    const nt = seqLoci(t.set, t.gu, [77]);
    expect(t.ctx.selection.close(nt, 'element')).toBe(nt);
    const closed = t.ctx.selection.close(nt, 'containing');
    expect(keyed(closed)).toEqual([
      ['gene/unit', 77],
      ['protein/unit', 11],
    ]);
    const unit = t.ctx.selection.widen(closed, G.unit) as SeqLoci;
    expect(unit.elements.map((e) => [e.unit, OrderedSet.size(e.indices)])).toEqual([
      [t.gu, 120],
      [t.pu, 20],
    ]);
    // A drag's release closes the same way (the select behavior), and the verbs store loci as given.
    t.ctx.selection.set(closed, 'pointer');
    expect(lit(t, t.ch, MarkBit.selection)).toEqual(HIS12);
  });

  it('the verbs take Coordinates, converted through ctx.systems; unregistered keys are dropped', async () => {
    const t = await geneOnChr();
    t.ctx.selection.set(
      { systems: ['seqstar:protein/unit', 'nothing'], runs: Int32Array.of(0, 11, 12, 1, 0, 5) },
      'api',
    );
    expect(keyed(t.ctx.selection.current.value.seq)).toEqual([['protein/unit', 11]]);
    t.ctx.selection.add({ systems: ['seqstar:gene/unit'], runs: Int32Array.of(0, 3, 5) }, 'api');
    expect(keyed(t.ctx.selection.current.value.seq)).toEqual([
      ['gene/unit', 3, 4],
      ['protein/unit', 11],
    ]);
    t.ctx.selection.set({ systems: ['nothing'], runs: Int32Array.of(0, 0, 1) }, 'api');
    expect(t.ctx.selection.current.value.seq).toBe(null);
  });
});

describe("the channels' marks", () => {
  it("a channel's change rebuilds only its set; a layout change none; a view change all", async () => {
    const t = await geneOnChr();
    const ch = t.ch;
    const counts = () => ({ ...ch.rebuilds });
    ch.sync();
    expect(counts()).toEqual({ hover: 1, selection: 1, highlight: 1 });
    for (let i = 0; i < 10; i++) t.ctx.hover.setFromHit(seqHit(t.set, t.gu, 40 + i), 'pointer');
    ch.sync();
    ch.sync();
    expect(counts()).toEqual({ hover: 2, selection: 1, highlight: 1 });
    t.ctx.selection.set(seqLoci(t.set, t.gu, [3]), 'api');
    ch.sync();
    t.ctx.highlight.set('host', seqLoci(t.set, t.pu, [2]));
    ch.sync();
    expect(counts()).toEqual({ hover: 2, selection: 2, highlight: 2 });
    // The layout: render props, sizes and a viewport's zoom rebuild sections, never the channels.
    const layout = t.ctx.layout.get('v.main');
    const version = layout?.version.value ?? 0;
    t.ctx.setRenderProps({ cellMin: 3 } as never);
    for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(640, 300, 2);
    t.scheduler.advance(16); // the frame that rebuilds the layouts
    ch.sync();
    expect(layout?.version.value).toBeGreaterThan(version); // the positive control: the layout was rebuilt
    expect(counts()).toEqual({ hover: 2, selection: 2, highlight: 2 });
    // A view change: every channel again, once.
    await t.ctx.applyPatch({ format: 'seq_view_spec_patch', ops: [{ op: 'remove', ref: 'v.chr' }] } as never);
    ch.sync();
    expect(counts()).toEqual({ hover: 3, selection: 3, highlight: 3 });
    expect(lit(t, ch, MarkBit.highlight)).toEqual({ gene: [[18, 21]], protein: [[2, 3]] });
  });

  it("a rows section: row residues and cells (a gap's too) are rectangles; the hover lights the cell's column", async () => {
    const t = await shown(rows());
    const ch = new ChannelMarks(t.ctx);
    const a = 0;
    t.ctx.selection.set(
      {
        kind: 'aligned-loci',
        set: t.set,
        alignment: a,
        elements: [{ row: 0, columns: OrderedSet.ofRange(1, 4) }],
      } as unknown as Loci,
      'api',
    );
    t.ctx.highlight.set('host', {
      kind: 'row-loci',
      set: t.set,
      alignment: a,
      elements: [{ row: 1, indices: OrderedSet.ofRange(0, 2) }],
    } as unknown as Loci);
    const m = ch.marks(0);
    const rects = (bit: number) =>
      Array.from({ length: m?.rows.count ?? 0 }, (_, r) => r)
        .filter((r) => ((m?.rows.bit[r] ?? 0) & ~MarkBit.row) === bit)
        .map((r) => [m?.rows.rowStart[r], m?.rows.rowStop[r], m?.rows.start[r], m?.rows.stop[r]]);
    // Row 0 "MK-VL": the cells [1, 4) with the gap at 2 are one rectangle (IN1 lit the gap as its column only).
    expect(rects(MarkBit.selection)).toEqual([[0, 1, 1, 4]]);
    expect(rects(MarkBit.highlight)).toEqual([[1, 2, 0, 2]]);
    // The row residues light their row whole (`MarkBit.row`, CS1b); without a rows layer, columns too.
    expect(m?.rows.bit[1]).toBe(MarkBit.highlight | MarkBit.row);
    expect([lit(t, ch, MarkBit.selection).a, lit(t, ch, MarkBit.highlight).a]).toEqual([[[1, 4]], [[0, 2]]]);
    const hit = Hit.create();
    Object.assign(hit, { section: 0, region: HitRegion.data, kind: HoverKind.aligned, block: 0, index: 3 });
    Object.assign(hit.aligned, { set: t.set, alignment: a, row: 0, column: 3 });
    Object.assign(hit.system, { kind: K.row, id: a, row: 0 });
    t.ctx.hover.setFromHit(hit, 'pointer');
    expect(lit(t, ch, MarkBit.hover).a).toEqual([[3, 4]]);
    expect(lit(t, ch, MarkBit.hover).x).toEqual([[4, 5]]); // row 0 is bound to x at offset 2 (V is x's 4)
  });
});

describe('coordinate systems at the edge', () => {
  const external = (key: string, members: string[] | null = null): SystemDescriptor => ({
    key,
    kind: members ? 'complex' : 'external',
    start: 0,
    stop: 200,
    refs: {},
    members,
    structure: null,
    numberings: [],
    owner: 'molstar',
  });

  it('own descriptors by key, built once per system and lineage (counter); coordinatesOf lists a complex', async () => {
    const t = await shown(case7(), 'complex-each');
    const sys = t.ctx.systems;
    const a = sys.descriptor('seqstar:1brs/A');
    expect(a).toMatchObject({ kind: 'unit', refs: { unit: '1brs/A', complex: '1brs' }, owner: 'seqstar' });
    expect(sys.descriptor('seqstar:1brs')).toMatchObject({
      kind: 'complex',
      members: expect.arrayContaining(['seqstar:1brs/A']),
    });
    const built = sys.built;
    for (let i = 0; i < 1000; i++) sys.descriptor('seqstar:1brs/A');
    expect(sys.descriptorOf({ kind: K.unit, id: t.unit('1brs/A'), row: -1 })).toBe(a);
    expect(sys.built).toBe(built);
    await t.ctx.applyPatch({ format: 'seq_view_spec_patch', ops: [{ op: 'remove', ref: 'complex-one' }] } as never);
    expect(sys.descriptor('seqstar:1brs/A')).toBe(a); // the same set's tables: not built again
    expect(sys.descriptor('seqstar:nothing')).toBe(null);

    const loci = {
      kind: 'seq-loci',
      set: t.set,
      elements: [
        { unit: t.unit('1brs/A'), indices: OrderedSet.ofRange(3, 6) },
        { unit: t.unit('1brs/D'), indices: OrderedSet.ofRange(0, 2) },
      ],
    } as unknown as Loci;
    const c = sys.coordinatesOf(loci);
    expect(c.systems).toEqual(['seqstar:1brs/A', 'seqstar:1brs', 'seqstar:1brs/D']);
    expect([...c.runs]).toEqual([0, 3, 6, 2, 0, 2]);
    // Passed back as a request, they light what the loci lights.
    t.ctx.highlight.set('host', c);
    const ch = new ChannelMarks(t.ctx);
    expect(lit(t, ch, MarkBit.highlight)).toMatchObject({ A: [[3, 6]], D: [[0, 2]], B: [] });
  });

  it('a Mol*-like structure registered: a highlight and a hover in it, a complex across chains, unload and reload', async () => {
    const t = await shown(case7(), 'complex-each');
    const sys = t.ctx.systems;
    const ch = new ChannelMarks(t.ctx);
    const unit = (key: string) => sys.system(`seqstar:${key}`) as { kind: number; id: number; row: number };
    const mapA = (s: { kind: number; id: number; row: number }) => [
      TranslateMaps.bySegments(s, unit('1brs/A'), Int32Array.of(5, 0, 105)),
    ];
    sys.register(external('molstar:3d/A'), mapA);
    sys.register(external('molstar:3d/D'), (s) => [
      TranslateMaps.bySegments(s, unit('1brs/D'), Int32Array.of(0, 0, 89)),
    ]);
    sys.register(external('molstar:3d', ['molstar:3d/A', 'molstar:3d/D']));
    const across = { systems: ['molstar:3d/A', 'molstar:3d/D'], runs: Int32Array.of(0, 10, 20, 1, 0, 3) };
    t.ctx.highlight.set('host', across);
    expect(lit(t, ch, MarkBit.highlight)).toEqual({ A: [[5, 15]], B: [], C: [], D: [[0, 3]], E: [], F: [] });
    t.ctx.hover.setLoci({ systems: ['molstar:3d/A'], runs: Int32Array.of(0, 30, 31) }, '3d');
    expect(lit(t, ch, MarkBit.hover).A).toEqual([[25, 26]]);
    expect(keyed(t.ctx.hover.toLoci())).toEqual([['1brs/A', 25]]);
    // Unloaded: nothing lights, the targets are kept by key; loaded again: they light again.
    sys.unregister('molstar:3d/A');
    expect(lit(t, ch, MarkBit.highlight)).toMatchObject({ A: [], D: [[0, 3]] });
    sys.register(external('molstar:3d/A'), mapA);
    expect(lit(t, ch, MarkBit.highlight)).toMatchObject({ A: [[5, 15]], D: [[0, 3]] });
  });

  it('a hover in a system not registered yet: toLoci and size resolve once it is, and again once it goes (LB-07)', async () => {
    const t = await shown(case7(), 'complex-each');
    const sys = t.ctx.systems;
    const unit = (key: string) => sys.system(`seqstar:${key}`) as { kind: number; id: number; row: number };
    const mapA = (s: { kind: number; id: number; row: number }) => [
      TranslateMaps.bySegments(s, unit('1brs/A'), Int32Array.of(5, 0, 105)),
    ];
    t.ctx.hover.setLoci({ systems: ['molstar:3d/A'], runs: Int32Array.of(0, 30, 31) }, '3d');
    expect([t.ctx.hover.toLoci(), t.ctx.hover.size()]).toEqual([null, null]);
    sys.register(external('molstar:3d/A'), mapA);
    expect(keyed(t.ctx.hover.toLoci())).toEqual([['1brs/A', 25]]);
    expect(t.ctx.hover.size()).toMatchObject({ length: 1 });
    sys.unregister('molstar:3d/A');
    expect(t.ctx.hover.toLoci()).toBeNull();
  });
});

describe('held targets follow their set', () => {
  const mol = (key: string, residues: string): MoleculeInput => ({ key, residues, ownUnit: true });
  const [p53, q] = [mol('p53', 'MEEPQSDPSVEPPLSQETFSDLWKLL'), mol('q', 'MKVLATGHEK')];
  const build = (ms: MoleculeInput[]) => SeqSets.build({ molecules: ms }).set;
  const units = (ctx: SeqViewerContext, l: Loci | null | undefined) =>
    l?.kind === 'seq-loci'
      ? l.elements.map((e) => [ctx.graph.assembled.value?.set.strings.get(l.set.units.key[e.unit]), e.unit])
      : l;

  it("a Coordinates hover resolves on the new set, not the old one's loci kept (LB-07)", async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    await ctx.setData({ set: build([p53, q]) });
    ctx.hover.setLoci({ systems: ['seqstar:q/unit'], runs: Int32Array.of(0, 1, 2) }, '3d');
    const before = ctx.hover.toLoci();
    expect(before?.set).toBe(ctx.graph.assembled.value?.set);
    await ctx.setData({ set: build([mol('r', 'MKV'), q]) });
    const after = ctx.hover.toLoci();
    expect(after).not.toBe(before);
    expect(after?.set).toBe(ctx.graph.assembled.value?.set);
    expect(units(ctx, after as Loci)).toEqual([['q/unit', expect.any(Number)]]);
  });

  it('a host highlight and an external hover remap across a new lineage; Coordinates are held by key', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    await ctx.setData({ set: build([p53, q]) });
    ctx.highlight.set('host', ctx.refs.toLoci('q/unit', { start: 2, stop: 4 }));
    ctx.hover.setLoci(ctx.refs.toLoci('q/unit', { start: 1, stop: 2 }), '3d');
    const coords = { systems: ['seqstar:q/unit'], runs: Int32Array.of(0, 0, 1) };
    ctx.highlight.set('coords', coords);
    await ctx.setData({ set: build([mol('r', 'MKV'), q]) }); // a new lineage: q's unit has another id
    const set = ctx.graph.assembled.value?.set;
    const host = ctx.highlight.channels.value.get('host')?.target as Loci;
    expect([host.set === set, units(ctx, host)]).toEqual([true, [['q/unit', 1]]]);
    const hovered = ctx.hover.live.current.external as Loci;
    expect([hovered.set === set, units(ctx, hovered), ctx.hover.live.current.owner]).toEqual([
      true,
      [['q/unit', 1]],
      '3d',
    ]);
    expect(ctx.highlight.channels.value.get('coords')?.target).toBe(coords);
    // Nothing left of a target: the channel goes; a null set changes nothing.
    ctx.highlight.set('p', ctx.refs.toLoci('r/unit'));
    await ctx.setData({ set: build([q]) });
    expect([...ctx.highlight.channels.value.keys()]).toEqual(['host', 'coords']);
  });
});
