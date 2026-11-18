/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { type SeqData, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { elements, rangesLocation } from '@molstar/seqx-core/query';
import { type PassName, recordDrawLists } from '@molstar/seqx-core/render';
import { nodeAddress } from '@molstar/seqx-core/resolve';
import { afterEach, describe, expect, it } from 'vitest';
import type { EditingBehavior } from '../behaviors/editing.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

// The selection's grammar through the viewer's input (ux-foundations' section 8): click sets (clicking the only
// selected element clears it), Shift+click extends from the anchor, Cmd/Ctrl+click toggles, a click on empty space or
// Esc clears, a drag on data selects a range (whole items when it began on one); the keyboard cursor; the selection and
// the highlights drawn, projected onto the gene's section, in the canvas's lists and an SVG export alike; the
// selection's dynamic nodes held.

const P = '/data/p';
const G = '/data/g';
const ProteinSection = '/views/main/protein';
const GeneSection = '/views/main/gene';

/** 10 residues on a 40-nucleotide gene (one exon, 5 to 35); domains A [0, 3) and B [5, 9) on the protein. */
function spec(): Node {
  const b = createSeqSpec({ title: 'selection' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [5], end: [35] }, strand: '+' });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: new Int32Array([0, 5]), end: new Int32Array([3, 9]), id: ['A', 'B'] },
  });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p });
    protein.track('letters').layer('l', { type: 'letters' });
    protein.track('domains').layer('d', { type: 'ranges', data: domains });
    v.section('gene', { of: g, follow: protein }).track('letters').layer('l', { type: 'letters' });
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

async function make(root: Node = spec(), config: readonly (readonly [any, unknown])[] = []) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000], ...config],
  });
  viewers.push(viewer);
  viewer.input.mac = false;
  await viewer.load(root);
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

/** The middle of a position's cell on a track (by index in its section), in canvas px. */
function at(v: SeqViewer, section: number, track: number, position: number, lane = 0): [number, number] {
  const s = v.frame!.layout.sections[section];
  const box = s.tracks[track].box;
  const x = s.axis!.x(position) + (s.axis!.x(position + 1) - s.axis!.x(position)) / 2;
  return [x, track === 1 && section === 0 ? box.y + 2 + lane * 12 + 4 : box.y + box.height / 2];
}

const click = (
  v: SeqViewer,
  [x, y]: [number, number],
  e: { shiftKey?: boolean; ctrlKey?: boolean; detail?: number } = {},
) => {
  v.input.move(x, y, e);
  v.input.down(x, y, e);
  v.input.up();
  return v.input.click(x, y, e);
};

const ranges = (v: SeqViewer, system = P) => [
  ...(v.selection.value.value.shown.parts.find((p) => p.system === system)?.ranges ?? []),
];

describe('clicks (ux-foundations’ section 8)', () => {
  it('click sets the selection to the cell; clicking it again clears it; empty space clears', async () => {
    const v = await make();
    expect(click(v, at(v, 0, 0, 2))).toBe(true);
    expect(ranges(v)).toEqual([2, 3]);
    expect(v.selection.cursor.value).toEqual({ section: ProteinSection, system: P, position: 2 });
    click(v, at(v, 0, 0, 2));
    expect(ranges(v)).toEqual([]);
    click(v, at(v, 0, 0, 4));
    // the label column is empty space
    click(v, [40, at(v, 0, 0, 0)[1]]);
    expect(ranges(v)).toEqual([]);
  });

  it('an item is selected by its key, drawn by its extent; Shift+click extends from it in one system', async () => {
    const v = await make();
    click(v, at(v, 0, 1, 6));
    const sel = v.selection.value.value;
    expect(sel.elements.keys.map((k) => k.key)).toEqual(['B']);
    expect(ranges(v)).toEqual([5, 9]);
    click(v, at(v, 0, 0, 2), { shiftKey: true });
    expect(ranges(v)).toEqual([2, 9]);
    expect(v.selection.value.value.elements.keys).toEqual([]);
    // extends from the same anchor replace one another
    click(v, at(v, 0, 0, 9), { shiftKey: true });
    expect(ranges(v)).toEqual([5, 10]);
  });

  it('Cmd/Ctrl+click toggles a cell or an item, keeping the rest', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 1), { ctrlKey: true });
    click(v, at(v, 0, 0, 3), { ctrlKey: true });
    expect(ranges(v)).toEqual([1, 2, 3, 4]);
    click(v, at(v, 0, 1, 0), { ctrlKey: true });
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['A']);
    expect(ranges(v)).toEqual([0, 4]);
    click(v, at(v, 0, 0, 1), { ctrlKey: true });
    click(v, at(v, 0, 1, 0), { ctrlKey: true });
    expect(ranges(v)).toEqual([3, 4]);
    expect(v.selection.value.value.elements.keys).toEqual([]);
  });

  it('the second click of a double-click selects nothing (the double-click focuses)', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 2));
    expect(click(v, at(v, 0, 0, 2), { detail: 2 })).toBe(false);
    expect(ranges(v)).toEqual([2, 3]);
  });
});

describe('drag-select', () => {
  it('a drag on data selects the range it covers, read past the data too; the click after it is none', async () => {
    const v = await make();
    const [x0, y] = at(v, 0, 0, 1);
    const start = v.view.get(ProteinSection)!.start;
    expect(v.input.down(x0, y)).toBe(true);
    // under the slop: nothing yet
    v.input.move(x0 + 2, y);
    expect(ranges(v)).toEqual([]);
    v.input.move(at(v, 0, 0, 4)[0], y);
    expect(ranges(v)).toEqual([1, 5]);
    // back past its start
    v.input.move(at(v, 0, 0, 0)[0], y);
    expect(ranges(v)).toEqual([0, 2]);
    // past the data's right edge: the last position
    v.input.move(10_000, y);
    expect(ranges(v)).toEqual([1, 10]);
    v.input.up();
    expect(v.input.click(10_000, y)).toBe(false);
    expect(ranges(v)).toEqual([1, 10]);
    expect(v.view.get(ProteinSection)!.start).toBe(start);
    // an extend from there starts at the drag's anchor
    click(v, at(v, 0, 0, 3), { shiftKey: true });
    expect(ranges(v)).toEqual([1, 4]);
  });

  it('begun on an item, it snaps to whole items', async () => {
    const v = await make();
    const [x0, y] = at(v, 0, 1, 6);
    v.input.down(x0, y);
    v.input.move(at(v, 0, 1, 1)[0], y);
    v.input.up();
    expect(ranges(v)).toEqual([0, 9]);
  });
});

describe('the keyboard cursor', () => {
  const key = (v: SeqViewer, k: string, e: { shiftKey?: boolean } = {}) => v.input.key({ key: k, ...e });

  it('the step keys place it, move it in display order, and Shift extends the selection from it', async () => {
    const v = await make();
    expect(key(v, 'ArrowRight')).toBe(true);
    expect(v.selection.cursor.value).toEqual({ section: ProteinSection, system: P, position: 0 });
    key(v, 'ArrowRight');
    key(v, 'ArrowRight');
    expect(v.selection.cursor.value?.position).toBe(2);
    key(v, 'ArrowLeft');
    expect(v.selection.cursor.value?.position).toBe(1);
    // clamped at the axis's start
    key(v, 'ArrowLeft');
    key(v, 'ArrowLeft');
    expect(v.selection.cursor.value?.position).toBe(0);
    key(v, 'ArrowRight', { shiftKey: true });
    key(v, 'ArrowRight', { shiftKey: true });
    expect(ranges(v)).toEqual([0, 3]);
    expect(v.selection.cursor.value?.position).toBe(2);
    // drawn on its section's tracks
    frames!.flush();
    expect(v.frame!.lists.passes.cursor.rects.n).toBe(2);
  });

  it('Space (a tap) toggles its cell; Enter focuses; down and up carry it through the codon map', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 3));
    v.input.setSpace(true);
    expect(v.input.spaceUp()).toBe(true);
    expect(ranges(v)).toEqual([]);
    v.input.setSpace(true);
    v.input.spaceUp();
    expect(ranges(v)).toEqual([3, 4]);
    expect(key(v, 'Enter')).toBe(true);
    const vp = v.view.get(ProteinSection)!;
    expect([vp.start, vp.end]).toEqual([3, 4]);
    v.view.fit(ProteinSection);
    // residue 3 is nucleotides 14 to 17 on the gene
    key(v, 'ArrowDown');
    expect(v.selection.cursor.value).toEqual({ section: GeneSection, system: G, position: 14 });
    key(v, 'ArrowUp');
    expect(v.selection.cursor.value).toEqual({ section: ProteinSection, system: P, position: 3 });
  });

  it('a cursor moved out of view pans its section to show it', async () => {
    const v = await make();
    v.view.setRange(ProteinSection, 0, 4);
    v.selection.setCursor({ section: ProteinSection, system: P, position: 3 });
    key(v, 'ArrowRight');
    const vp = v.view.get(ProteinSection)!;
    expect(vp.end).toBe(5);
    expect(vp.end - vp.start).toBeCloseTo(4);
  });

  it('Esc cancels a gesture’s preview first, then clears the selection', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 2));
    (v.behaviors.get('editing') as EditingBehavior).preview({
      ops: [{ op: 'set', at: ['views', 'main', 'protein', 'letters'], key: 'height', value: 30 }],
    });
    expect(key(v, 'Escape')).toBe(true);
    expect(v.spec.preview.value).toBeUndefined();
    expect(ranges(v)).toEqual([2, 3]);
    expect(key(v, 'Escape')).toBe(true);
    expect(ranges(v)).toEqual([]);
    // nothing left to do: the key is the page's
    expect(key(v, 'Escape')).toBe(false);
  });
});

/** The rects of a pass within a section's band, as recorded. */
function band(v: SeqViewer, pass: PassName, section: number): string[] {
  const s = v.frame!.layout.sections[section];
  return recordDrawLists(v.frame!.lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .filter((l) => {
      const y = Number(l.split(' ')[2]);
      return y >= s.y && y < s.y + s.height;
    });
}

describe('drawn, projected onto linked sections', () => {
  it('a selection on the protein shows on the gene, lighter; the SVG export draws it alike', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 2));
    frames!.flush();
    const protein = band(v, 'selection', 0);
    const gene = band(v, 'selection', 1);
    expect(protein.some((l) => l.includes('#8f4c062e'))).toBe(true);
    expect(gene.some((l) => l.includes('#8f4c0617'))).toBe(true);
    // the export of the same spec: the selection's Location, projected by the export itself
    const svg = await v.renderSvg();
    expect(recordDrawLists(svg.lists, ['selection'])).toBe(recordDrawLists(v.frame!.lists, ['selection']));
  });

  it('a hover on the protein is outlined on the gene’s codon, lighter, and not again on its own section', async () => {
    const v = await make();
    v.input.move(...at(v, 0, 0, 2));
    frames!.flush();
    const gene = band(v, 'hover', 1);
    expect(gene.length).toBeGreaterThan(0);
    expect(gene.every((l) => l.includes('stroke #2d271e80'))).toBe(true);
    expect(band(v, 'hover', 0).every((l) => l.includes('stroke #2d271eff'))).toBe(true);
  });

  it('highlight channels are drawn as tints in their colors, projected, in the canvas and the export', async () => {
    const v = await make();
    v.addresses.highlightAddress('/data/p/(range :from 2 :to 3)');
    frames!.flush();
    // `link`: the default color, `token:cat-5` (#7aa3a6) at 0.3
    expect(band(v, 'highlight', 0).some((l) => l.includes('#7aa3a64d'))).toBe(true);
    expect(band(v, 'highlight', 1).some((l) => l.includes('#7aa3a64d'))).toBe(true);
    const svg = await v.renderSvg();
    expect(recordDrawLists(svg.lists, ['highlight'])).toBe(recordDrawLists(v.frame!.lists, ['highlight']));
    v.addresses.unhover();
    v.highlight.clear('link');
    frames!.flush();
    expect(v.frame!.lists.passes.highlight.rects.n).toBe(0);
  });

  it('the tint under the data (a config item): drawn in the background pass, the base redrawn on a change', async () => {
    const { SeqCoreConfig } = await import('@molstar/seqx-core/context');
    const v = await make(spec(), [[SeqCoreConfig.render, { emphasis: { selectionPlace: 'under' } }]]);
    click(v, at(v, 0, 0, 2));
    frames!.flush();
    expect(v.frame!.base).toBe(true);
    expect(band(v, 'background', 0).some((l) => l.includes('#8f4c062e'))).toBe(true);
    expect(band(v, 'selection', 0).every((l) => l.includes('stroke'))).toBe(true);
  });
});

describe('the selection holds its nodes', () => {
  const Genes = '/data/genes';

  /** A gene region whose one answer brings a transcript (a dynamic node), drawn by a template. */
  function dynamicSpec(): Node {
    const transcript: Node = {
      kind: 'ranges',
      params: { target: '/data/gene', rows: { start: [10], end: [40], id: ['T1'] } },
    };
    const data: SeqData = {
      coverage: { systems: [{ name: 'g', length: 100 }], shapes: [{ system: 'g', shape: 'letters' }], levels: [1] },
      pairs: [
        {
          request: { system: 'g', start: 0, end: 100, level: 1, shapes: ['letters'] },
          answer: {
            ok: true,
            tables: { letters: { form: 'letters', start: 0, letters: 'A'.repeat(100) } },
            nodes: [{ key: 'T1', node: transcript }],
          },
        },
      ],
    };
    const b = createSeqSpec();
    const src = b.data.source('genes', { format: 'seq', value: seqToValue(data) });
    const gene = b.data.sequence('gene', { alphabet: 'dna', source: src, name: 'g' });
    b.views.view('main', (v) =>
      v
        .section('gene', { of: gene })
        .track('t', { for_each: { source: src, node_kind: 'ranges', depth: 1 } })
        .layer('ranges', { type: 'ranges' }),
    );
    return b.build();
  }

  it('a selected transcript is held while it is selected, and let go when it is not', async () => {
    const v = await make(dynamicSpec());
    const owner = nodeAddress(Genes, 'T1');
    v.selection.set({
      elements: elements([{ owner, key: 'T1' }]),
      location: rangesLocation('/data/gene', [10, 40]),
      shown: rangesLocation('/data/gene', [10, 40]),
    });
    expect(v.store.sessionHolds('selection')).toEqual([`${Genes}#T1`]);
    v.selection.clear();
    expect(v.store.sessionHolds('selection')).toEqual([]);
  });
});
