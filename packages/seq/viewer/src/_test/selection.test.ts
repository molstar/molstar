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
import { CanvasBackend } from '../canvas/backend.js';
import { FakeCanvas } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

// The selection's grammar through the viewer's input (the interaction grammar): click sets (clicking the only
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

describe('clicks (the interaction grammar)', () => {
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

  it('a double-click on the only selected item keeps it: its first click’s clear is undone', async () => {
    const v = await make();
    const item = at(v, 0, 1, 6);
    click(v, item);
    // the double-click's first click clears it, its second (detail 2) gives it back
    click(v, item, { detail: 1 });
    expect(ranges(v)).toEqual([]);
    click(v, item, { detail: 2 });
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['B']);
    // the same through the platform's dblclick
    click(v, item, { detail: 1 });
    v.input.doubleClick(...item);
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['B']);
    // a single click still clears it
    click(v, item, { detail: 1 });
    click(v, at(v, 0, 0, 0), { detail: 1 });
    expect(ranges(v)).toEqual([0, 1]);
  });
});

type Keys = { shiftKey?: boolean; ctrlKey?: boolean; altKey?: boolean };

/** A drag on the protein's letters from one position to another, as the canvas gets it (the button left down). */
function dragOver(v: SeqViewer, from: number, to: number, e: Keys = {}) {
  const [x0, y] = at(v, 0, 0, from);
  v.input.move(x0, y, e);
  v.input.down(x0, y, e);
  for (let i = 1; i <= 4; i++) v.input.move(x0 + ((at(v, 0, 0, to)[0] - x0) * i) / 4, y, e);
}

function dragged(v: SeqViewer, from: number, to: number, e: Keys = {}) {
  dragOver(v, from, to, e);
  v.input.up();
}

describe('several ranges in one selection (the interaction grammar)', () => {
  it('Cmd/Ctrl+drag adds a range, keeping the rest; a plain drag replaces them all', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 1));
    dragged(v, 4, 6, { ctrlKey: true });
    expect(ranges(v)).toEqual([1, 2, 4, 7]);
    // and an item stays an item
    click(v, at(v, 0, 1, 0), { ctrlKey: true });
    dragged(v, 8, 9, { ctrlKey: true });
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['A']);
    expect(ranges(v)).toEqual([0, 3, 4, 7, 8, 10]);
    dragged(v, 2, 3);
    expect(ranges(v)).toEqual([2, 4]);
    expect(v.selection.value.value.elements.keys).toEqual([]);
  });

  it('Alt/Option+drag removes a range; an item it cuts becomes the positions left of it', async () => {
    const v = await make();
    dragged(v, 0, 8);
    dragged(v, 3, 4, { altKey: true });
    expect(ranges(v)).toEqual([0, 3, 5, 9]);
    // domain B (5 to 9) selected by its key, cut in its middle: its positions left; domain A, missed, stays a key
    click(v, at(v, 0, 1, 6));
    click(v, at(v, 0, 1, 0), { ctrlKey: true });
    dragged(v, 6, 7, { altKey: true });
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['A']);
    expect(ranges(v)).toEqual([0, 3, 5, 6, 8, 9]);
  });

  it('Shift+click and Shift+drag extend only the range the anchor is in, keeping the others', async () => {
    const v = await make();
    dragged(v, 0, 2);
    dragged(v, 6, 7, { ctrlKey: true });
    expect(ranges(v)).toEqual([0, 3, 6, 8]);
    // the anchor is 6, where the added range began
    click(v, at(v, 0, 0, 9), { shiftKey: true });
    expect(ranges(v)).toEqual([0, 3, 6, 10]);
    // extends from it replace one another, the first range kept
    click(v, at(v, 0, 0, 5), { shiftKey: true });
    expect(ranges(v)).toEqual([0, 3, 5, 7]);
    // a Shift+drag grows the anchor's range wherever it starts
    dragged(v, 9, 8, { shiftKey: true });
    expect(ranges(v)).toEqual([0, 3, 6, 9]);
    // Cmd/Ctrl+click toggles, and the anchor moves there
    click(v, at(v, 0, 0, 1), { ctrlKey: true });
    expect(ranges(v)).toEqual([0, 1, 2, 3, 6, 9]);
    click(v, at(v, 0, 0, 0), { shiftKey: true });
    expect(ranges(v)).toEqual([0, 3, 6, 9]);
    click(v, at(v, 0, 0, 4), { shiftKey: true });
    expect(ranges(v)).toEqual([0, 5, 6, 9]);
  });

  it('Esc during an additive or a subtractive drag gives back the selection it began with', async () => {
    const v = await make();
    dragged(v, 0, 2);
    dragged(v, 6, 7, { ctrlKey: true });
    for (const e of [{ ctrlKey: true }, { altKey: true }, { shiftKey: true }]) {
      dragOver(v, 1, 8, e);
      expect(ranges(v)).not.toEqual([0, 3, 6, 8]);
      expect(v.input.key({ key: 'Escape' })).toBe(true);
      v.input.up();
      expect(ranges(v)).toEqual([0, 3, 6, 8]);
    }
    // the anchor came back too: an extend still grows the range at 6
    click(v, at(v, 0, 0, 9), { shiftKey: true });
    expect(ranges(v)).toEqual([0, 3, 6, 10]);
  });

  it('Shift+arrows extend the anchor’s range; Cmd/Ctrl+Shift+arrows start a new range at the cursor', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 1));
    v.input.key({ key: 'ArrowRight', shiftKey: true });
    expect(ranges(v)).toEqual([1, 3]);
    v.selection.setCursor({ section: ProteinSection, system: P, position: 6 });
    v.input.key({ key: 'ArrowRight', shiftKey: true, ctrlKey: true });
    expect(ranges(v)).toEqual([1, 3, 6, 8]);
    v.input.key({ key: 'ArrowRight', shiftKey: true });
    expect(ranges(v)).toEqual([1, 3, 6, 9]);
    v.input.key({ key: 'ArrowLeft', shiftKey: true });
    v.input.key({ key: 'ArrowLeft', shiftKey: true });
    v.input.key({ key: 'ArrowLeft', shiftKey: true });
    expect(ranges(v)).toEqual([1, 3, 5, 7]);
  });

  it('selection.add and selection.subtract over a range or an address, for the keyboard and a host', async () => {
    const v = await make();
    expect(await v.execute('selection.add', { range: { system: P, start: 2, end: 4 } })).toBe('Added');
    expect(await v.execute('selection.add', { address: '/data/domains/(item :id B)' })).toBe('Added');
    expect(ranges(v)).toEqual([2, 4, 5, 9]);
    expect(await v.execute('selection.subtract', { address: '/data/p/(range :from 7 :to 7)' })).toBe('Removed');
    expect(ranges(v)).toEqual([2, 4, 5, 6, 7, 9]);
    const add = v.commands.find((c) => c.id === 'selection.add')!;
    expect(add.when!({} as never)).toBe('No range, address or position given');
  });

  it('with selection.extendScope: selection, an extend makes the whole selection the anchor’s range', async () => {
    const v = await make(spec(), [[viewerItem('selection.extendScope'), 'selection']]);
    dragged(v, 0, 2);
    dragged(v, 6, 7, { ctrlKey: true });
    click(v, at(v, 0, 0, 9), { shiftKey: true });
    expect(ranges(v)).toEqual([6, 10]);
  });
});

describe('the review’s smaller points', () => {
  it('an item selected by a link clears on the first click on it (compared without its address)', async () => {
    const v = await make();
    expect(v.addresses.select('/data/domains/(item :id B)').done).toBe(true);
    click(v, at(v, 0, 1, 6));
    expect(ranges(v)).toEqual([]);
  });

  it('Esc during a drag-select gives back the selection before it', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 7));
    const [x0, y] = at(v, 0, 0, 1);
    v.input.down(x0, y);
    v.input.move(at(v, 0, 0, 4)[0], y);
    expect(ranges(v)).toEqual([1, 5]);
    expect(v.input.key({ key: 'Escape' })).toBe(true);
    expect(ranges(v)).toEqual([7, 8]);
    v.input.up();
    expect(v.input.click(at(v, 0, 0, 4)[0], y)).toBe(false);
    expect(ranges(v)).toEqual([7, 8]);
  });

  it('a viewer without the select behavior has no selection actions: its bindings to them are reported', async () => {
    const { SelectProvider } = await import('../behaviors/providers.js');
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      ui: undefined,
      behaviors: base.behaviors.filter((b) => b.provider !== SelectProvider),
    });
    viewers.push(v);
    expect(v.bindings.issues.map((i) => i.code)).toContain('seq.viewer.binding.action');
    expect(v.bindings.resolve('click', {}, 'cell')).toBeUndefined();
    expect((await make()).bindings.issues).toEqual([]);
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
    expect(v.selection.cursor.value).toEqual({ section: ProteinSection, system: P, position: 0, keyboard: true });
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
    // drawn once over its section's column
    frames!.flush();
    expect(v.frame!.lists.passes.cursor.rects.n).toBe(1);
  });

  it('a click places the cursor unseen; the first key shows it where it is (focus-visible)', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 2));
    frames!.flush();
    expect(v.selection.cursor.value?.keyboard).toBeUndefined();
    expect(v.frame!.lists.passes.cursor.rects.n).toBe(0);
    key(v, 'ArrowRight');
    frames!.flush();
    expect(v.selection.cursor.value).toMatchObject({ position: 3, keyboard: true });
    expect(v.frame!.lists.passes.cursor.rects.n).toBe(1);
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
    expect(v.selection.cursor.value).toEqual({ section: GeneSection, system: G, position: 14, keyboard: true });
    key(v, 'ArrowUp');
    expect(v.selection.cursor.value).toEqual({ section: ProteinSection, system: P, position: 3, keyboard: true });
  });

  it('Enter on a follower’s cursor focuses its leader, as a focus on an address does', async () => {
    const v = await make();
    // nucleotide 15 is in residue 3's codon (14 to 17)
    v.selection.setCursor({ section: GeneSection, system: G, position: 15 });
    expect(key(v, 'Enter')).toBe(true);
    const vp = v.view.get(ProteinSection)!;
    expect([vp.start, vp.end]).toEqual([3, 4]);
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

  it('a cursor on a follower moved out of view moves no section, though the follower leads a pan', async () => {
    const v = await make();
    v.view.setRange(ProteinSection, 0, 4);
    const gene = v.view.get(GeneSection)!;
    v.selection.setCursor({ section: GeneSection, system: G, position: Math.floor(gene.end) - 1 });
    key(v, 'ArrowRight');
    key(v, 'ArrowRight');
    expect([v.view.get(ProteinSection)!.start, v.view.get(ProteinSection)!.end]).toEqual([0, 4]);
    expect(v.view.pan(GeneSection, 3)).toBe(true);
    expect(v.view.get(ProteinSection)!.start).toBeCloseTo(1, 9);
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

  it('a shortcut with modifiers is consumed even when its commands refuse (Alt+ArrowLeft on a follower is no Back)', async () => {
    const v = await make();
    v.input.move(...at(v, 1, 0, 20));
    const start = v.view.get(GeneSection)!.start;
    expect(v.input.key({ key: 'ArrowLeft', altKey: true })).toBe(true);
    expect(v.view.get(GeneSection)!.start).toBe(start);
    // a key no command has stays the page's
    expect(v.input.key({ key: 'q', altKey: true })).toBe(false);
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
      return y >= s.y - 4 && y < s.y + s.height;
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

describe('on the Canvas2D backend', () => {
  async function mounted(config: readonly (readonly [any, unknown])[] = []) {
    const v = await make(spec(), config);
    const base = new FakeCanvas();
    const overlay = new FakeCanvas();
    const labels = new FakeCanvas();
    v.useBackend(new CanvasBackend(base, overlay, labels));
    frames!.flush();
    for (const c of [base, overlay, labels]) c.context.take();
    return { v, base, overlay };
  }

  it('a selection redraws the overlay alone: its tints and outlines, on both sections', async () => {
    const { v, base, overlay } = await mounted();
    const layouts = v.counts.layouts;
    click(v, at(v, 0, 0, 2));
    frames!.flush();
    expect(base.context.take()).toEqual([]);
    // a band on the protein's two tracks (a squircle, a path), one on the gene: a tint and an outline each
    expect(overlay.context.calls('fill').length + overlay.context.calls('fillRect').length).toBeGreaterThanOrEqual(2);
    expect(overlay.context.calls('stroke').length + overlay.context.calls('strokeRect').length).toBeGreaterThanOrEqual(
      2,
    );
    expect(v.counts.layouts).toBe(layouts);
  });

  it('with the tint under the data, a selection redraws the base too', async () => {
    const { SeqCoreConfig } = await import('@molstar/seqx-core/context');
    const { v, base } = await mounted([[SeqCoreConfig.render, { emphasis: { selectionPlace: 'under' } }]]);
    click(v, at(v, 0, 0, 2));
    frames!.flush();
    expect(base.context.calls('fillRect').length).toBeGreaterThan(0);
  });
});

describe('across a cyclic origin', () => {
  /** A cyclic system of 10, its section showing [7, 10) then [0, 3); an item W wrapping the origin, 9 to 1. */
  function cyclic(): Node {
    const b = createSeqSpec({ title: 'cyclic' });
    const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
    const items = b.data.ranges('items', {
      target: c,
      rows: { start: new Int32Array([9]), end: new Int32Array([1]), id: ['W'] },
    });
    b.views.view('main', (v) => {
      const s = v.section('c', { of: c, range: [7, 3] });
      s.track('letters').layer('l', { type: 'letters' });
      s.track('items').layer('d', { type: 'ranges', data: items });
    });
    return b.build();
  }
  const C = '/data/c';

  it('the axis wraps: [7, 10) then [0, 3)', async () => {
    const v = await make(cyclic());
    expect([...v.view.resolved.value!.sections[0].axis.ranges!]).toEqual([7, 10, 0, 3]);
  });

  it('a drag from 9 to 1 selects [9, 10) and [0, 2)', async () => {
    const v = await make(cyclic());
    const [x0, y] = at(v, 0, 0, 9);
    v.input.down(x0, y);
    v.input.move(at(v, 0, 0, 1)[0], y);
    v.input.up();
    expect(ranges(v, C)).toEqual([0, 2, 9, 10]);
  });

  it('Shift+arrows extend over the origin', async () => {
    const v = await make(cyclic());
    v.selection.setCursor({ section: '/views/main/c', system: C, position: 9 });
    v.input.key({ key: 'ArrowRight', shiftKey: true });
    v.input.key({ key: 'ArrowRight', shiftKey: true });
    expect(v.selection.cursor.value?.position).toBe(1);
    expect(ranges(v, C)).toEqual([0, 2, 9, 10]);
  });

  it('an item over the origin is selected whole, and an extend from it goes on in display order', async () => {
    const v = await make(cyclic());
    const [x, y0] = at(v, 0, 1, 9);
    const y = [0, 1, 2].map((l) => y0 + l * 12).find((ly) => v.hit(x, ly)?.element?.key === 'W')!;
    click(v, [x, y]);
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['W']);
    expect(ranges(v, C)).toEqual([0, 1, 9, 10]);
    click(v, at(v, 0, 0, 2), { shiftKey: true });
    expect(ranges(v, C)).toEqual([0, 3, 9, 10]);
  });
});

describe('a gesture’s selection follows the spec', () => {
  const rename = (at: string[], id: string) => ({
    ops: [{ op: 'move' as const, at, parent: at.slice(0, -1), place: { last: true as const }, rename: id }],
  });

  it('a cell and its anchor follow a rename of their system; an extend after it works', async () => {
    const v = await make();
    click(v, at(v, 0, 0, 2));
    const o = await v.spec.patch('Rename', rename(['data', 'p'], 'q'));
    expect(o.ok).toBe(true);
    expect(ranges(v, '/data/q')).toEqual([2, 3]);
    expect(v.selection.anchor?.system).toBe('/data/q');
    expect(v.selection.cursor.value?.system).toBe('/data/q');
    frames!.flush();
    click(v, at(v, 0, 0, 5), { shiftKey: true });
    expect(ranges(v, '/data/q')).toEqual([2, 6]);
  });

  it('an item follows a rename of its owner, and is dropped with an issue when its node is removed', async () => {
    const v = await make();
    click(v, at(v, 0, 1, 6));
    await v.spec.patch('Rename', rename(['data', 'domains'], 'd2'));
    expect(v.selection.value.value.elements.keys).toEqual([{ owner: '/data/d2', key: 'B' }]);
    expect(ranges(v)).toEqual([5, 9]);
    const o = await v.spec.patch('Remove', {
      ops: [
        { op: 'remove', at: ['views', 'main', 'protein', 'domains'] },
        { op: 'remove', at: ['data', 'd2'] },
      ],
    });
    expect(o.ok, JSON.stringify(o)).toBe(true);
    expect(v.selection.value.value.elements.keys).toEqual([]);
    expect(ranges(v)).toEqual([]);
    expect(v.addresses.issues.value.map((i) => i.code)).toEqual(['seq.address.missing']);
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
