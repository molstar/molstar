/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import {
  addRefPatch,
  createSeqSpec,
  modified as modifiedTable,
  type RefInput,
  SeqLanguage,
  theme,
} from '@molstar/seqx-core/lang';
import {
  covariationGrid,
  ehzPairs,
  pairRows,
  ModifiedCompounds,
  ModifiedResidues,
  MoleculeFile,
  SeedFile,
  ScaleSize,
  scaleGrid,
  scaleLetters,
  seedRowPairs,
  SystemName,
} from './data.js';
import { DerivedData } from './derived-data.js';

// The tRNA views: PDB 1EHZ's chain A with its modified residues and 27 base pairs, the four stems
// nested and the six tertiary pairs among them (four cross the stems' pairs); the pairs as arcs, split by crossing, as chords on a ring and as a dot
// plot (the family's covariation, computed on 1EHZ's positions, above the diagonal, the structure's pairs below); the
// seed row 1EHZ is placed on with its `SS_cons` pairs; the seed's rows with their logo; and the dot plot at 4,000 ×
// 4,000. A modified residue is its parent's letter with its component in the sequence's `modified` table (a dot at the
// cell's top right; its code and name on hover), its code also in a track of ranges, a label where it fits; the
// modified residues found by query (`(elements :modified true)`), a bookmark; the arcs' legend lists the pair kinds.

export const title = 'tRNA-Phe';

type Builder = ReturnType<typeof createSeqSpec>;

/** The bookmark the builder writes: the modified residues, found by query. */
export const ModifiedRef: RefInput = {
  id: 'modified',
  to: '/data/trna/(elements :modified true)',
  label: 'Modified residues',
};

/** 1EHZ: its letters and modified residues from the committed file, its components, its base pairs. */
function molecule(b: Builder) {
  const file = b.data.source('ehz', { format: 'seq', uri: MoleculeFile });
  for (const [code, c] of Object.entries(ModifiedCompounds))
    b.data.compound(code, { name: c.name, parents: [c.parent] });
  // the letters from the file: a grid over them is read once the file gives their length
  const trna = b.data.sequence('trna', {
    alphabet: 'rna',
    source: file,
    name: SystemName,
    modified: modifiedTable(ModifiedResidues.map(([position, component]) => ({ position, component }))),
  });
  const modified = b.data.ranges('modified', { target: trna, source: file });
  const pairs = b.data.pairs('pairs', { a: trna, rows: ehzPairs() });
  return { trna, modified, pairs };
}

export interface TrnaOptions {
  readonly crossing?: 'same' | 'split';
  readonly layout?: 'rings';
  readonly dotPlot?: boolean;
}

/** 1EHZ's letters, its modified residues, its pairs, and a view of them. */
export function build(options: TrnaOptions = {}): Node {
  const b = createSeqSpec({ title });
  const { trna, modified, pairs } = molecule(b);
  const covariation = options.dotPlot ? b.data.grid('covariation', { x: trna, values: covariationGrid() }) : undefined;
  b.views.view('main', (v) => {
    const s = v.section('trna', { of: trna, ...(options.layout ? { layout: options.layout } : {}) });
    const ring = options.layout === 'rings';
    const letters = () => s.track('letters', { label: '1EHZ' }).layer('letters', { type: 'letters' });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    if (!ring)
      s.track('modified', { label: 'Modified' }).layer('modified', {
        type: 'ranges',
        data: modified,
        color: theme.field('category', {
          palette: theme.categorical({ colors: 'categories', values: ['A', 'C', 'G', 'U'] }),
        }),
      });
    // the arcs over the letters, their feet on them (on a ring, inside: chords through the hole)
    if (ring) letters();
    s.track('pairs', { label: 'Base pairs' }).layer('arcs', {
      type: 'pairs',
      data: pairs,
      ...(options.crossing ? { crossing: options.crossing } : {}),
      ...(ring ? { side: 'below' } : {}),
    });
    if (!ring) letters();
    if (!covariation) return;
    // the dot plot: x and y 1EHZ's, y following x; the covariation above the diagonal, the pairs below it (mirrored:
    // the builder sets `mirror`, a plain `false` otherwise)
    v.section('dots', { layout: 'grid', of: trna, follow: s.address }).track('plot', { label: 'Dot plot' }, (t) => {
      t.layer('covariation', {
        type: 'grid',
        data: covariation,
        triangle: 'upper',
        color: theme.field('value', { domain: [0, 1.2], palette: theme.continuous({ colors: 'ramp' }) }),
      });
      t.layer('structure', { type: 'pairs', data: pairs, mirror: true, triangle: 'lower' });
    });
  });
  const root = b.build();
  const r = applyPatch(root, addRefPatch(root, ModifiedRef), { language: SeqLanguage });
  if (!r.ok) throw new Error(`The bookmark did not apply: ${r.issues.map((i) => i.code).join(', ')}`);
  return r.root;
}

/** Crossing pairs on the other side: the largest nested set above, the tertiary pairs that cross it below. */
export const split = () => build({ crossing: 'split' });

/** The tRNA on a ring: its pairs as chords through the hole. */
export const ring = () => build({ layout: 'rings' });

/** The arcs over the letters, and the base-pair dot plot under them following their window. */
export const dotPlot = () => build({ dotPlot: true });

/** The arcs' data with the pairs a field has a value on left out: what a patch of the layer's `data` writes (the editing view). */
export const without = (pairs: string, field: 'type' | 'tertiary', value: string) =>
  `${pairs}/(except :of (item :${field} ${value}))`;

/**
 * The seed row 1EHZ is placed on (the row of highest identity, computed: `DerivedData.row`): its letters and the pairs
 * of the seed's `SS_cons` projected onto it, nested only, as arcs.
 */
export function seedRow(): Node {
  const b = createSeqSpec({ title: `Rfam RF00005 seed: ${DerivedData.row}` });
  const row = b.data.sequence('row', { alphabet: 'rna', letters: DerivedData.rowLetters, name: DerivedData.row });
  const pairs = b.data.pairs('pairs', { a: row, rows: seedRowPairs() });
  b.views.view('main', (v) => {
    const s = v.section('row', { of: row });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    s.track('pairs', { label: 'SS_cons pairs' }).layer('arcs', { type: 'pairs', data: pairs });
    s.track('letters', { label: DerivedData.row }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/**
 * 1EHZ's positions on the seed's columns: the runs of its computed placement (`DerivedData.columns`, through the seed
 * row of highest identity), each a `segments` map's segment (the stopgap of an `alignment` map: a chain onto a seed
 * row's columns).
 */
export function ehzOnColumns() {
  const from: number[] = [];
  const to: number[] = [];
  const length: number[] = [];
  DerivedData.columns.forEach((c, i) => {
    if (c < 0) return;
    const k = length.length - 1;
    if (k >= 0 && from[k] + length[k] === i && to[k] + length[k] === c) length[k]++;
    else {
      from.push(i);
      to.push(c);
      length.push(1);
    }
  });
  return { from, to, length };
}

/**
 * The Rfam RF00005 seed (a Stockholm file as published): its `SS_cons` pairs as arcs over its columns, its rows under a
 * conservation logo; with `chain`, 1EHZ's letters below, mapped onto the columns (`ehzOnColumns`), so an arc's hover
 * reaches its residues and a residue's its column.
 */
export function seed(options: { readonly chain?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'Rfam RF00005 seed' });
  const file = b.data.source('seed-file', { format: 'stockholm', uri: SeedFile });
  const aln = b.data.alignment('seed', { alphabet: 'rna', source: file, reference: DerivedData.row });
  const ss = b.data.pairs('ss-cons', {
    a: aln,
    rows: pairRows(DerivedData.seedPairs.map(([i, j]) => [i, j, 'canonical'])),
  });
  const ehz = options.chain ? molecule(b).trna : undefined;
  if (ehz) b.data.map('onto-seed', { from: ehz, to: aln, method: 'segments', segments: ehzOnColumns() });
  b.views.view('main', (v) => {
    const s = v.section('family', { of: aln });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    s.track('ss-cons', { label: 'SS_cons' }).layer('arcs', { type: 'pairs', data: ss });
    s.track('logo', { label: 'Conservation' }).layer('logo', { type: 'logo', data: aln });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
    if (ehz) v.section('ehz', { of: ehz }).track('letters', { label: '1EHZ' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The seed with 1EHZ below it, on its columns. */
export const seedWithChain = () => seed({ chain: true });

/** The dot plot at `size` × `size` (4,000 by default): generated letters and values from seeds, the plot alone. */
export function scale(size: number = ScaleSize): Node {
  const b = createSeqSpec({
    title: `tRNA dot plot, ${size.toLocaleString('en-US')} × ${size.toLocaleString('en-US')}`,
  });
  const seq = b.data.sequence('rna', { alphabet: 'rna', letters: scaleLetters(size) });
  const grid = b.data.grid('probability', { x: seq, values: scaleGrid(size) });
  b.views.view('main', (v) => {
    const s = v.section('rna', { of: seq });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    v.section('dots', { layout: 'grid', of: seq, follow: s.address }).track('plot', { label: 'Probability' }, (t) =>
      t.layer('probability', {
        type: 'grid',
        data: grid,
        // zoomed out, the strongest cell of a sample: the helices stay in view
        aggregate: 'max',
        color: theme.field('value', { domain: [0, 1], palette: theme.continuous({ colors: 'ramp' }) }),
      }),
    );
  });
  return b.build();
}
