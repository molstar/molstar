/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import {
  covariationGrid,
  ehzPairs,
  MoleculeFile,
  SeedFile,
  ScaleSize,
  scaleGrid,
  scaleLetters,
  seedRowPairs,
  SystemName,
  TrnaLetters,
} from './data.js';
import { DerivedData } from './derived-data.js';

// The tRNA views (examples.md, tRNA): PDB 1EHZ's chain A with its modified residues and 27 base pairs, the four stems
// nested and the tertiary pairs crossing them; the pairs as arcs, split by crossing, as chords on a ring and as a dot
// plot (the family's covariation, computed on 1EHZ's positions, above the diagonal, the structure's pairs below); the
// seed row 1EHZ is placed on with its `SS_cons` pairs; the seed's rows with their logo; and the dot plot at 4,000 ×
// 4,000. A modified residue is drawn as its parent letter, its code in a track of ranges (the mark is later).

export const title = 'tRNA-Phe';

type Builder = ReturnType<typeof createSeqSpec>;

/** 1EHZ: its letters, its modified residues from the committed file, its base pairs. */
function molecule(b: Builder) {
  const file = b.data.source('ehz', { format: 'seq', uri: MoleculeFile });
  // the letters inline (compared with the file's by the fixtures test): pairs on a sequence whose letters come from a
  // source are not placed (the section's units are cached before the source opens)
  const trna = b.data.sequence('trna', { alphabet: 'rna', letters: TrnaLetters, name: SystemName });
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
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
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
        color: theme.field('value', { domain: [0, 1.2], palette: theme.continuous({ colors: 'viridis' }) }),
      });
      t.layer('structure', { type: 'pairs', data: pairs, mirror: true, triangle: 'lower' });
    });
  });
  return b.build();
}

/** Crossing pairs on the other side: the largest nested set above, the tertiary pairs that cross it below. */
export const split = () => build({ crossing: 'split' });

/** The tRNA on a ring: its pairs as chords through the hole. */
export const ring = () => build({ layout: 'rings' });

/** The arcs over the letters, and the base-pair dot plot under them following their window. */
export const dotPlot = () => build({ dotPlot: true });

/** The arcs' data with the pairs of a kind left out: what a patch of the layer's `data` writes (the editing view). */
export const withoutKind = (pairs: string, kind: string) => `${pairs}/(not :of (item :type ${kind}))`;

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
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    s.track('pairs', { label: 'SS_cons pairs' }).layer('arcs', { type: 'pairs', data: pairs });
    s.track('letters', { label: DerivedData.row }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The Rfam RF00005 seed (a Stockholm file as published): its rows over its columns under a conservation logo. */
export function seed(): Node {
  const b = createSeqSpec({ title: 'Rfam RF00005 seed' });
  const file = b.data.source('seed-file', { format: 'stockholm', uri: SeedFile });
  const aln = b.data.alignment('seed', { alphabet: 'rna', source: file, reference: DerivedData.row });
  b.views.view('main', (v) => {
    const s = v.section('family', { of: aln });
    s.track('ruler', { label: 'Columns' }).layer('ruler', { type: 'ruler' });
    s.track('logo', { label: 'Conservation' }).layer('logo', { type: 'logo', data: aln });
    s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
  });
  return b.build();
}

/** The dot plot at `size` × `size` (4,000 by default): generated letters and values from seeds, the plot alone. */
export function scale(size: number = ScaleSize): Node {
  const b = createSeqSpec({
    title: `tRNA dot plot, ${size.toLocaleString('en-US')} × ${size.toLocaleString('en-US')}`,
  });
  const seq = b.data.sequence('rna', { alphabet: 'rna', letters: scaleLetters(size) });
  const grid = b.data.grid('probability', { x: seq, values: scaleGrid(size) });
  b.views.view('main', (v) => {
    const s = v.section('rna', { of: seq });
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    v.section('dots', { layout: 'grid', of: seq, follow: s.address }).track('plot', { label: 'Probability' }, (t) =>
      t.layer('probability', {
        type: 'grid',
        data: grid,
        // zoomed out, the strongest cell of a sample: the helices stay in view
        aggregate: 'max',
        color: theme.field('value', { domain: [0, 1], palette: theme.continuous({ colors: 'viridis' }) }),
      }),
    );
  });
  return b.build();
}
