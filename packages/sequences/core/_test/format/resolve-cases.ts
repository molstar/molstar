/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Resolve's per-row cases, for their heap guard and its quick twin (resolve.heap.ts, resolve.quick.ts): a step is one row.
// A 300-residue protein over its 900-nucleotide gene, and `ROWS` one-residue rows on it with categories and values:
// `cover` (what `positions` and an annotation focus mark) through a codon map, the filter's test of a row's value (a
// number column compared as numbers, no double boxed), and a row's color through the layer's one categorical palette.
// A step's index wraps at `ROWS`, so a longer run than the guard's keeps measuring.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqSet } from '@molstar/seq-core/model';
import { expect } from 'vitest';
import { compile, type End, layerPalette, paletteValue, resolveView, supportedAll } from '../../src/format/index.js';
import { cover, matcher } from '../../src/format/resolve/reach.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

/** The rows: four times the guard's steps in its shorter run, the index its longer run reaches. */
export const ROWS = 80_000;

export interface ResolveCases {
  /** `cover`: a row's codons marked through the map. */
  cover: (i: number) => void;
  /** The mask `cover` marks, one byte a nucleotide. */
  mask: Uint8Array;
  path: { span: number };
  /** The filter's test of a row, a string and a number field. */
  filter: (i: number) => void;
  kept: () => number;
  /** A row's color through the categorical palette. */
  color: (i: number) => void;
  colored: () => number;
}

let loaded: Promise<ResolveCases> | undefined;

/** The set and view built once for every case. */
export function resolveCases(): Promise<ResolveCases> {
  loaded ??= build();
  return loaded;
}

async function build(): Promise<ResolveCases> {
  const reg = registries();
  const b = createSeqViewSpec();
  const protein = Array.from({ length: 300 }, (_, i) => 'ACDEFGHIKLMNPQRSTVWY'[i % 20]).join('');
  const gene = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: 'ATG'.repeat(300) });
  b.data.sequence({ ref: 'p', alphabet: 'protein', value: protein });
  b.data.mapping({ ref: 'm', from: 'p', to: 'gene', method: 'codon', cds: [[0, 900]] } as never);
  const start = Array.from({ length: ROWS }, (_, i) => i % 300);
  const category = start.map((i) => ['helix', 'strand', 'turn'][i % 3]);
  const value = start.map((i) => i / 300);
  b.data.annotation({ ref: 'a', geometry: 'ranges', target: 'p', value: { start, category, value } });
  b.view({ ref: 'v', label: 'V' }).block({ sequence: gene }).track({ label: 'A' }).ranges({ data: 'a' });
  const doc = b.getDocument();
  const r = await compile(doc, reg, fixtureLoaders());
  expect(r.issues).toEqual([]);
  const set: SeqSet = r.set;
  const view = resolveView(doc, set, 'v', reg, supportedAll(reg)).view;
  const a = set.annotations[view.data[0]] as NonNullable<SeqSet['annotations'][number]>;

  const path = view.paths[view.placements.path[0]];
  const ends: End[] = Array.from({ length: ROWS }, (_, i) => ({ system: path.from, slot: 0, side: 0, rows: [i] }));
  const mask = new Uint8Array(900);

  const keep = matcher(['helix', 0.5, '1']);
  let kept = 0;
  // A continuous palette's color of a row is the viewer's per row (LY3 reads the number column): a double passed to
  // the palette's function is boxed, so it is not guarded here.
  const f = layerPalette(set, [a], { kind: 'categorical' }, 'category', [0, 1]);
  let colored = 0;

  return {
    cover: (i) => cover(a, ends[i % ROWS], path, mask),
    mask,
    path,
    filter: (i) => {
      const k = i % ROWS;
      kept += +keep(a, set, k, 'category') + +keep(a, set, k, 'value');
    },
    kept: () => kept,
    color: (i) => {
      colored += f(paletteValue(a, set, i % ROWS, 'category')) === undefined ? 0 : 1;
    },
    colored: () => colored,
  };
}
