/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Region reads' render fixtures (_spec/README.md "Render", The grid layer: by region; The matrix layer): a symmetric
// contact map at two levels, genotypes with a long deletion over items and bins, generated, in the `seq` format.

import type { GridDense, MatrixBins, MatrixItems, SeqData } from '@molstar/seqx-core/data';
import { seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** A contact map's chromosome length. */
export const HicLength = 2048;

/** A contact map over chr (2,048 positions): a band along the diagonal, served symmetric at levels 1 and 4. */
export function hic(): SeqData {
  const dense = (level: number): GridDense => {
    const n = HicLength / level;
    const v = new Float32Array(n * n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) v[i * n + j] = (level * level * 100) / (1 + Math.abs(i - j));
    return { form: 'dense', x: 0, y: 0, nx: n, ny: n, value: v };
  };
  const request = (level: number) => ({
    system: 'chr',
    start: 0,
    end: HicLength,
    level,
    shapes: ['grid' as const],
    y: { system: 'chr', start: 0, end: HicLength },
  });
  return {
    coverage: {
      systems: [{ name: 'chr', length: HicLength }],
      shapes: [{ system: 'chr', shape: 'grid', symmetric: true, aggregate: 'sum' }],
      levels: [1, 4],
    },
    pairs: [1, 4].map((level) => ({
      request: request(level),
      answer: { ok: true as const, tables: { grid: dense(level) } },
    })),
  };
}

/** A top track section and a grid section over chr, its grid read by region. */
export function gridRegionSpec() {
  const b = createSeqSpec({ title: 'Hi-C by region' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(hic()) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const g = b.data.grid('contacts', { x: chr, source: src } as never);
  b.views.view('main', (v) => {
    const top = v.section('chr', { of: chr });
    top.track('ruler').layer('ruler', { type: 'ruler' });
    v.section('plot', { layout: 'grid', of: chr, follow: top.address }).track('plot', {}, (t) => {
      t.layer('grid', { type: 'grid', data: g });
    });
  });
  return b.build();
}

/** A genotype matrix's chromosome length. */
export const GenotypesLength = 262_144;
const N = GenotypesLength;
export const States = ['NA12878', 'NA12891', 'NA12892', 'HG00096'];
/** A ring's 40 states (a heatmap ring of samples, Circos). */
export const Many = Array.from({ length: 40 }, (_, i) => `s${i}`);

/** Genotypes over chr: a SNP every 100 positions in 0..2,000, a deletion over 500..1,500 in the first state. */
export function genotypes(names: readonly string[] = States): SeqData {
  const S = names.length;
  const start: number[] = [];
  const end: number[] = [];
  for (let p = 0; p < 2000; p += 100) start.push(p), end.push(p + 1);
  start.push(500), end.push(1500);
  const n = start.length;
  const codes = new Uint8Array(n * S);
  for (let i = 0; i < n; i++) for (let s = 0; s < S; s++) codes[i * S + s] = (i + s) % 3;
  // the deletion: every copy in the first state, none in the others; one state uncalled at the first SNP
  codes[(n - 1) * S] = 2;
  codes[1] = 255;
  const items: MatrixItems = {
    form: 'items',
    start: Int32Array.from(start),
    end: Int32Array.from(end),
    states: S,
    codes,
  };
  const binsN = N / 256;
  const classes = new Uint8Array(binsN * S);
  const shares = new Uint8Array(binsN * S).fill(200);
  for (let b = 0; b < binsN; b++) for (let s = 0; s < S; s++) classes[b * S + s] = (b + s) % 3;
  const bins: MatrixBins = { form: 'bins', start: 0, states: S, classes, shares };
  const request = (level: number) => ({ system: 'chr', start: 0, end: N, level, shapes: ['matrix' as const] });
  return {
    coverage: {
      systems: [{ name: 'chr', length: N }],
      shapes: [{ system: 'chr', shape: 'matrix', states: [...names] }],
      levels: [1, 256],
    },
    pairs: [
      { request: request(1), answer: { ok: true as const, tables: { matrix: items } } },
      { request: request(256), answer: { ok: true as const, tables: { matrix: bins } } },
    ],
  };
}

/** A genotype matrix read by region, in a tracks or a rings section. */
export function sourced(names: readonly string[] = States, layout: 'tracks' | 'rings' = 'tracks') {
  const b = createSeqSpec({ title: 'Genotypes' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(genotypes(names)) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const m = b.data.matrix('gt', { target: chr, source: src } as never);
  b.views.view('main', (v) => {
    v.section('chr', { of: chr, layout }).track('gt', { height: 80 }, (t) => {
      t.layer('gt', { type: 'matrix', data: m });
    });
  });
  return b.build();
}
