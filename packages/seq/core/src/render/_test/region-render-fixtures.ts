/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Region reads' render fixtures (_spec/README.md "Render", The grid layer: by region; The matrix layer): a symmetric
// contact map at two levels, links between hotspots at three levels, genotypes with a long deletion over items and bins, generated, in the `seq` format.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { GridDense, MatrixBins, MatrixItems, PairsList, SeqData } from '@molstar/seqx-core/data';
import { binPairs, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** A contact map's chromosome length. */
export const HicLength = 2048;

/** A contact map over chr (2,048 positions): a band along the diagonal, served symmetric at levels 1 and 4. */
export function hic(scale = 1): SeqData {
  const dense = (level: number): GridDense => {
    const n = HicLength / level;
    const v = new Float32Array(n * n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) v[i * n + j] = (scale * level * level * 100) / (1 + Math.abs(i - j));
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

/**
 * A top track section and a grid section over chr, its grid read by region; with `second`, a second source's grid on
 * the same pair (EC-59: the upper triangle one map, the lower the other).
 */
export function gridRegionSpec(second = false) {
  const b = createSeqSpec({ title: 'Hi-C by region' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(hic()) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const g = b.data.grid('contacts', { x: chr, source: src } as never);
  const other = second ? b.data.source('other', { format: 'seq', value: seqToValue(hic(3)) as never }) : undefined;
  const h = other ? b.data.grid('treated', { x: chr, source: other } as never) : undefined;
  b.views.view('main', (v) => {
    const top = v.section('chr', { of: chr });
    top.track('ruler').layer('ruler', { type: 'ruler' });
    v.section('plot', { layout: 'grid', of: chr, follow: top.address }).track('plot', {}, (t) => {
      t.layer('grid', { type: 'grid', data: g, ...(h ? { triangle: 'upper' } : {}) });
      if (h) t.layer('treated', { type: 'grid', data: h, triangle: 'lower' });
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

/** A links chromosome's length. */
export const LinksLength = 1_000_000;
const Length = LinksLength;
const LinksLevels = [1, 1024, 16_384];

/** `n` links between 20 hotspots of a chromosome (each end jittered 2 kb), a level of bins at each coarse level. */
export function linksSource(n: number): { data: SeqData; list: PairsList } {
  const rng = new PCG(9);
  const spots = Array.from({ length: 20 }, (_, i) => 20_000 + i * 48_000);
  const as = new Int32Array(n);
  const bs = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    const a = spots[Math.floor(rng.float() * spots.length)];
    const b = spots[Math.floor(rng.float() * spots.length)];
    as[i] = a + Math.floor(rng.float() * 2000);
    bs[i] = b + Math.floor(rng.float() * 2000);
  }
  const list: PairsList = {
    form: 'list',
    aStart: as,
    aEnd: as.map((s) => s + 10),
    bStart: bs,
    bEnd: bs.map((s) => s + 10),
    id: categoryColumn(Array.from({ length: n }, (_, i) => `l${i}`)),
  };
  const whole = (level: number) => ({ system: 'chr', start: 0, end: Length, level, shapes: ['pairs' as const] });
  return {
    list,
    data: {
      coverage: {
        systems: [{ name: 'chr', length: Length }],
        shapes: [{ system: 'chr', shape: 'pairs' }],
        levels: LinksLevels,
      },
      pairs: LinksLevels.map((level) => ({
        request: whole(level),
        answer: { ok: true as const, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
      })),
    },
  };
}

/** A pairs layer over chr, its links read by region. */
export function linksSpec(n: number) {
  const b = createSeqSpec({ title: 'links by region' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(linksSource(n).data) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const links = b.data.pairs('links', { a: chr, source: src } as never);
  b.views.view('main', (v) => {
    const s = v.section('chr', { of: chr });
    s.track('links').layer('links', { type: 'pairs', data: links });
  });
  return b.build();
}
