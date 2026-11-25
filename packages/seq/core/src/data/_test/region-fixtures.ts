/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type {
  Coverage,
  GridDense,
  GridSparse,
  MatrixItems,
  PairsList,
  RangesList,
  RegionRequest,
  RowsKeyed,
  SystemInfo,
} from '../protocol.js';
import type { SeqData, SeqPair } from '../seq-format.js';
import { binPairs } from '../summaries-pairs.js';
import { gridPlane, halveGrid } from '../summaries-grid.js';

// Small generated region data (the region reads' tests): pairs over two systems found from either end, sided pairs
// between two genomes, a symmetric contact grid with two levels and its one-level twin, a matrix of codes with a long
// deletion, values with missing ones, keyed rows, transcripts as pieces with item extents. Every source answers its
// system whole at each level it lists, so the runtime cuts any region from it.

export const Chr1 = 5000;
export const Chr2 = 3000;
export const Systems: SystemInfo[] = [
  { name: 'chr1', length: Chr1 },
  { name: 'chr2', length: Chr2 },
];

/** `n` pairs within and between chr1 and chr2: strands, and (unless `light`) types, ids and a score field. */
export function pairsOf(n: number, seed = 7, light = false): PairsList {
  const rng = new PCG(seed);
  const aSys: string[] = [];
  const bSys: string[] = [];
  const as = new Int32Array(n);
  const ae = new Int32Array(n);
  const bs = new Int32Array(n);
  const be = new Int32Array(n);
  const aStrand = new Int8Array(n);
  const bStrand = new Int8Array(n);
  const type: string[] = [];
  const score = new Float32Array(n);
  const id: string[] = [];
  for (let i = 0; i < n; i++) {
    const a = rng.float() < 0.7 ? 'chr1' : 'chr2';
    const b = rng.float() < 0.5 ? 'chr1' : 'chr2';
    const al = a === 'chr1' ? Chr1 : Chr2;
    const bl = b === 'chr1' ? Chr1 : Chr2;
    aSys.push(a);
    bSys.push(b);
    as[i] = Math.floor(rng.float() * (al - 40));
    ae[i] = as[i] + 1 + Math.floor(rng.float() * 30);
    bs[i] = Math.floor(rng.float() * (bl - 40));
    be[i] = bs[i] + 1 + Math.floor(rng.float() * 30);
    aStrand[i] = rng.float() < 0.5 ? 1 : -1;
    bStrand[i] = rng.float() < 0.5 ? 1 : -1;
    type.push(rng.float() < 0.8 ? 'contact' : 'other');
    score[i] = Math.floor(rng.float() * 100);
    id.push(`p${i}`);
  }
  return {
    form: 'list',
    aSystem: categoryColumn(aSys),
    bSystem: categoryColumn(bSys),
    aStart: as,
    aEnd: ae,
    bStart: bs,
    bEnd: be,
    aStrand,
    bStrand,
    ...(!light && { type: categoryColumn(type) }),
    ...(!light && { id: categoryColumn(id), fields: [{ name: 'score', column: score }] }),
  };
}

/** The rows of a list with an end on `system` (a source's whole read of it). */
function withEndOn(t: PairsList, system: string): PairsList {
  const rows: number[] = [];
  for (let i = 0; i < t.aStart.length; i++)
    if (t.aSystem!.get(i) === system || t.bSystem!.get(i) === system) rows.push(i);
  const g = <A extends Int32Array | Int8Array | Float32Array>(c: A) =>
    (c.constructor as unknown as { from(a: number[]): A }).from(rows.map((i) => c[i]));
  const cat = (c: ReturnType<typeof categoryColumn>) => categoryColumn(rows.map((i) => c.get(i)));
  return {
    form: 'list',
    aSystem: cat(t.aSystem!),
    bSystem: cat(t.bSystem!),
    aStart: g(t.aStart),
    aEnd: g(t.aEnd),
    bStart: g(t.bStart),
    bEnd: g(t.bEnd),
    aStrand: g(t.aStrand!),
    bStrand: g(t.bStrand!),
    ...(t.type && { type: cat(t.type) }),
    ...(t.id && { id: cat(t.id), fields: [{ name: 'score', column: g(t.fields![0].column as Float32Array) }] }),
  };
}

const whole = (system: string, length: number, level: number, shapes: RegionRequest['shapes']): RegionRequest => ({
  system,
  start: 0,
  end: length,
  level,
  shapes,
});

/** Pairs served at level 1 only (the store's summaries give them levels): each system's pairs whole. */
export function pairsSource(n = 400, light = false): SeqData {
  const all = pairsOf(n, 7, light);
  const coverage: Coverage = {
    systems: Systems,
    shapes: Systems.map((s) => ({
      system: s.name,
      shape: 'pairs' as const,
      binBy: 'type',
      fields: [{ name: 'score', type: 'float32' as const }],
    })),
    levels: [1],
  };
  const pairs: SeqPair[] = Systems.map((s) => ({
    request: whole(s.name, s.length!, 1, ['pairs']),
    answer: { ok: true, tables: { pairs: withEndOn(all, s.name) } },
  }));
  return { coverage, pairs };
}

/** Pairs with their own levels: items at 1, bins at 16 and 256 (each system's bins with an end on it). */
export function leveledPairsSource(n = 400): SeqData {
  const all = pairsOf(n);
  const levels = [1, 16, 256];
  const coverage: Coverage = {
    systems: Systems,
    shapes: Systems.map((s) => ({ system: s.name, shape: 'pairs' as const, binBy: 'type' })),
    levels,
  };
  const pairs: SeqPair[] = [];
  for (const s of Systems)
    for (const level of levels) {
      const own = withEndOn(all, s.name);
      const table = level === 1 ? own : binPairs(own, level, { category: true });
      pairs.push({
        request: whole(s.name, s.length!, level, ['pairs']),
        answer: { ok: true, tables: { pairs: table } },
      });
    }
  return { coverage, pairs };
}

/** Sided pairs between two genomes (a synteny): `chr1` on both sides, read by side. */
export function sidedPairsSource(n = 120): SeqData {
  const rng = new PCG(11);
  const aStart = new Int32Array(n);
  const aEnd = new Int32Array(n);
  const bStart = new Int32Array(n);
  const bEnd = new Int32Array(n);
  const row = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    aStart[i] = Math.floor(rng.float() * (Chr1 - 100));
    aEnd[i] = aStart[i] + 50;
    bStart[i] = Math.floor(rng.float() * (Chr1 - 100));
    bEnd[i] = bStart[i] + 50;
    row[i] = Math.floor(i / 3);
  }
  const list: PairsList = {
    form: 'list',
    aSystem: categoryColumn(Array(n).fill('chr1')),
    bSystem: categoryColumn(Array(n).fill('chr1')),
    aStart,
    aEnd,
    bStart,
    bEnd,
    row,
  };
  const coverage: Coverage = {
    systems: [
      { name: 'chr1', length: Chr1, side: 'a' },
      { name: 'chr1', length: Chr1, side: 'b' },
    ],
    shapes: [{ system: 'chr1', shape: 'pairs' }],
    levels: [1],
  };
  const pairs: SeqPair[] = (['a', 'b'] as const).map((side) => ({
    request: { ...whole('chr1', Chr1, 1, ['pairs']), side },
    answer: { ok: true, tables: { pairs: list } },
  }));
  return { coverage, pairs };
}

/** A symmetric contact map over chr1 × chr1 and chr1 × chr2 (a band along the diagonal), counts. */
export function contactCells(nx: number, ny: number, diagonal: boolean, seed = 3): Float32Array {
  const rng = new PCG(seed);
  const v = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      const d = Math.abs(i - j);
      v[i * ny + j] = diagonal ? Math.floor(1000 / (1 + d)) : Math.floor(rng.float() * 4);
    }
  if (diagonal) for (let i = 0; i < nx; i++) for (let j = 0; j < i && j < ny; j++) v[i * ny + j] = v[j * ny + i];
  return v;
}

/** The grid shapes: chr1 with itself and chr2, symmetric, counts. */
function gridCoverage(levels: number[]): Coverage {
  return {
    systems: Systems,
    shapes: [{ system: 'chr1', shape: 'grid', y: ['chr1', 'chr2'], symmetric: true, aggregate: 'sum' }],
    levels,
  };
}

/** A symmetric contact grid with one level (the summaries give it more), at `scale` positions a cell. */
export function oneLevelGrid(length1 = 600, length2 = 400): SeqData {
  const sys: SystemInfo[] = [
    { name: 'chr1', length: length1 },
    { name: 'chr2', length: length2 },
  ];
  const self: GridDense = {
    form: 'dense',
    x: 0,
    y: 0,
    nx: length1,
    ny: length1,
    value: contactCells(length1, length1, true),
  };
  const cells = contactCells(length1, length2, false);
  const xBin: number[] = [];
  const yBin: number[] = [];
  const value: number[] = [];
  for (let i = 0; i < length1; i++)
    for (let j = 0; j < length2; j++)
      if (cells[i * length2 + j] > 0) {
        xBin.push(i);
        yBin.push(j);
        value.push(cells[i * length2 + j]);
      }
  const other: GridSparse = {
    form: 'sparse',
    x: 0,
    y: 0,
    xBin: Int32Array.from(xBin),
    yBin: Int32Array.from(yBin),
    value: Float32Array.from(value),
    fill: 0,
  };
  return {
    coverage: { ...gridCoverage([1]), systems: sys },
    pairs: [
      {
        request: { ...whole('chr1', length1, 1, ['grid']), y: { system: 'chr1', start: 0, end: length1 } },
        answer: { ok: true, tables: { grid: self } },
      },
      {
        request: { ...whole('chr1', length1, 1, ['grid']), y: { system: 'chr2', start: 0, end: length2 } },
        answer: { ok: true, tables: { grid: other } },
      },
    ],
  };
}

/** The same grid served with two levels of its own (1 and 4: a source with levels). */
export function twoLevelGrid(length1 = 600, length2 = 400): SeqData {
  const one = oneLevelGrid(length1, length2);
  const pairs: SeqPair[] = [...one.pairs];
  for (const p of one.pairs) {
    const t = (p.answer as { tables: { grid: GridDense | GridSparse } }).tables.grid;
    const ny = p.request.y!.end;
    let plane = gridPlane(t, p.request.end, ny);
    plane = halveGrid(halveGrid(plane, 'sum'), 'sum');
    const grid: GridDense = { form: 'dense', x: 0, y: 0, nx: plane.nx, ny: plane.ny, value: plane.value };
    pairs.push({ request: { ...p.request, level: 4 }, answer: { ok: true, tables: { grid } } });
  }
  return { coverage: { ...one.coverage, levels: [1, 4] }, pairs };
}

/** Genotype codes over chr1: SNPs and one long deletion across them, 4 states, item-major. */
export function genotypes(states = 4, n = 60): { matrix: MatrixItems; deletion: number } {
  const rng = new PCG(5);
  const start = new Int32Array(n + 1);
  const end = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) {
    start[i] = i * 80 + 10;
    end[i] = start[i] + 1;
  }
  // the long deletion: one item across the SNPs from 1,000 to 2,000 (EC-222)
  start[n] = 1000;
  end[n] = 2000;
  const order = Array.from({ length: n + 1 }, (_, i) => i).sort((a, b) => start[a] - start[b]);
  const s = Int32Array.from(order, (i) => start[i]);
  const e = Int32Array.from(order, (i) => end[i]);
  const codes = new Uint8Array((n + 1) * states);
  for (let k = 0; k < order.length; k++)
    for (let j = 0; j < states; j++)
      codes[k * states + j] =
        order[k] === n ? (j === 0 ? 2 : 0) : rng.float() < 0.1 ? 255 : Math.floor(rng.float() * 3);
  const id = categoryColumn(order.map((i) => (i === n ? 'del1' : `rs${i}`)));
  return { matrix: { form: 'items', start: s, end: e, states, codes, id }, deletion: order.indexOf(n) };
}

/** A matrix source over chr1 at level 1. */
export function matrixSource(states = 4): SeqData {
  const { matrix } = genotypes(states);
  return {
    coverage: {
      systems: [Systems[0]],
      shapes: [{ system: 'chr1', shape: 'matrix', states: Array.from({ length: states }, (_, i) => `S${i}`) }],
      levels: [1],
    },
    pairs: [{ request: whole('chr1', Chr1, 1, ['matrix']), answer: { ok: true, tables: { matrix } } }],
  };
}

/** Keyed rows on a records system of 50. */
export function rowsSource(): SeqData {
  const n = 50;
  const rows: RowsKeyed = {
    form: 'keyed',
    start: 0,
    key: categoryColumn(Array.from({ length: n }, (_, i) => `gene${i}`)),
    fields: [{ name: 'q', column: Float32Array.from({ length: n }, (_, i) => i / n) }],
  };
  return {
    coverage: {
      systems: [{ name: 'results', length: n }],
      shapes: [{ system: 'results', shape: 'rows' }],
      levels: [1],
    },
    pairs: [{ request: whole('results', n, 1, ['rows']), answer: { ok: true, tables: { rows } } }],
  };
}

/** Transcripts as exon pieces with their whole extent (BED12-like): one with a long intron. */
export function transcripts(): RangesList {
  // t1: [100, 200) and [3000, 3100), its extent [100, 3100); t2: [400, 500)
  return {
    form: 'list',
    start: Int32Array.from([100, 3000, 400]),
    end: Int32Array.from([200, 3100, 500]),
    id: categoryColumn(['t1', 't1', 't2']),
    itemStart: Int32Array.from([100, 100, 400]),
    itemEnd: Int32Array.from([3100, 3100, 500]),
  };
}
