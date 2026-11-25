/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Coverage, GridDense, SeqData, SeqPair, ShapeInfo } from '@molstar/seqx-core/data';

// A generated Hi-C over a toy assembly (Features/Region reads): a symmetric contact
// map as counts, the upper triangle's block pairs served (a chromosome with itself whole, a pair of chromosomes once),
// with two levels of its own (`hicSource`) or at the finest level alone (`hicOneLevel`: the store's summaries give it
// the others, and the whole-assembly level). Deterministic on core's `PCG`, small enough to be generated at open.

export interface HicChromosome {
  readonly name: string;
  /** Positions. */
  readonly length: number;
}

/**
 * A toy assembly's six chromosomes, 1,376 positions in all. A position is one bin of the map (100 kb of a real genome's
 * Hi-C): level 1 is a cell, so the store's summaries give a one-level source its other levels (summaries are for sources
 * whose one level is 1), and a whole-assembly level.
 */
export const HicChromosomes: readonly HicChromosome[] = [
  { name: 'chr1', length: 384 },
  { name: 'chr2', length: 320 },
  { name: 'chr3', length: 256 },
  { name: 'chr4', length: 192 },
  { name: 'chr5', length: 128 },
  { name: 'chrX', length: 96 },
];

/** The finest level's bin and the coarser one's (positions a cell a side). */
export const HicBin = 1;
export const HicCoarseBin = 4;
/** The whole-assembly system's name. */
export const HicWhole = 'all';
/** The assembly's name as a view writes it. */
export const HicAssemblyName = 'ToyHiC';

export interface HicOptions {
  readonly chromosomes?: readonly HicChromosome[];
  readonly seed?: number;
}

/** A chromosome's cells at the finest level: counts falling with the distance, domains along the diagonal. */
function selfPlane(n: number, rng: PCG): Float32Array {
  const v = new Float32Array(n * n);
  // domains (TADs): runs of 8 to 24 bins whose contacts are stronger inside
  const domain = new Int32Array(n);
  for (let i = 0, d = 0; i < n; d++) {
    const len = 8 + Math.floor(rng.float() * 17);
    for (let k = 0; k < len && i < n; k++) domain[i++] = d;
  }
  for (let i = 0; i < n; i++)
    for (let j = i; j < n; j++) {
      const d = j - i;
      const base = 600 / (1 + d) ** 1.1;
      const inside = domain[i] === domain[j] ? 1.8 : 1;
      const c = Math.floor(base * inside * (0.7 + 0.6 * rng.float()));
      v[i * n + j] = c;
      v[j * n + i] = c;
    }
  return v;
}

/** Two chromosomes' cells: sparse noise, and a hot spot (a translocation) between chr2's and chr5's arms. */
function crossPlane(nx: number, ny: number, hot: boolean, rng: PCG): Float32Array {
  const v = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      let c = Math.floor(rng.float() * 3);
      if (hot) {
        const dx = (i - nx * 0.7) / 6;
        const dy = (j - ny * 0.4) / 6;
        c += Math.floor(40 * Math.exp(-(dx * dx + dy * dy)) * (0.7 + 0.6 * rng.float()));
      }
      v[i * ny + j] = c;
    }
  return v;
}

/** A plane at `factor` bins of the finest level a side, each cell the sum of the cells below. */
function sumPlane(v: Float32Array, nx: number, ny: number, factor: number): GridDense {
  const mx = Math.ceil(nx / factor);
  const my = Math.ceil(ny / factor);
  const out = new Float32Array(mx * my);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) out[Math.floor(i / factor) * my + Math.floor(j / factor)] += v[i * ny + j];
  return { form: 'dense', x: 0, y: 0, nx: mx, ny: my, value: out };
}

/**
 * The symmetric contact map of a toy assembly. `levels: 2` serves the finest level (`HicBin`) and a coarser one
 * (`HicCoarseBin`, four bins to one, the cells summed); `1` serves the finest alone. The coverage lists a whole-assembly
 * system (`all`, the chromosomes end to end) that the store's summaries fill in.
 */
function hic(levels: 1 | 2, options: HicOptions): SeqData {
  const chromosomes = options.chromosomes ?? HicChromosomes;
  const rng = new PCG(options.seed ?? 5);
  const total = chromosomes.reduce((a, c) => a + c.length, 0);
  const factors = levels === 2 ? [1, HicCoarseBin / HicBin] : [1];
  const shapes: ShapeInfo[] = chromosomes.map((c, i) => ({
    system: c.name,
    shape: 'grid' as const,
    y: chromosomes.slice(i).map((o) => o.name),
    symmetric: true,
    aggregate: 'sum' as const,
    ...(i === 0 && { whole: { system: HicWhole, sequences: chromosomes.map((o) => o.name) } }),
  }));
  shapes.push({ system: HicWhole, shape: 'grid', y: [HicWhole], symmetric: true, aggregate: 'sum' });
  const coverage: Coverage = {
    systems: [...chromosomes.map((c) => ({ name: c.name, length: c.length })), { name: HicWhole, length: total }],
    shapes,
    levels: factors.map((f) => f * HicBin),
  };
  const pairs: SeqPair[] = [];
  chromosomes.forEach((cx, i) => {
    for (let j = i; j < chromosomes.length; j++) {
      const cy = chromosomes[j];
      const nx = Math.ceil(cx.length / HicBin);
      const ny = Math.ceil(cy.length / HicBin);
      const cells = i === j ? selfPlane(nx, rng) : crossPlane(nx, ny, cx.name === 'chr2' && cy.name === 'chr5', rng);
      for (const f of factors) {
        const grid =
          f === 1 ? ({ form: 'dense', x: 0, y: 0, nx, ny, value: cells } as GridDense) : sumPlane(cells, nx, ny, f);
        pairs.push({
          request: {
            system: cx.name,
            start: 0,
            end: cx.length,
            level: f * HicBin,
            shapes: ['grid'],
            y: { system: cy.name, start: 0, end: cy.length },
          },
          answer: { ok: true, tables: { grid } },
        });
      }
    }
  });
  return { coverage, pairs };
}

/** A Hi-C with two levels of its own (`HicBin` and `HicCoarseBin`). */
export function hicSource(options: HicOptions = {}): SeqData {
  return hic(2, options);
}

/** The same Hi-C at its finest level alone: its other levels are the store's summaries, its twin by the same seed. */
export function hicOneLevel(options: HicOptions = {}): SeqData {
  return hic(1, options);
}
