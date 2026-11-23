/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Coverage, GridDense, GridSparse, SeqData, SeqPair, ShapeInfo } from '@molstar/seqx-core/data';
import type { ToyAssembly } from './assemblies.js';
import { Concern, stream } from './streams.js';

// Step 8 of the Genomes generator (examples.md, Genomes; EC-223): a Hi-C over toy A's chr1 to chr3, symmetric, counts
// in cells of 4,096 and 16,384 positions a side: within a chromosome falling with the distance and stronger inside
// domains (dense, the whole square), between two chromosomes sparse noise and one hotspot between chr1 and chr3 (a
// translocation's mark; sparse, the rest 0). The upper triangle's sequence pairs are served; the lower is the
// mirror. Why not binned from the ring's links: the links are the ring extension's step (7) and over GRCh38, not toy A.

/** The contact map's levels (positions a cell a side). */
export const ContactLevels: readonly number[] = [4096, 16_384];
/** The hotspot between chr1 and chr3, in positions. */
export const Hotspot: readonly [chr1: number, chr3: number] = [230_000, 40_000];

/** A chromosome's cells at the finest level: counts by distance, domains of 6 to 16 cells along the diagonal. */
function selfPlane(n: number, seed: number, k: number): Float32Array {
  const rng = stream(seed, Concern.grids, k);
  const v = new Float32Array(n * n);
  const domain = new Int32Array(n);
  for (let i = 0, d = 0; i < n; d++) {
    const len = 6 + Math.floor(rng.float() * 11);
    for (let j = 0; j < len && i < n; j++) domain[i++] = d;
  }
  for (let i = 0; i < n; i++)
    for (let j = i; j < n; j++) {
      const c = Math.floor(
        (400 / (1 + j - i) ** 1.1) * (domain[i] === domain[j] ? 1.8 : 1) * (0.7 + 0.6 * rng.float()),
      );
      v[i * n + j] = c;
      v[j * n + i] = c;
    }
  return v;
}

/** Two chromosomes' cells: sparse noise, and the hotspot when `hot` gives its cell. */
function crossPlane(nx: number, ny: number, hot: readonly [number, number] | undefined, seed: number, k: number) {
  const rng = stream(seed, Concern.grids, k);
  const v = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      let c = rng.float() < 0.15 ? 1 : 0;
      if (hot) c += Math.floor(60 * Math.exp(-((i - hot[0]) ** 2 + (j - hot[1]) ** 2) / 4));
      v[i * ny + j] = c;
    }
  return v;
}

/** A plane summed to `factor` cells a side. */
function sumPlane(v: Float32Array, nx: number, ny: number, factor: number) {
  const mx = Math.ceil(nx / factor);
  const my = Math.ceil(ny / factor);
  const out = new Float32Array(mx * my);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) out[Math.floor(i / factor) * my + Math.floor(j / factor)] += v[i * ny + j];
  return { nx: mx, ny: my, value: out };
}

/** A plane as a sparse grid (its non-zero cells). */
function sparse(nx: number, ny: number, value: Float32Array): GridSparse {
  const keys: number[] = [];
  for (let k = 0; k < nx * ny; k++) if (value[k] !== 0) keys.push(k);
  return {
    form: 'sparse',
    x: 0,
    y: 0,
    xBin: Int32Array.from(keys, (k) => Math.floor(k / ny)),
    yBin: Int32Array.from(keys, (k) => k % ny),
    value: Float32Array.from(keys, (k) => value[k]),
    fill: 0,
  };
}

/** The contact map of toy A's chr1 to chr3 at both levels. */
export function contactMap(a: ToyAssembly, seed: number): SeqData {
  const chrs = a.sequences.filter((s) => /^chr[123]$/.test(s.name));
  const [fine, coarse] = ContactLevels;
  const shapes: ShapeInfo[] = chrs.map((c, i) => ({
    system: c.name,
    shape: 'grid' as const,
    y: chrs.slice(i).map((o) => o.name),
    symmetric: true,
    aggregate: 'sum' as const,
  }));
  const coverage: Coverage = {
    systems: chrs.map((c) => ({ name: c.name, length: c.length })),
    shapes,
    levels: [...ContactLevels],
  };
  const pairs: SeqPair[] = [];
  let k = 0;
  chrs.forEach((cx, i) => {
    for (let j = i; j < chrs.length; j++) {
      const cy = chrs[j];
      const nx = Math.ceil(cx.length / fine);
      const ny = Math.ceil(cy.length / fine);
      const hot =
        cx.name === 'chr1' && cy.name === 'chr3'
          ? ([Math.floor(Hotspot[0] / fine), Math.floor(Hotspot[1] / fine)] as const)
          : undefined;
      const cells = i === j ? selfPlane(nx, seed, k++) : crossPlane(nx, ny, hot, seed, k++);
      for (const level of ContactLevels) {
        const p = level === fine ? { nx, ny, value: cells } : sumPlane(cells, nx, ny, coarse / fine);
        const grid: GridDense | GridSparse =
          i === j ? { form: 'dense', x: 0, y: 0, nx: p.nx, ny: p.ny, value: p.value } : sparse(p.nx, p.ny, p.value);
        pairs.push({
          request: {
            system: cx.name,
            start: 0,
            end: cx.length,
            level,
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
