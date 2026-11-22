/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { MatrixBins, MatrixItems, SeqData } from '@molstar/seqx-core/data';

// A generated genotype matrix (Features/Region reads): 40 samples over one
// chromosome's variants, a code each (`0` no copy, `1` one, `2` both, `255` uncalled), with one long deletion that a
// few samples carry across the SNVs inside it. Served at the items' level and at a coarse one (a bin's majority class
// and its share). Deterministic on core's `PCG`.

export const GenotypeSystem = 'chr';
export const GenotypeLength = 1_048_576;
export const GenotypeSamples: readonly string[] = Array.from(
  { length: 40 },
  (_, i) => `HG${String(i + 1).padStart(5, '0')}`,
);
/** The coarse level's bin width (positions). */
export const GenotypeBin = 2048;
/** The long deletion's extent and the samples (indices) that carry it. */
export const GenotypeDeletion = { start: 200_000, end: 260_000, carriers: [3, 4, 5, 11, 12, 27] } as const;
const Variants = 800;

/** The variants (sorted, the deletion among them) and a code per variant and sample, item-major. */
export function genotypeItems(seed = 21): MatrixItems & { readonly deletion: number } {
  const rng = new PCG(seed);
  const S = GenotypeSamples.length;
  const rows: { start: number; end: number; id: string }[] = [];
  for (let i = 0; i < Variants; i++) {
    const start = 1000 + Math.floor((i / Variants) * (GenotypeLength - 2000)) + Math.floor(rng.float() * 600);
    rows.push({ start, end: start + 1, id: `rs${100_000 + i}` });
  }
  rows.push({ start: GenotypeDeletion.start, end: GenotypeDeletion.end, id: 'del1' });
  rows.sort((a, b) => a.start - b.start || b.end - a.end);
  const n = rows.length;
  const codes = new Uint8Array(n * S);
  for (let k = 0; k < n; k++) {
    // a SNV's allele frequency, the samples drawn from it; a few uncalled
    const af = rng.float() ** 2 * 0.6;
    for (let s = 0; s < S; s++) {
      const r = rng.float();
      codes[k * S + s] = r < 0.02 ? 255 : r < 0.02 + af * af ? 2 : r < 0.02 + af * af + 2 * af * (1 - af) ? 1 : 0;
    }
    if (rows[k].id === 'del1')
      for (let s = 0; s < S; s++) codes[k * S + s] = GenotypeDeletion.carriers.includes(s as never) ? 1 : 0;
  }
  const deletion = rows.findIndex((r) => r.id === 'del1');
  return {
    form: 'items',
    start: Int32Array.from(rows, (r) => r.start),
    end: Int32Array.from(rows, (r) => r.end),
    states: S,
    codes,
    id: categoryColumn(rows.map((r) => r.id)),
    deletion,
  };
}

/** A matrix of codes at a coarser level: each bin's majority class per sample and its share (0 to 255). */
export function binCodes(items: MatrixItems, length: number, width: number): MatrixBins {
  const S = items.states;
  const n = Math.ceil(length / width);
  const counts = new Float32Array(n * S * 4);
  const called = new Float32Array(n * S);
  const codes = items.codes!;
  for (let k = 0; k < items.start.length; k++) {
    // an item counts in the bin holding its start (a long deletion once, where it starts)
    const b = Math.floor(items.start[k] / width);
    for (let s = 0; s < S; s++) {
      const c = codes[k * S + s];
      if (c > 3) continue;
      counts[(b * S + s) * 4 + c]++;
      called[b * S + s]++;
    }
  }
  const classes = new Uint8Array(n * S).fill(255);
  const shares = new Uint8Array(n * S);
  for (let i = 0; i < n * S; i++) {
    if (called[i] === 0) continue;
    let best = 0;
    for (let c = 1; c < 4; c++) if (counts[i * 4 + c] > counts[i * 4 + best]) best = c;
    classes[i] = best;
    shares[i] = Math.round((255 * counts[i * 4 + best]) / called[i]);
  }
  return { form: 'bins', start: 0, states: S, classes, shares };
}

/** The genotypes as a source: the items (level 1) and the bins (`GenotypeBin`). */
export function genotypesSource(seed = 21): SeqData {
  const items = genotypeItems(seed);
  const { deletion: _d, ...table } = items;
  const request = (level: number) => ({
    system: GenotypeSystem,
    start: 0,
    end: GenotypeLength,
    level,
    shapes: ['matrix' as const],
  });
  return {
    coverage: {
      systems: [{ name: GenotypeSystem, length: GenotypeLength }],
      shapes: [{ system: GenotypeSystem, shape: 'matrix', states: [...GenotypeSamples] }],
      levels: [1, GenotypeBin],
    },
    pairs: [
      { request: request(1), answer: { ok: true, tables: { matrix: table } } },
      {
        request: request(GenotypeBin),
        answer: { ok: true, tables: { matrix: binCodes(table, GenotypeLength, GenotypeBin) } },
      },
    ],
  };
}
