/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { MatrixSummary, MatrixItems, SeqData } from '@molstar/seqx-core/data';

// Generated methylation (Features/Region reads): a fraction in 0 to 1 per CpG and
// sample, a CpG island where it falls to near zero, and a share of values missing (a sample that did not cover a CpG:
// NaN), served as the items and as bins of their mean. Deterministic on core's `PCG`.

export const MethylationSystem = 'chr';
export const MethylationLength = 524_288;
export const MethylationSamples: readonly string[] = [
  'liver 1',
  'liver 2',
  'brain 1',
  'brain 2',
  'blood 1',
  'blood 2',
  'lung 1',
  'lung 2',
];
export const MethylationBin = 1024;
/** The island's extent (CpGs there are unmethylated). */
export const MethylationIsland = { start: 240_000, end: 262_000 } as const;
const Cpgs = 2000;

/** Each CpG a position, its values a sample each (item-major), about a tenth missing. */
export function methylationItems(seed = 31): MatrixItems & { readonly missing: number } {
  const rng = new PCG(seed);
  const S = MethylationSamples.length;
  const pos = new Int32Array(Cpgs);
  for (let i = 0; i < Cpgs; i++) {
    const inIsland = i % 5 === 0;
    pos[i] = inIsland
      ? MethylationIsland.start + Math.floor(rng.float() * (MethylationIsland.end - MethylationIsland.start))
      : Math.floor(rng.float() * (MethylationLength - 1));
  }
  pos.sort();
  const value = new Float32Array(Cpgs * S);
  let missing = 0;
  for (let i = 0; i < Cpgs; i++) {
    const island = pos[i] >= MethylationIsland.start && pos[i] < MethylationIsland.end;
    // brain differs from the others at the island's edge: a tissue's mark
    for (let s = 0; s < S; s++) {
      const base = island ? 0.05 : 0.8;
      const v = Math.min(1, Math.max(0, base + (rng.float() - 0.5) * 0.3 - (s >> 1 === 1 && !island ? 0.25 : 0)));
      const gone = rng.float() < 0.1;
      if (gone) missing++;
      value[i * S + s] = gone ? Number.NaN : v;
    }
  }
  return { form: 'items', start: pos, end: pos.map((p) => p + 2), states: S, value, missing };
}

/** The values at a coarser level: each bin's mean per sample (NaN where no CpG was covered), and the counts. */
export function binValues(items: MatrixItems, length: number, width: number): MatrixSummary {
  const S = items.states;
  const n = Math.ceil(length / width);
  const sum = new Float64Array(n * S);
  const count = new Float32Array(n * S);
  for (let k = 0; k < items.start.length; k++) {
    const b = Math.floor(items.start[k] / width);
    for (let s = 0; s < S; s++) {
      const v = items.value![k * S + s];
      if (Number.isNaN(v)) continue;
      sum[b * S + s] += v;
      count[b * S + s]++;
    }
  }
  const mean = new Float32Array(n * S);
  for (let i = 0; i < n * S; i++) mean[i] = count[i] === 0 ? Number.NaN : sum[i] / count[i];
  return { form: 'summary', start: 0, states: S, mean, count };
}

/** The methylation as a source: items (level 1) and bins (`MethylationBin`). */
export function methylationSource(seed = 31): SeqData {
  const { missing: _m, ...items } = methylationItems(seed);
  const request = (level: number) => ({
    system: MethylationSystem,
    start: 0,
    end: MethylationLength,
    level,
    shapes: ['matrix' as const],
  });
  return {
    coverage: {
      systems: [{ name: MethylationSystem, length: MethylationLength }],
      shapes: [{ system: MethylationSystem, shape: 'matrix', states: [...MethylationSamples] }],
      levels: [1, MethylationBin],
    },
    pairs: [
      { request: request(1), answer: { ok: true, tables: { matrix: items } } },
      {
        request: request(MethylationBin),
        answer: { ok: true, tables: { matrix: binValues(items, MethylationLength, MethylationBin) } },
      },
    ],
  };
}
