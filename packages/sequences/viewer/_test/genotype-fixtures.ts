/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS2's genotype matrix: samples as a matrix's states over a sequence's positions, genotype codes (0, 1, 2) as values,
// and a sample table keyed by state (`population`, `depth`), for the dense-matrix tests and heap guard.

import { createSeqViewSpec } from '@molstar/seq-core/builder';

/** Sample q's genotype at position i: mostly 0, some 1 and 2, varying along both axes. */
export const genotype = (i: number, q: number): number =>
  (i * 7 + q * 13) % 11 === 0 ? 2 : (i * 3 + q * 5) % 4 === 0 ? 1 : 0;
export const POPULATIONS = ['AFR', 'EUR', 'EAS'];
/** Heterozygous and homozygous alternative genotypes colored; the reference genotype (0) not. */
const CARRIERS = {
  theme: 'field_palette',
  params: { palette: { kind: 'categorical', colors: { '1': '#e8a33d', '2': '#c2412d' } } },
};

/** `holes`: genotype 0 uncolored (a categorical palette for 1 and 2 only), so reduced rows fade by their share. */
export function genotypeDoc(o: {
  samples: number;
  positions: number;
  layer?: Record<string, unknown>;
  holes?: boolean;
  /** Sample q's value at position i (default `genotype`). */
  value?: (i: number, q: number) => number;
  /** Letters tracks after the matrix, so the section's content is taller than its box and scrolls. */
  below?: number;
  /** A second region, stacked under the first, with its own matrix of the same states in the same layer. */
  stacked?: boolean;
}) {
  const b = createSeqViewSpec({ title: 'Genotypes' });
  const s = b.data.sequence({ ref: 's', value: 'ACGT'.repeat(Math.ceil(o.positions / 4)).slice(0, o.positions) });
  const states = Array.from({ length: o.samples }, (_, q) => `n${q}`);
  const index = Array.from({ length: o.positions }, (_, i) => i);
  const at = o.value ?? genotype;
  const value = { index, value: index.map((i) => states.map((_, q) => at(i, q))) };
  const m = b.data.annotation({ ref: 'g', geometry: 'matrix', target: s, states, value } as never);
  b.data.annotation({
    ref: 'samples',
    geometry: 'rows',
    value: {
      state: states,
      population: states.map((_, q) => POPULATIONS[q % 3]),
      depth: states.map((_, q) => (q * 37) % 101),
    },
  } as never);
  const view = b.view({ ref: 'v', label: 'V', ...(o.stacked ? { arrangement: 'stacked' } : {}) } as never);
  const v = view.block({ sequence: s });
  const data: unknown[] = [m];
  if (o.stacked) {
    const s2 = b.data.sequence({ ref: 's2', value: 'TGCA'.repeat(Math.ceil(o.positions / 4)).slice(0, o.positions) });
    data.push(b.data.annotation({ ref: 'g2', geometry: 'matrix', target: s2, states, value } as never));
    view.block({ sequence: s2 });
  }
  v.track({ label: 'G' }).matrix({
    ref: 'geno',
    data: o.stacked ? data : m,
    ...(o.holes ? { color: CARRIERS } : {}),
    ...(o.layer ?? {}),
  } as never);
  for (let k = 0; k < (o.below ?? 0); k++) v.track({ label: `L${k}` }).letters();
  return b.getDocument();
}
