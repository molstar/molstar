/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hash3, PCG } from '@molstar/core/data/util/hash-functions';

// The Genomes dataset's randomness (examples.md, Genomes, The generator: outline): core's `PCG` alone, a stream per
// `(seed, concern, index)` seeded by core's `hash3`, so adding a concern, or drawing more from one, never moves another
// concern's bytes. No clock, no `Math.random`.

/** The generator's concerns, each its own streams (7, the links, is the ring extension's). */
export const Concern = {
  letters: 1,
  alignment: 2,
  genes: 3,
  signal: 4,
  variants: 5,
  methylation: 6,
  links: 7,
  grids: 8,
} as const;

/** The default seed (`genomes.seed`): one fixed value, so the fixtures are byte-identical on every run. */
export const GenomesSeed = 2026;

/** The stream of `(seed, concern, index)`: a `PCG` seeded by their hash. */
export function stream(seed: number, concern: number, index: number): PCG {
  return new PCG(hash3(seed, concern, index) >>> 0);
}

/** An integer in `[lo, hi)`. */
export function between(rng: PCG, lo: number, hi: number): number {
  return lo + Math.floor(rng.float() * (hi - lo));
}

/** A number written with `digits` decimals at most, no trailing zeros (`1.50` → `1.5`, `2.00` → `2`). */
export function decimal(v: number, digits: number): string {
  const s = v.toFixed(digits);
  return s.includes('.') ? s.replace(/0+$/, '').replace(/\.$/, '') : s;
}
