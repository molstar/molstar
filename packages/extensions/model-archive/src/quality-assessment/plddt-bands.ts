/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// pLDDT's confidence bands, AlphaFold DB's: the color theme's and the PAE plot's strip.

export const PlddtNoScore = 0xaaaaaa;

/** From the highest: a score above `above` is in the band. */
export const PlddtBands = [
  { label: 'Very high', range: '> 90', above: 90, color: 0x0053d6 },
  { label: 'Confident', range: '70–90', above: 70, color: 0x65cbf3 },
  { label: 'Low', range: '50–70', above: 50, color: 0xffdb13 },
  { label: 'Very low', range: '< 50', above: Number.NEGATIVE_INFINITY, color: 0xff7d45 },
] as const;

/** A score's band color (0–100); a negative or missing score has none. */
export function plddtColor(score: number): number {
  if (!(score >= 0)) return PlddtNoScore;
  for (const b of PlddtBands) if (score > b.above) return b.color;
  return PlddtBands[PlddtBands.length - 1].color;
}
