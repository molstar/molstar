/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ToyAssembly } from './assemblies.js';
import { NGap } from './assemblies.js';
import { Concern, decimal, stream } from './streams.js';

// Step 4 of the Genomes generator (examples.md, Genomes; EC-215): one bedGraph over toy A, 100 bp bins, a
// TPM-like value from 0.1 to 10^4 (smooth waves in log space plus noise, so a linear scale shows a few peaks and a log
// scale the rest), and stretches with no line (missing values: chr1's run of N, and 10 kb of A's chr2). Sorted by
// sequence, then start: what tabix indexes. Why not over toy B too (the design's "both toy genomes"): B's `1` is A's
// `chr1` by the names rule's `chr` prefix, so one file naming both would join their values onto one sequence.

/** The signal's bin (positions a line). */
export const SignalBin = 100;
/** A's chr2 stretch with no value (missing, never 0). */
export const SignalGap: readonly [number, number] = [80_000, 90_000];
/** The value's bound: TPM-like, 0 to 10^4. */
export const SignalMax = 10_000;

/** Whether `[s, e)` of `chr` has no value. */
function missing(assembly: string, chr: string, s: number, e: number): boolean {
  if (assembly !== 'Toy A') return false;
  const inside = (r: readonly [number, number]) => s < r[1] && e > r[0];
  return (chr === 'chr1' && inside(NGap)) || (chr === 'chr2' && inside(SignalGap));
}

/** The bedGraph of `assemblies`' sequences, each from its own stream. */
export function bedGraph(assemblies: readonly ToyAssembly[], seed: number): string {
  let out = '';
  let k = 0;
  for (const a of assemblies)
    for (const s of a.sequences) {
      const rng = stream(seed, Concern.signal, k++);
      const phase = rng.float() * 2 * Math.PI;
      for (let start = 0; start < s.length; start += SignalBin) {
        const end = Math.min(s.length, start + SignalBin);
        if (missing(a.name, s.name, start, end)) continue;
        // log10(1 + v) in [0, 4]: two waves and noise, the high peaks rare
        const w =
          1.6 +
          1.1 * Math.sin((2 * Math.PI * start) / 37_000 + phase) +
          0.9 * Math.sin((2 * Math.PI * start) / 9_100 + 2 * phase) ** 7 +
          0.35 * (rng.float() - 0.5);
        // at least 0.1: a log scale reads every value (a zero would be missing there)
        const v = Math.min(SignalMax, Math.max(0.1, 10 ** Math.max(0, w) - 1));
        out += `${s.name}\t${start}\t${end}\t${decimal(v, 1)}\n`;
      }
    }
  return out;
}

/** A bedGraph line's sequence and extent for the index (columns 1 to 3, 0-based half-open). */
export const bedExtent = (c: readonly string[]) =>
  c.length >= 3 ? { seq: c[0], begin: Number(c[1]), end: Number(c[2]) } : undefined;
