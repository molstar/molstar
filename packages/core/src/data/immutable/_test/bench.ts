/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small timing loop for the heavy lane's ratios (README §8): runs a case until it has taken `minMs`, then takes the
// fastest of `rounds` such runs, in nanoseconds an operation. The fastest, not the median: a major collection landing in
// a round of an allocating case (100,000 persistent inserts, after the lane's model runs filled the heap) doubled
// a median once (2026-11-18); both sides are timed alike.

export interface Timed {
  nsPerOp: number;
}

/** `run()` does `ops` operations; returns the fastest round's ns an operation. */
export function time(run: () => unknown, ops: number, minMs = 25, rounds = 7): number {
  let reps = 1;
  // warm and size the repetitions
  for (;;) {
    const t0 = performance.now();
    for (let r = 0; r < reps; r++) sink(run());
    const dt = performance.now() - t0;
    if (dt >= minMs) break;
    reps = Math.max(reps * 2, Math.ceil((reps * minMs) / Math.max(dt, 0.01)));
  }
  const samples: number[] = [];
  for (let i = 0; i < rounds; i++) {
    const t0 = performance.now();
    for (let r = 0; r < reps; r++) sink(run());
    samples.push(((performance.now() - t0) * 1e6) / (reps * ops));
  }
  return Math.min(...samples);
}

let sunk: unknown;
function sink(x: unknown) {
  sunk = x;
}

export function sinkValue() {
  return sunk;
}
