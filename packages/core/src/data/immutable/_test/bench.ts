/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small timing loop for the heavy lane's ratios (README §8): runs a case until it has taken `minMs`, then takes the
// median of `rounds` such runs, in nanoseconds an operation.

export interface Timed {
  nsPerOp: number;
}

/** `run()` does `ops` operations; returns the median ns an operation. */
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
  samples.sort((a, b) => a - b);
  return samples[rounds >> 1];
}

let sunk: unknown;
function sink(x: unknown) {
  sunk = x;
}

export function sinkValue() {
  return sunk;
}
