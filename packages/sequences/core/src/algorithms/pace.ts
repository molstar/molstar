/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A long run's pacing: it counts basic steps and yields to the event loop every `SLICE` of them, checking an abort signal on
// both sides, so a newer request stops it at the next slice (a synchronous loop could not observe the signal).

/** Basic steps (one pair of cells compared, one matrix cell scanned) between yields to the event loop. */
const SLICE = 1 << 21;

/** Counts basic steps; `due` says when to yield, `yieldNow` yields and checks the signal on both sides. */
export class Pacer {
  private work = 0;
  private readonly signal: AbortSignal | undefined;
  constructor(signal: AbortSignal | undefined) {
    this.signal = signal;
    signal?.throwIfAborted();
  }
  due(steps: number): boolean {
    this.work += steps;
    return this.work >= SLICE;
  }
  async yieldNow(): Promise<void> {
    this.work = 0;
    this.signal?.throwIfAborted();
    await new Promise<void>((r) => setTimeout(r, 0));
    this.signal?.throwIfAborted();
  }
}
