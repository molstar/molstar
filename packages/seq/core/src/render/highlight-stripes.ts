/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LocationPart } from '@molstar/seqx-core/query';

// Overlapping highlight channels stay apart (_spec/README.md "Render", Emphasis): the runs of positions over which the
// set of channels lighting is the same, each channel's with its stripe: the slot among the channels lighting there (the
// last stripe shared past the maximum) and how many stripes the band is split into. Shared by the tracks' bands and the
// grid's rows. Allocates nothing once warm: the channels' sorted boundaries are merged, never sorted.

/** Per channel: the next boundary of its ranges, whether it lights, its open run's start, slot, stripes and kind. */
let Cursor = new Int32Array(8);
let Lit = new Uint8Array(8);
let Start = new Int32Array(8);
let Slot = new Int32Array(8);
let Stripes = new Int32Array(8);
let First = new Uint8Array(8);
let Opened = new Uint8Array(8);
let Closed = new Uint8Array(8);

/** What `walkStripes` reports for a run: the channel (its index in `parts`), its `[start, end)`, its slot and count. */
export type StripeRun = (
  channel: number,
  start: number,
  end: number,
  slot: number,
  count: number,
  /** 1: the run begins where one of the channel's ranges does (not where its slot or count changed). */
  first: number,
  /** 1: the run ends where one of the channel's ranges does (its real end: a squircle's round end there). */
  last: number,
) => void;

/**
 * The runs of `n` channels' parts (half-open sorted ranges), at most `max` stripes: `emit` is given each channel's run
 * where its slot and the stripe count stay the same, in the order the runs end (channel order within a position).
 */
export function walkStripes(parts: readonly (LocationPart | undefined)[], n: number, max: number, emit: StripeRun) {
  if (Cursor.length < n) {
    Cursor = new Int32Array(n * 2);
    Lit = new Uint8Array(n * 2);
    Start = new Int32Array(n * 2);
    Slot = new Int32Array(n * 2);
    Stripes = new Int32Array(n * 2);
    First = new Uint8Array(n * 2);
    Opened = new Uint8Array(n * 2);
    Closed = new Uint8Array(n * 2);
  }
  for (let c = 0; c < n; c++) {
    Cursor[c] = 0;
    Slot[c] = -1;
  }
  for (;;) {
    let at = 0x7fffffff;
    let any = false;
    for (let c = 0; c < n; c++) {
      const r = parts[c]!.ranges;
      if (Cursor[c] < r.length && r[Cursor[c]] <= at) {
        at = r[Cursor[c]];
        any = true;
      }
    }
    if (!any) break;
    // who lights [at, next): a channel lights after it consumed a range's start, not after its end
    let present = 0;
    for (let c = 0; c < n; c++) {
      const r = parts[c]!.ranges;
      Opened[c] = 0;
      Closed[c] = 0;
      while (Cursor[c] < r.length && r[Cursor[c]] === at) {
        if ((Cursor[c] & 1) === 0) Opened[c] = 1;
        else Closed[c] = 1;
        Cursor[c]++;
      }
      const lit = (Cursor[c] & 1) === 1;
      Lit[c] = lit ? 1 : 0;
      if (lit) present++;
    }
    const count = Math.min(present, max);
    let rank = 0;
    for (let c = 0; c < n; c++) {
      const slot = Lit[c] ? Math.min(rank++, max - 1) : -1;
      if (Slot[c] >= 0 && (slot !== Slot[c] || count !== Stripes[c])) {
        emit(c, Start[c], at, Slot[c], Stripes[c], First[c], Closed[c] === 1 && Lit[c] === 0 ? 1 : 0);
        Slot[c] = -1;
      }
      if (slot >= 0 && Slot[c] < 0) {
        Start[c] = at;
        Slot[c] = slot;
        Stripes[c] = count;
        First[c] = Opened[c];
      }
    }
  }
}
