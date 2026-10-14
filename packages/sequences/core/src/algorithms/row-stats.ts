/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Row statistics and redundancy over an alignment's rows (seqspec:format#computed-nodes, seqspec:model#providers): the
// kernel providers `row_stats` and `redundancy`. Pure, over typed arrays, generic over any alphabet (a cell is a letter
// or the gap code). Identity between two rows is equal letters over the columns where both have a letter, NaN with
// none. `redundancy` is greedy clustering, rows × representatives × columns at most, so it is bounded
// (`Redundancy.MaxCells`) and takes an abort signal, which it honors between slices of work (`Pacer`, pace.ts).

import { Pacer } from "./pace.js";

/** An alignment's cells, row-major `rows × columns`; `gap` is the code of a cell with no letter. */
export interface RowCodes {
  readonly codes: Uint8Array;
  readonly rows: number;
  readonly columns: number;
  readonly gap: number;
}

/** The size bound of `redundancy`: rows × rows × columns cells compared at most. */
export const Redundancy = { MaxCells: 2 ** 31 } as const;

/** Thrown before any work when `redundancy` would compare more than `Redundancy.MaxCells` cells. */
export class RedundancySizeError extends Error {
  readonly cells: number;
  constructor(rows: number, columns: number) {
    const cells = rows * rows * columns;
    super(
      `clustering ${rows} rows of ${columns} columns compares up to ${cells} cells; the limit is ${Redundancy.MaxCells}.`,
    );
    this.name = "RedundancySizeError";
    this.cells = cells;
  }
}

/** Per row: its identity to row `reference`, the share of its columns that hold a letter, and its letter count. */
export interface RowStats {
  readonly identity: Float64Array;
  readonly coverage: Float64Array;
  readonly length: Float64Array;
}

/** Identity of rows `a` and `b` over the columns where both have a letter (NaN with none). */
export function identityOf(c: RowCodes, a: number, b: number): number {
  const [oa, ob] = [a * c.columns, b * c.columns];
  let both = 0;
  let same = 0;
  for (let k = 0; k < c.columns; k++) {
    const x = c.codes[oa + k];
    const y = c.codes[ob + k];
    if (x !== c.gap && y !== c.gap) {
      both++;
      if (x === y) same++;
    }
  }
  return both === 0 ? Number.NaN : same / both;
}

/** `identity` to row `reference`, `coverage` (letters over all columns) and `length` (letters) of every row. */
export function rowStats(c: RowCodes, reference: number): RowStats {
  const identity = new Float64Array(c.rows);
  const coverage = new Float64Array(c.rows);
  const length = new Float64Array(c.rows);
  for (let r = 0; r < c.rows; r++) {
    let letters = 0;
    for (let k = 0; k < c.columns; k++) if (c.codes[r * c.columns + k] !== c.gap) letters++;
    length[r] = letters;
    coverage[r] = c.columns === 0 ? Number.NaN : letters / c.columns;
    identity[r] = identityOf(c, r, reference);
  }
  return { identity, coverage, length };
}

/** Per row: the row index of its cluster's representative, and 1 where it is one. */
export interface Clusters {
  readonly cluster: Int32Array;
  readonly representative: Uint8Array;
}

/**
 * Greedy clustering at an identity `threshold`: rows are taken in `order` (row indices; `null`: row order, rows left out
 * of a given order are taken after it, in row order); a row joins the first representative so far it is at least
 * `threshold` identical to, else it is one. Past `Redundancy.MaxCells` it throws `RedundancySizeError`.
 */
export async function redundancy(
  c: RowCodes,
  threshold: number,
  order: ArrayLike<number> | null,
  signal?: AbortSignal,
): Promise<Clusters> {
  if (c.rows * c.rows * c.columns > Redundancy.MaxCells) throw new RedundancySizeError(c.rows, c.columns);
  const pace = new Pacer(signal);
  const cluster = new Int32Array(c.rows).fill(-1);
  const representative = new Uint8Array(c.rows);
  const taken = new Uint8Array(c.rows);
  const sequence: number[] = [];
  for (let k = 0; order && k < order.length; k++) {
    const r = order[k];
    if (r >= 0 && r < c.rows && !taken[r]) {
      taken[r] = 1;
      sequence.push(r);
    }
  }
  for (let r = 0; r < c.rows; r++) if (!taken[r]) sequence.push(r);
  const reps: number[] = [];
  for (const r of sequence) {
    let home = -1;
    for (let k = 0; k < reps.length && home < 0; k++) {
      const x = identityOf(c, r, reps[k]);
      if (x >= threshold) home = reps[k];
      if (pace.due(c.columns)) await pace.yieldNow();
    }
    if (home < 0) {
      home = r;
      representative[r] = 1;
      reps.push(r);
    }
    cluster[r] = home;
  }
  return { cluster, representative };
}
