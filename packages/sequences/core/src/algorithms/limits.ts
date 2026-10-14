/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The size bound of the aligners that keep a traceback for the full matrix (seqspec:format#mappings): needleman_wunsch
// and biomsa's complete DP. `myers` keeps linear memory and has no bound.

/** The most DP cells, (n + 1) × (m + 1), a full-matrix traceback may cover: 2^26, about 32 MB for biomsa's half byte a cell. */
export const MAX_CELLS = 2 ** 26;

/** Thrown before anything is allocated when two sequences need more than `MAX_CELLS` cells. Not an options error. */
export class AlignerSizeError extends Error {
  readonly cells: number;
  constructor(aligner: string, n: number, m: number) {
    const cells = (n + 1) * (m + 1);
    super(
      `${aligner}: aligning ${n} against ${m} letters takes ${cells} DP cells; the limit is ${MAX_CELLS}. ` +
        "The myers aligner has no such limit.",
    );
    this.name = "AlignerSizeError";
    this.cells = cells;
  }
}

/** Throws `AlignerSizeError` when aligning `n` against `m` letters would pass `MAX_CELLS`. */
export function checkCells(aligner: string, n: number, m: number): void {
  if ((n + 1) * (m + 1) > MAX_CELLS) throw new AlignerSizeError(aligner, n, m);
}
