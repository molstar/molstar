/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's integer geometry and kind masks, shared by its draw, `tiny` and hit
// (seqspec:viewer#extension-points): rows and stripes in whole device pixels, ints in and out, so no double crosses a
// call and the hit finds the pixels the draw filled. The rows split the entry a device pixel inside its top and bottom
// into equal rows two device pixels apart where a row is 6 or more (the last ending where the track's bottom border
// begins); the kinds present split each row into equal stripes, the first on top. A row with fewer than 2 device pixels
// a kind is compact.

import type { AxisView } from "@molstar/seq-viewer/draw";
import { Kind, Kinds } from "../vocabulary.js";

/** The stripes' order: these kinds first, then the others in the vocabulary's (the legend's) order. */
const FIRST = [Kind.hbond, Kind.pi_stacking, Kind.cation_pi, Kind.hydrophobic];
export const STRIPES: readonly number[] = [...FIRST, ...Kinds.map((_, k) => k).filter((k) => !FIRST.includes(k))];

/** Per mask, its first kind in stripe order (255: none): a compact row's color and `part`. Built once. */
export const LEAD = new Uint8Array(1 << Kinds.length).fill(255);
for (let m = 1; m < LEAD.length; m++) LEAD[m] = STRIPES.find((k) => m & (1 << k)) ?? 255;

/** The frame `rowFrame` fills, read by the draw and the hit: the rows' top and height, the gap between rows, the rows. */
export const geo = new Int32Array(4);
const A = 0;
const H = 1;
const G = 2;
const S = 3;
export const Geo = { a: A, height: H, gap: G, rows: S } as const;

/** The rows' box of the entry `v` (a device pixel inside its top and bottom) for `rows` rows, into `geo`. */
export function rowFrame(v: AxisView, rows: number): void {
  const a = Math.round((v.top + 1) * v.dpr) | 0;
  const h = (Math.round((v.top + v.height - 1) * v.dpr) | 0) - a;
  geo[A] = a;
  geo[H] = h;
  geo[G] = h >= 6 * rows ? 2 : 0;
  geo[S] = rows;
}

/** Row q's top edge (q = S: the end of the last gap); `Math.round`'s half up, in ints. */
export function rowTop(a: number, h: number, g: number, rows: number, q: number): number {
  return a + (((2 * q * (h + g) + rows) / (2 * rows)) | 0);
}

/** Stripe j's top edge in a row from `y0` of height `h` split into `k` stripes (j = k: the row's end). */
export function stripeTop(y0: number, h: number, k: number, j: number): number {
  return y0 + (((2 * j * h + k) / (2 * k)) | 0);
}

/** Whether `kinds` stripes do not fit the frame's rows (under 2 device pixels a kind): one decision for the layer. */
export function compactRows(kinds: number): boolean {
  return (((geo[H] + geo[G]) / geo[S]) | 0) - geo[G] < 2 * kinds;
}
