/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Zoomed-out reductions by the theme's kind (seqspec:sources#levels-and-the-zoom-gate;
// seqspec:viewer#extension-points): the indexes a draw reduces device columns from, built once per prepared state and
// theme (a class or a slot is a theme color), never per frame. A `classes` theme's cells reduce by the plurality of
// their (color, alpha) classes, a block's own (no class is shared across blocks, none past 255 is folded into another:
// such a block draws through the raw per-column path); a `mean` theme's by the mean value (the prepared values'
// sum-and-count pyramid, prepare.ts); items by their counts per color slot, ranges as coverage and points once at their
// center, never over a rows track's union.

import type { Features } from "./prepare.ts";
import { type ClassPyramid, type CountPyramid, classPyramidFrom, countPyramid } from "./pyramids.ts";

/** What the class build reads per cell of a block (the draw's theme at a block index; RN1 fills it). */
export interface CellColoring {
  /** 0xRRGGBB, or -1 (none). */
  color(block: number, index: number): number;
  /** The theme's alpha as a byte, 0..255 (its own byte beside the color). */
  alpha(block: number, index: number): number;
  /** Equal keys ≥ 0 have equal colors and alphas (one theme call per run); -1: none. */
  key(block: number, index: number): number;
}
/** A block's cells as classes 1..255 (0: none): their pyramid, and per class its color and alpha byte. */
export interface CellClasses {
  readonly pyramid: ClassPyramid;
  readonly color: Int32Array;
  readonly alpha: Uint8Array;
  /** Classes used. */
  readonly count: number;
}

const job = {
  c: null as CellColoring | null,
  b: 0,
  base: 0,
  color: new Int32Array(256),
  alpha: new Uint8Array(256),
  count: 0,
  pairs: new Map<number, number>(),
};
/** Fills a chunk of classes; false past 255 classes. */
function fillClasses(start: number, out: Uint8Array): boolean {
  const c = job.c as CellColoring;
  let lastKey = -1;
  let lastPair = -2;
  let last = 0;
  for (let k = 0; k < out.length; k++) {
    const i = job.base + start + k;
    const key = c.key(job.b, i);
    if (key < 0 || key !== lastKey) {
      const color = c.color(job.b, i);
      const pair = color < 0 ? -1 : color * 256 + (c.alpha(job.b, i) & 255);
      let cls = pair === lastPair ? last : pair < 0 ? 0 : (job.pairs.get(pair) ?? -1);
      if (cls < 0) {
        if (job.count === 255) return false;
        cls = ++job.count;
        job.pairs.set(pair, cls);
        job.color[cls] = color;
        job.alpha[cls] = pair & 255;
      }
      last = cls;
      lastKey = key;
      lastPair = pair;
    }
    out[k] = last;
  }
  return true;
}
/**
 * Block b's cell classes over its indices [base, base + length) under a coloring, or null past 255 distinct (color,
 * alpha) pairs (such a block draws through the raw per-column path). Edge-rate: once per prepared state, theme and
 * block, on the first aggregate draw.
 */
export function cellClasses(c: CellColoring, b: number, base: number, length: number): CellClasses | null {
  Object.assign(job, { c, b, base, count: 0, color: new Int32Array(256).fill(-1), alpha: new Uint8Array(256) });
  job.pairs.clear();
  const pyramid = classPyramidFrom(length, fillClasses);
  job.c = null;
  return pyramid && { pyramid, color: job.color, alpha: job.alpha, count: job.count };
}

/** Color slots for item density: a layer's first 8 colors; later ones share the last. */
export const SLOTS = 8;
/** Bins of 2^6 positions at the finest: below, the items in view are few enough to count one by one. */
export const ITEM_BASE = 6;
/**
 * A feature layer's counts per color slot (stacked or majority density), per block on first use:
 * with `point` each item once at its center, summed per column (points, variants); else as coverage, a column taking
 * its bins' maximum (ranges).
 */
export class ItemDensity {
  /** Per feature its slot; per slot its color, in first-seen (palette) order. */
  readonly slot: Uint8Array;
  readonly color = new Int32Array(SLOTS).fill(-1);
  readonly slots: number;
  readonly point: boolean;
  private readonly f: Features;
  private readonly blocks: (CountPyramid | null | undefined)[] = [];
  constructor(f: Features, colorOf: (k: number) => number, point: boolean) {
    this.f = f;
    this.point = point;
    this.slot = new Uint8Array(f.count);
    let slots = 0;
    for (let k = 0; k < f.count; k++) {
      const c = colorOf(k);
      let s = 0;
      while (s < slots && this.color[s] !== c) s++;
      if (s === slots && slots < SLOTS) this.color[slots++] = c;
      this.slot[k] = Math.min(s, SLOTS - 1);
    }
    this.slots = Math.max(1, slots);
  }
  /** Block b's count pyramid over its indices [base, base + length); null on a rows track. */
  pyramid(b: number, base: number, length: number): CountPyramid | null {
    const f = this.f;
    if (f.hasRows || b * f.lanes >= f.seg.length - 1) return null;
    let p = this.blocks[b];
    if (p === undefined) {
      const items = { start: f.start, stop: f.stop, slot: this.slot, offset: base, point: this.point };
      p = countPyramid(length, ITEM_BASE, this.slots, {
        ...items,
        from: f.seg[b * f.lanes],
        to: f.seg[(b + 1) * f.lanes],
      });
      this.blocks[b] = p;
    }
    return p;
  }
}

/**
 * One mode per layer per frame: density past the type's threshold of items in view, items again
 * only below 0.8 of it (a resize near the threshold does not flip it frame by frame).
 */
export const densityOn = (count: number, threshold: number, was: boolean): boolean =>
  was ? count >= 0.8 * threshold : count > threshold;
