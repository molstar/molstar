/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { blockSpan, type ResolvedBlock } from '@molstar/seqx-core/view';

// Switching chromosomes (_spec/README.md "Behaviors", A section of blocks, Switching blocks): what window a jump to
// another block of a row lands on, one rule for a click on the navigator, `block.next` and `block.previous`, and the
// square's click (each axis by its own section).

/** The block of a row holding a unit (the last starting at or before it). */
export function blockOfUnit(blocks: readonly ResolvedBlock[], u: number): number {
  let lo = 0;
  let hi = blocks.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if ((blocks[mid].offset ?? 0) <= u) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/**
 * The window a switch from the view's block to block `to` lands on (IGV's: a picked chromosome whole): from a view
 * showing its block whole (its ends within `tolerance` units, a pixel's) the new block whole; zoomed inside a block, the
 * same width centred at the same relative place in the new block, clamped to its ends; a width past the new block's
 * length, the new block whole. The view's block is the one holding its centre.
 */
export function blockSwitchRange(
  blocks: readonly ResolvedBlock[],
  view: { readonly start: number; readonly end: number },
  to: number,
  tolerance = 0,
): readonly [number, number] {
  const t = blocks[to];
  const lo = t.offset ?? 0;
  const len = blockSpan(t);
  const whole = [lo, lo + len] as const;
  const c = (view.start + view.end) / 2;
  const from = blocks[blockOfUnit(blocks, c)];
  const fo = from.offset ?? 0;
  const fl = blockSpan(from);
  const w = view.end - view.start;
  const tol = Math.max(tolerance, 1e-9);
  if (Math.abs(view.start - fo) <= tol && Math.abs(view.end - (fo + fl)) <= tol) return whole;
  if (!(fl > 0) || w >= len) return whole;
  const rel = Math.min(1, Math.max(0, (c - fo) / fl));
  const start = Math.min(Math.max(lo + rel * len - w / 2, lo), lo + len - w);
  return [start, start + w];
}
