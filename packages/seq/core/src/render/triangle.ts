/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid's or a dot plot's `triangle` on a strip (_spec/README.md "Render", The grid layer: triangles): it halves the
// plot of one system against itself, along the diagonal as drawn. Why not on raw positions everywhere: over two systems
// (chr1 × chr4, a block pair off the diagonal, a grid over two sections) positions of one axis are no positions of the
// other, and comparing them cut a pair in half.

import type { LayoutContext } from './layer-type.js';

/** A cell kept by a triangle (`upper`: y before x, `lower`: y after x), by the two display ranks. */
export type Triangle = 'both' | 'upper' | 'lower';

/** Whether ranks (i, j) are kept by the triangle (`upper`: y before x, `lower`: y after x). */
export function keeps(t: Triangle, i: number, j: number): boolean {
  return t === 'both' || (t === 'upper' ? j < i : j > i);
}

/** How a triangle keeps a strip: every cell, none, or each cell by its display ranks. */
export const StripKeep = { All: 0, None: 1, Cells: 2 } as const;
export type StripKeep = (typeof StripKeep)[keyof typeof StripKeep];

/**
 * A triangle on a strip: over two sections of blocks every cell (not one assembly against itself); over one section of
 * blocks a pair of two blocks whole on the side the triangle keeps, nothing on the other, a block against itself by its
 * cells; on two systems every cell; else by cells (a categorical axis's by display rank).
 */
export function stripKeep(t: Triangle, ctx: LayoutContext): StripKeep {
  if (t === 'both' || ctx.section.y?.blocks) return StripKeep.All;
  const st = ctx.section.blocks && ctx.place?.strips;
  if (st && ctx.strip !== undefined) {
    const bx = st.block[ctx.strip];
    const by = st.yBlock[ctx.strip];
    if (bx !== by) return keeps(t, bx, by) ? StripKeep.All : StripKeep.None;
  }
  if (ctx.axisY && ctx.axis.system !== ctx.axisY.system) return StripKeep.All;
  return StripKeep.Cells;
}
