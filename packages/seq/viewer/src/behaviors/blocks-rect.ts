/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Location, location } from '@molstar/seqx-core/query';
import { type Axis, axisTotal, cellOffset, elementAt } from '@molstar/seqx-core/view';

// The rect of a drag on a grid over blocks (_spec/README.md "Behaviors", A rect over blocks): internal to the select
// behavior, not re-exported by the behaviors' entry (its tests import this module).

/** A block a grid over blocks is laid on: its system, its axis and its first display unit (internal: the rect drag's). */
export interface RectBlock {
  readonly system: string;
  readonly axis: Axis;
  readonly offset: number;
}

/** A cell of a grid over blocks: its x and y blocks (their indices in the section) and its positions on them. */
export interface BlocksCell {
  readonly block: number;
  readonly position: number;
  readonly yBlock: number;
  readonly yPosition: number;
}

/**
 * The rect from one cell to another over a section's blocks (a grid over a complex):
 * each cell's display units on its own blocks (the hit strip's, so a copy split is told from its twin), the span between
 * them on each axis split at the blocks' boundaries, a part a block, a rect a pair of parts (copies of one system merge).
 * Undefined when a cell is on no block.
 */
export function blocksRect(blocks: readonly RectBlock[], a: BlocksCell, b: BlocksCell): Location | undefined {
  const unit = (k: number, position: number) => {
    const block = blocks[k];
    const o = block ? cellOffset(block.axis, Math.floor(position)) : undefined;
    return o === undefined ? undefined : block.offset + o;
  };
  const ax = unit(a.block, a.position);
  const bx = unit(b.block, b.position);
  const ay = unit(a.yBlock, a.yPosition);
  const by = unit(b.yBlock, b.yPosition);
  if (ax === undefined || bx === undefined || ay === undefined || by === undefined) return undefined;
  // each block's part of the units [u0, u1]: its positions there (a reversed block's from its last unit)
  const parts = (u0: number, u1: number) => {
    const out: { system: string; ranges: number[] }[] = [];
    for (const k of blocks) {
      const n = axisTotal(k.axis) ?? 0;
      const s = Math.max(u0, k.offset);
      const e = Math.min(u1, k.offset + n - 1);
      if (e < s) continue;
      const p0 = elementAt(k.axis, s - k.offset);
      const p1 = elementAt(k.axis, e - k.offset);
      if (p0 === undefined || p1 === undefined) continue;
      out.push({ system: k.system, ranges: [Math.min(p0, p1), Math.max(p0, p1) + 1] });
    }
    return out;
  };
  const xs = parts(Math.min(ax, bx), Math.max(ax, bx));
  const ys = parts(Math.min(ay, by), Math.max(ay, by));
  return location({ rects: xs.flatMap((x) => ys.map((y) => ({ a: x, b: y }))) });
}
