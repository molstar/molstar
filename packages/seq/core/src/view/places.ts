/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// From a system and a position to its places in a section of blocks, and back (_spec/README.md "Resolved view",
// Places): a block's units are its positions times its scale from its first display unit, each block in its own
// direction; an index by system of the blocks showing it, sorted by their first shown position, built once per resolved
// section, so `placesOf` is a binary search and a scan of the blocks over the position (overlapping blocks: every one).

import { axisTotal, elementAt, offsetOf, pieceOffsets, positionAt } from './axis.js';
import { blockUnits } from './blocks.js';
import type { ResolvedSection } from './view.js';

type Block = NonNullable<ResolvedSection['blocks']>[number];

/** A system and a position on it: what a viewport's end, a history step and a saved view carry over blocks. */
export interface SystemPosition {
  readonly system: string;
  readonly position: number;
}

/** A place of a position in a section: the block showing it and the section's display unit there. */
export interface BlockPlace {
  readonly block: number;
  readonly unit: number;
}

/** A block's display units a position: its `scale` (absent: 1); a `whole` chip's one unit shared by its positions. */
export function blockScaleOf(block: Pick<Block, 'scale' | 'form' | 'axis'>): number {
  if (block.form === 'whole' && block.axis.ranges) return 1 / Math.max(1, axisTotal(block.axis) ?? 1);
  return block.scale ?? 1;
}

/** The section unit of a continuous position on a block, or undefined when the block does not show it. */
export function blockUnitOf(block: Block, position: number): number | undefined {
  const o = offsetOf(block.axis, position);
  if (o === undefined) return undefined;
  return (block.offset ?? 0) + o * blockScaleOf(block);
}

/** The continuous position at a section unit on a block (clamped to the block), or undefined while it is pending. */
export function blockPositionAt(block: Block, unit: number, side: 'start' | 'end' = 'end'): number | undefined {
  return positionAt(block.axis, (unit - (block.offset ?? 0)) / blockScaleOf(block), side);
}

/** The element (a whole position) at a section unit on a block, or undefined while it is pending. */
export function blockElementAt(block: Block, unit: number): number | undefined {
  return elementAt(block.axis, (unit - (block.offset ?? 0)) / blockScaleOf(block));
}

/** A block's display units in its section (its positions times its scale; a chip one), 0 while pending. */
export function blockSpan(block: Block): number {
  const u = blockUnits(block);
  return u === undefined ? 0 : block.form === 'whole' ? u : u * (block.scale ?? 1);
}

/** The index of a section's blocks by system: each system's blocks sorted by their least shown position. */
interface PlacesIndex {
  readonly bySystem: ReadonlyMap<string, { readonly blocks: Int32Array; readonly lo: Float64Array }>;
}

const Indexes = new WeakMap<ResolvedSection, PlacesIndex>();

function indexOf(section: ResolvedSection): PlacesIndex {
  let ix = Indexes.get(section);
  if (ix) return ix;
  const lists = new Map<string, { b: number; lo: number }[]>();
  const blocks = section.blocks ?? [];
  for (let b = 0; b < blocks.length; b++) {
    const r = blocks[b].axis.ranges;
    if (!r || blocks[b].offset === undefined) continue;
    let lo = Infinity;
    for (let k = 0; k < r.length; k += 2) lo = Math.min(lo, r[k]);
    for (const system of new Set(blocks[b].systems.length ? blocks[b].systems : [blocks[b].axis.system])) {
      let l = lists.get(system);
      if (!l) lists.set(system, (l = []));
      l.push({ b, lo });
    }
  }
  const bySystem = new Map<string, { blocks: Int32Array; lo: Float64Array }>();
  for (const [system, l] of lists) {
    l.sort((x, y) => x.lo - y.lo || x.b - y.b);
    bySystem.set(system, { blocks: Int32Array.from(l, (e) => e.b), lo: Float64Array.from(l, (e) => e.lo) });
  }
  ix = { bySystem };
  Indexes.set(section, ix);
  return ix;
}

/**
 * Every place showing `position` on `system` in a section of blocks (_spec/README.md "Resolved view", Places): each
 * block that shows it (overlapping blocks: every one, the merged-copies rule) and the section unit there, in block
 * order; none when no block shows it. Into `out` when given (cleared first). A section without blocks has its own axis:
 * one place on block −1 when the axis is on the system and shows the position.
 */
export function placesOf(
  section: ResolvedSection,
  system: string,
  position: number,
  out: BlockPlace[] = [],
): BlockPlace[] {
  out.length = 0;
  if (!section.blocks) {
    if (section.axis.system !== system) return out;
    const o = offsetOf(section.axis, position);
    if (o !== undefined) out.push({ block: -1, unit: o });
    return out;
  }
  const list = indexOf(section).bySystem.get(system);
  if (!list) return out;
  // the blocks whose least position is at most `position`: a binary search, then each of them asked
  let lo = 0;
  let hi = list.lo.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (list.lo[mid] <= position) lo = mid + 1;
    else hi = mid;
  }
  for (let k = 0; k < lo; k++) {
    const b = list.blocks[k];
    const unit = blockUnitOf(section.blocks[b], position);
    if (unit !== undefined) out.push({ block: b, unit });
  }
  out.sort((x, y) => x.block - y.block);
  return out;
}

/** The block holding a section unit (the last starting at or before it; a gap's: the one before), or −1. */
export function blockAtUnit(section: ResolvedSection, unit: number): number {
  const blocks = section.blocks;
  if (!blocks?.length) return -1;
  let lo = 0;
  let hi = blocks.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if ((blocks[mid].offset ?? 0) <= unit) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/**
 * A section unit as a system and a position (_spec/README.md "Resolved view", Places): the block holding it, the
 * continuous position there (`side` on a block's edge: `start` the next block's first, `end` the one before's last), or
 * undefined while pending. What a viewport's end, a history step and a saved view carry: never units.
 */
export function placeOfUnit(
  section: ResolvedSection,
  unit: number,
  side: 'start' | 'end' = 'end',
): SystemPosition | undefined {
  const blocks = section.blocks;
  if (!blocks) {
    const p = positionAt(section.axis, unit, side);
    return p === undefined ? undefined : { system: section.axis.system, position: p };
  }
  let b = blockAtUnit(section, unit);
  if (b < 0) return undefined;
  // on an edge as a start: the next block's first
  if (side === 'start' && b + 1 < blocks.length && unit >= (blocks[b].offset ?? 0) + blockSpan(blocks[b])) b++;
  // on a block's first unit as an end: the one before's last
  if (side === 'end' && b > 0 && unit <= (blocks[b].offset ?? 0)) b--;
  const block = blocks[b];
  const p = blockPositionAt(block, unit, side);
  return p === undefined ? undefined : { system: block.axis.system, position: p };
}

/**
 * The section unit of a system and a position (`placeOfUnit`'s inverse): the first block showing it, else, a position
 * no block shows (a block no longer shown, a changed `include`), the nearest edge of a block on that system, else of
 * the nearest block in order (`near`, a unit) or the section's start; undefined while pending.
 */
export function unitOfPlace(section: ResolvedSection, at: SystemPosition, near?: number): number | undefined {
  const shown = placesOf(section, at.system, at.position);
  if (shown.length) {
    if (near === undefined || shown.length === 1) return shown[0].unit;
    let best = shown[0].unit;
    for (const p of shown) if (Math.abs(p.unit - near) < Math.abs(best - near)) best = p.unit;
    return best;
  }
  const blocks = section.blocks;
  if (!blocks) return undefined;
  // the nearest shown edge on the same system: clamped to what is drawn
  let best: number | undefined;
  let gap = Infinity;
  for (const block of blocks) {
    const r = block.axis.ranges;
    if (block.axis.system !== at.system || !r || block.offset === undefined) continue;
    const offsets = pieceOffsets(r);
    for (let k = 0; k < r.length; k += 2) {
      for (const edge of [r[k], r[k + 1]]) {
        const d = Math.abs(edge - at.position);
        if (d >= gap) continue;
        gap = d;
        const o = block.axis.dir === -1 ? offsets[k >> 1] + r[k + 1] - edge : offsets[k >> 1] + edge - r[k];
        best = block.offset + o * blockScaleOf(block);
      }
    }
  }
  if (best !== undefined) return best;
  if (!blocks.length || blocks[0].offset === undefined) return undefined;
  return near === undefined ? 0 : near;
}

/**
 * Whether two sections show equal block lists: the same systems, ranges, directions and scales in order (_spec/README.md
 * "Resolved view", Places). Such sections follow by identity in units (two sections of one assembly, a ring of the
 * genome and its linear detail); unequal lists do not follow (`seq.render.follow.blocks`).
 */
export function sameBlocks(a: ResolvedSection, b: ResolvedSection): boolean {
  const x = a.blocks;
  const y = b.blocks;
  if (!x || !y || x.length !== y.length) return false;
  for (let i = 0; i < x.length; i++) {
    const p = x[i];
    const q = y[i];
    if (p.axis.system !== q.axis.system || p.axis.dir !== q.axis.dir || blockScaleOf(p) !== blockScaleOf(q))
      return false;
    if (p.offset !== q.offset || p.form !== q.form) return false;
    const r = p.axis.ranges;
    const s = q.axis.ranges;
    if (r === s) continue;
    if (!r || !s || r.length !== s.length) return false;
    for (let k = 0; k < r.length; k++) if (r[k] !== s[k]) return false;
  }
  return true;
}
