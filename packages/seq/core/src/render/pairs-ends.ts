/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Where a pairs layer's ends are in a section (_spec/README.md "Render", The pairs layer, Ends on blocks): each end's
// centre as the section's display units, on the section's axis or on a block. An end is on the first block showing
// its system: the block of its system or of its row's unit, the block a merged copy of it is merged into, or a block of
// a unit of it (a sequence or a molecule reaches every unit over it, EC-48). A pair whose two ends are on one block of
// one system is drawn again on each copy of that block (the resolver's copy groups, `ResolvedBlock.copy`): split copies
// each show it, a unit named by the data shows it alone.

import type { PairsData, UnitEntry } from '@molstar/seqx-core/resolve';
import { cellOffset, type ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';

/** The ends drawn in a section: a pair's two ends' centres as display units (NaN: not on it), and its pair, by entry. */
export interface PairEnds {
  /** The axis pieces they were placed on (kept by them: a source's coverage fills the layout's). */
  readonly pieces: readonly number[] | undefined;
  /** Entries: a pair once, and once again on each copy it is drawn on. */
  readonly n: number;
  /** Two a entry: `a`'s centre, `b`'s. */
  readonly u: Float64Array;
  /** Each entry's pair. */
  readonly pair: Int32Array;
}

/** What placing a pair's ends reads. */
export interface PairEndsInput {
  readonly data: PairsData;
  readonly aSystem: string;
  readonly bSystem: string;
  /** A complex's links: the units their ends are on (`PairsData.aUnit`, `bUnit` index them). */
  readonly endUnits?: readonly string[];
  /** The model's units: a block of a unit shows what it instantiates. */
  readonly molecularUnits: ReadonlyMap<string, UnitEntry>;
}

type Axis = ResolvedSection['axis'];

/** The centre of `[start, end)` on an axis as a display unit, NaN when not shown. */
function centreOf(axis: Parameters<typeof cellOffset>[0], start: number, end: number): number {
  const a = cellOffset(axis, start);
  const b = cellOffset(axis, end - 1);
  if (a === undefined || b === undefined) return Number.NaN;
  return (Math.min(a, b) + Math.max(a, b) + 1) / 2;
}

/**
 * What each block shows, flat: for block `i`, entries `reach[i]` to `reach[i + 1]` of `systems`, each a system with
 * the shift and length of its positions on the block (a unit of a molecule: its start and length; its own: 0, ∞).
 */
interface Reach {
  readonly at: Int32Array;
  readonly systems: string[];
  readonly shift: number[];
  readonly length: number[];
}

function reachOf(section: ResolvedSection, units: ReadonlyMap<string, UnitEntry>): Reach {
  const blocks = section.blocks!;
  const at = new Int32Array(blocks.length + 1);
  const systems: string[] = [];
  const shift: number[] = [];
  const length: number[] = [];
  let from = 0;
  const add = (system: string, s: number, n: number) => {
    for (let j = from; j < systems.length; j++) if (systems[j] === system) return;
    systems.push(system);
    shift.push(s);
    length.push(n);
  };
  for (let i = 0; i < blocks.length; i++) {
    at[i] = from = systems.length;
    const b = blocks[i];
    add(b.axis.system, 0, Number.POSITIVE_INFINITY);
    for (const s of b.systems) add(s, 0, Number.POSITIVE_INFINITY);
    for (const s of b.systems) {
      const u = units.get(s);
      if (u) add(u.of, u.start, u.end === undefined ? Number.POSITIVE_INFINITY : u.end - u.start);
    }
  }
  at[blocks.length] = systems.length;
  return { at, systems, shift, length };
}

/** `[s, e)` on `system` on block `i` as the section's display units, NaN when the block does not show it. */
function onBlock(
  section: ResolvedSection,
  reach: Reach,
  shape: (axis: Axis) => Axis,
  i: number,
  system: string,
  s: number,
  e: number,
): number {
  const b = section.blocks![i];
  if (b.offset === undefined) return Number.NaN;
  for (let j = reach.at[i]; j < reach.at[i + 1]; j++) {
    if (reach.systems[j] !== system) continue;
    const from = s - reach.shift[j];
    const to = e - reach.shift[j];
    if (from < 0 || to > reach.length[j]) return Number.NaN;
    // a chip is one unit: every residue of its molecule's ends at its middle
    if (b.form === 'whole') return b.offset + 0.5;
    const c = centreOf(shape(b.axis), from, to);
    return Number.isNaN(c) ? c : b.offset + c;
  }
  return Number.NaN;
}

/** Placing: the first block holding an end, into `At`. */
const At = { block: -1, u: Number.NaN };

function firstOn(
  section: ResolvedSection,
  reach: Reach,
  shape: (axis: Axis) => Axis,
  system: string | undefined,
  s: number,
  e: number,
) {
  At.block = -1;
  At.u = Number.NaN;
  if (system === undefined) return;
  const n = section.blocks!.length;
  for (let i = 0; i < n; i++) {
    const u = onBlock(section, reach, shape, i, system, s, e);
    if (!Number.isNaN(u)) {
      At.block = i;
      At.u = u;
      return;
    }
  }
}

/**
 * Each pair's ends in a section (the module's comment): on the section's axis, or on its blocks, a pair whose ends are
 * on one block of one system again on each of that block's copies.
 */
export function pairEnds(input: PairEndsInput, section: ResolvedSection, layout: AxisTransform): PairEnds {
  // the pieces the layout placed the axis on: the resolved axis has none while its letters come from a source
  const pieces = section.axis.ranges ?? layout.pieces;
  const d = input.data;
  const shape = (axis: Axis): Axis =>
    axis.ranges || axis.system !== layout.system ? axis : { ...axis, ranges: layout.pieces };
  const units = input.endUnits;
  const aOf = (k: number) => (units && d.aUnit ? units[d.aUnit[k]] : input.aSystem);
  const bOf = (k: number) => (units && d.bUnit ? units[d.bUnit[k]] : input.bSystem);
  const blocks = section.blocks;
  if (!blocks) {
    const u = new Float64Array(2 * d.n);
    const pair = new Int32Array(d.n);
    const on = (system: string | undefined, s: number, e: number) =>
      section.axis.system === system ? centreOf(shape(section.axis), s, e) : Number.NaN;
    for (let k = 0; k < d.n; k++) {
      u[2 * k] = on(aOf(k), d.aStart[k], d.aEnd[k]);
      u[2 * k + 1] = on(bOf(k), d.bStart[k], d.bEnd[k]);
      pair[k] = k;
    }
    return { pieces, n: d.n, u, pair };
  }
  const reach = reachOf(section, input.molecularUnits);
  const u: number[] = [];
  const pair: number[] = [];
  for (let k = 0; k < d.n; k++) {
    const as = aOf(k);
    const bs = bOf(k);
    firstOn(section, reach, shape, as, d.aStart[k], d.aEnd[k]);
    const ia = At.block;
    u.push(At.u);
    firstOn(section, reach, shape, bs, d.bStart[k], d.bEnd[k]);
    const ib = At.block;
    u.push(At.u);
    pair.push(k);
    // one system, one block (or one end not placed): each copy of that block holding the end shows the pair too
    const i = ia >= 0 ? ia : ib;
    if (i < 0 || as !== bs || (ia >= 0 && ib >= 0 && ia !== ib)) continue;
    const copy = blocks[i].copy;
    for (let j = i + 1; j < blocks.length; j++) {
      if (blocks[j].copy !== copy) continue;
      const ua = ia < 0 ? Number.NaN : onBlock(section, reach, shape, j, as!, d.aStart[k], d.aEnd[k]);
      const ub = ib < 0 ? Number.NaN : onBlock(section, reach, shape, j, bs!, d.bStart[k], d.bEnd[k]);
      if (Number.isNaN(ua) && Number.isNaN(ub)) continue;
      u.push(ua, ub);
      pair.push(k);
    }
  }
  return { pieces, n: pair.length, u: Float64Array.from(u), pair: Int32Array.from(pair) };
}
