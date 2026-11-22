/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Where a pairs layer's ends are in a section (_spec/README.md "Render", The pairs layer, Ends on blocks): each end's
// centre as the section's display units, on the section's axis or on a block. An end is on the first block showing
// its system: the block of its system or of its row's unit, the block a merged copy of it is merged into, or a block of
// a unit of it (a sequence or a molecule reaches every unit over it, EC-48, by the resolver's step down to the unit).
// A pair whose two ends are on one block of one system is drawn again on each copy of that block (the resolver's copy
// groups, `ResolvedBlock.copy`). A pair between two systems with an end on split copies its data names no unit of is
// on no copy (never a guessed one): flagged, said by the layer (`seq.render.pairs.copy`).

import { composeEdges, type Path, projectRanges, type PairsData, type StaticModel } from '@molstar/seqx-core/resolve';
import { cellOffset, type ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';

/** The ends drawn in a section: each entry's two ends' centres as display units (NaN: not on it), and its pair. */
export interface PairEnds {
  /** The axis pieces they were placed on (kept by them: a source's coverage fills the layout's). */
  readonly pieces: readonly number[] | undefined;
  /** Entries: a pair once, and once again on each copy it is drawn on. */
  readonly n: number;
  /** Two a entry: `a`'s centre, `b`'s. */
  readonly u: Float64Array;
  /** Each entry's pair. */
  readonly pair: Int32Array;
  /**
   * Each entry's ends on split copies its data names no unit of (bit 1: `a`, 2: `b`): drawn on no copy; bit 4, its
   * other end is on no block either.
   */
  readonly unnamed: Uint8Array;
  /** The copies shown of `a`'s system and of `b`'s where an entry is `unnamed` there (0: none). */
  readonly aCopies: number;
  readonly bCopies: number;
}

/** What placing a pair's ends reads. */
export interface PairEndsInput {
  readonly data: PairsData;
  readonly aSystem: string;
  readonly bSystem: string;
  /** A complex's links: the units their ends are on (`PairsData.aUnit`, `bUnit` index them). */
  readonly endUnits?: readonly string[];
  /** The resolver's steps: a unit's step down from what it instantiates. */
  readonly model: Pick<StaticModel, 'edgesFrom'>;
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
 * What each block shows, flat: for block `i`, entries `at[i]` to `at[i + 1]`, each a system and the resolver's path
 * from it onto the block's (undefined: the block's own, or a merged copy's).
 */
interface Reach {
  readonly at: Int32Array;
  readonly systems: string[];
  readonly paths: (Path | undefined)[];
}

function reachOf(section: ResolvedSection, model: PairEndsInput['model']): Reach {
  const blocks = section.blocks!;
  const at = new Int32Array(blocks.length + 1);
  const systems: string[] = [];
  const paths: (Path | undefined)[] = [];
  let from = 0;
  const add = (system: string, path: Path | undefined) => {
    for (let j = from; j < systems.length; j++) if (systems[j] === system) return;
    systems.push(system);
    paths.push(path);
  };
  for (let i = 0; i < blocks.length; i++) {
    at[i] = from = systems.length;
    const b = blocks[i];
    add(b.axis.system, undefined);
    for (const s of b.systems) add(s, undefined);
    // a unit's molecule (or sequence) by the unit's step down from it
    for (const s of b.systems)
      for (const up of model.edgesFrom(s)) {
        if (up.unit !== 'up') continue;
        const down = model.edgesFrom(up.to).find((e) => e.unit === 'down' && e.to === s);
        const path = down && composeEdges(up.to, [down]);
        if (path) add(up.to, path);
      }
  }
  at[blocks.length] = systems.length;
  return { at, systems, paths };
}

const Projected: number[] = [];

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
    const path = reach.paths[j];
    let from = s;
    let to = e;
    if (path) {
      // onto the unit: the whole end within its range, else not on it
      projectRanges(path, [s, e], Projected);
      if (Projected.length !== 2 || Projected[1] - Projected[0] !== e - s) return Number.NaN;
      from = Projected[0];
      to = Projected[1];
    }
    // a chip is one unit: every residue of its molecule's ends at its middle
    if (b.form === 'whole') return b.offset + 0.5;
    // a block's units are its positions times its scale (`ResolvedBlock.scale`)
    const c = centreOf(shape(b.axis), from, to);
    return Number.isNaN(c) ? c : b.offset + c * (b.scale ?? 1);
  }
  return Number.NaN;
}

/** Placing an end: its first block, its centre there, its copies holding it (blocks of that block's copy group). */
const At = { block: -1, u: Number.NaN, copies: 0 };

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
  At.copies = 0;
  if (system === undefined) return;
  const blocks = section.blocks!;
  for (let i = 0; i < blocks.length; i++) {
    const u = onBlock(section, reach, shape, i, system, s, e);
    if (Number.isNaN(u)) continue;
    At.block = i;
    At.u = u;
    At.copies = 1;
    for (let j = i + 1; j < blocks.length; j++)
      if (blocks[j].copy === blocks[i].copy && !Number.isNaN(onBlock(section, reach, shape, j, system, s, e)))
        At.copies++;
    return;
  }
}

/** Each pair's ends in a section (the module's comment). */
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
    return { pieces, n: d.n, u, pair, unnamed: new Uint8Array(d.n), aCopies: 0, bCopies: 0 };
  }
  const reach = reachOf(section, input.model);
  const u: number[] = [];
  const pair: number[] = [];
  const unnamed: number[] = [];
  let aCopies = 0;
  let bCopies = 0;
  for (let k = 0; k < d.n; k++) {
    const as = aOf(k);
    const bs = bOf(k);
    firstOn(section, reach, shape, as, d.aStart[k], d.aEnd[k]);
    const ia = At.block;
    const ua = At.u;
    const ca = At.copies;
    firstOn(section, reach, shape, bs, d.bStart[k], d.bEnd[k]);
    const ib = At.block;
    const ub = At.u;
    const cb = At.copies;
    pair.push(k);
    // one system, one block (or one end on none): the pair is on each copy of that block
    const i = ia >= 0 ? ia : ib;
    if (i >= 0 && as === bs && (ia < 0 || ib < 0 || ia === ib)) {
      u.push(ua, ub);
      unnamed.push(0);
      const copy = blocks[i].copy;
      for (let j = i + 1; j < blocks.length; j++) {
        if (blocks[j].copy !== copy) continue;
        const ja = ia < 0 ? Number.NaN : onBlock(section, reach, shape, j, as!, d.aStart[k], d.aEnd[k]);
        const jb = ib < 0 ? Number.NaN : onBlock(section, reach, shape, j, bs!, d.bStart[k], d.bEnd[k]);
        if (Number.isNaN(ja) && Number.isNaN(jb)) continue;
        u.push(ja, jb);
        pair.push(k);
        unnamed.push(0);
      }
      continue;
    }
    // an end on split copies its data names no unit of: on no copy, never a guessed one
    const flag = (ca > 1 ? 1 : 0) | (cb > 1 ? 2 : 0) | (ia < 0 || ib < 0 ? 4 : 0);
    if (flag & 3) {
      u.push(Number.NaN, Number.NaN);
      if (ca > 1) aCopies = Math.max(aCopies, ca);
      if (cb > 1) bCopies = Math.max(bCopies, cb);
    } else u.push(ua, ub);
    unnamed.push(flag & 3 ? flag : 0);
  }
  return {
    pieces,
    n: pair.length,
    u: Float64Array.from(u),
    pair: Int32Array.from(pair),
    unnamed: Uint8Array.from(unnamed),
    aCopies,
    bCopies,
  };
}
