/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge's reads of a source's pairs (_spec/README.md "Render", Bridges: what is read, and at which level): each side
// its own window at its own level, on the systems it shows (over assemblies, each block's own); the near side (the
// finer) its `R`, the drawn-once rule's ranges by system; each side's pending stretches the hatch along its edge
// (EC-169). Part of the bridge chunk.

import { BridgeHooks } from './bridge-hooks.js';
import { widenedToGrid } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import type { BridgeSide } from './layout/bridge.js';
import type { PairsPrep } from './pairs-lay.js';
import type { PairsRegion } from './pairs-region.js';
import type { BridgeLayoutContext, RibbonsLay } from './pairs-ribbons-lay.js';

/** The near side's window and the far side's: what the drawn-once rule reads (`ownership`), by data system. */
export interface Ownership {
  /** The near side (0: the bridge's `a`). */
  readonly near: number;
  /**
   * The near side's window widened outward to the far side's bins (shrunk to them while its `R` read is not whole);
   * the far side's window (sorted, merged).
   */
  readonly r: ReadonlyMap<string, readonly number[]>;
  readonly far: ReadonlyMap<string, readonly number[]>;
  /** The far level (`R`'s grid). */
  readonly step: number;
  /**
   * While `R`'s read is not whole: the near side's windows, where it also draws a pair whose far end's far bin is not
   * in the far window (the far side draws none of those).
   */
  win?: ReadonlyMap<string, readonly number[]>;
}

/**
 * A side's region: the layer's on the side's system (its name in the source; a sided source's side the data's end the
 * side shows), kept a side while that is the same.
 */
function sideRegion(
  lay: RibbonsLay,
  i: number,
  region: PairsRegion,
  name: string,
  system: string,
  side: 'a' | 'b' | undefined,
): PairsRegion {
  const kept = lay.regions[i];
  const key = `${name}\u0000${system}\u0000${side}`;
  if (kept[0] === key && kept[1]) return kept[1];
  if (region.src.system === name && region.src.dataSystem === system && region.src.side === side) {
    kept[1] = region;
  } else {
    const length = lay.prep?.model.systems.get(system)?.length;
    const { length: _, side: __, ...rest } = region.src;
    kept[1] = {
      ...region,
      src: {
        ...rest,
        system: name,
        dataSystem: system,
        ...(side ? { side } : {}),
        ...(length !== undefined ? { length } : {}),
      },
    };
  }
  kept[0] = key;
  return kept[1];
}

/**
 * The frame's reads of a layer in a bridge: its held pairs; a source's pairs by region on each side, each on the
 * system(s) it shows at its own level (over assemblies a read a block of its own; a sided source on the end the side
 * shows; `mirror` over one genome's sided source, both ends a side, each read drawing the ribbons whose end on its side
 * is the end it read), with the drawn-once rule's ranges; the pending stretches of each side's edge.
 */
export function readSides(
  lay: RibbonsLay,
  prep: PairsPrep,
  ctx: BridgeLayoutContext,
  data: StoreReader,
): Ownership | undefined {
  const reads = lay.reads;
  reads.length = 0;
  lay.readSide.length = 0;
  lay.readOrient.length = 0;
  lay.hatch.length = 0;
  lay.prep = prep;
  const region = prep.region;
  const read = BridgeHooks.region;
  const state = BridgeHooks.state;
  if (!region || prep.data || !read || !state) {
    if (prep.data) {
      reads.push(prep);
      lay.readSide.push(0);
      lay.readOrient.push(0);
    }
    return undefined;
  }
  const { a: sa, b: sb } = ctx.bridge.placement;
  // a side of blocks over no assembly: one read, on `a`'s side (its blocks' systems are not the data's)
  const perSide = !!region.assembly || (!sa.section.blocks && !sb.section.blocks);
  const sided = !!region.src.side || !!data.coverage(region.src.source)?.systems.some((x) => x.side);
  const both = perSide && sided && prep.mirror && prep.aSystem === prep.bSystem;
  const regions = Regions;
  const sides = Sides;
  Jobs.length = 0;
  for (let j = 0; j < (perSide ? (both ? 4 : 2) : 1); j++) {
    const i = j & 1;
    const s = i ? sb : sa;
    regions[j] = undefined;
    if (!s.ready) continue;
    let end: 'a' | 'b' | undefined = both ? (j < 2 ? 'a' : 'b') : undefined;
    let r = region;
    if (perSide && !region.assembly) {
      const sys = s.section.axis.system;
      // the data's end the side shows (one system on both: the bridge's side)
      end ??= sys === prep.aSystem && (sys !== prep.bSystem || !i) ? 'a' : sys === prep.bSystem ? 'b' : undefined;
      if (end === undefined) continue;
      r = sideRegion(lay, j, region, end === 'a' ? region.aName : region.bName, sys, sided ? end : undefined);
    }
    regions[j] = r;
    Ends[j] = end ?? (i ? 'b' : 'a');
    const side = (sides[j] = { ...ctx, section: s.section, axis: s.layout.axis!, place: s.layout.placement });
    const got = read(prep, r, side, data, lay.sides[j], undefined, Ends[j]);
    // a side with nothing read (its tiles pending or failed) takes no part: the other draws all it reads
    if (got?.data?.n) {
      reads.push(got);
      lay.readSide.push(i);
      lay.readOrient.push(both ? (!i === (end === 'a') ? 1 : -1) : 0);
      Jobs.push(j);
    }
  }
  if (reads.length) lay.prep = reads[0];
  const own = lay.readSide.includes(0) && lay.readSide.includes(1) ? ownership(lay, state) : undefined;
  if (own) {
    // the near side reads `R` at its own level (its window and the widening to the far bins: a pair anchored there is
    // its alone), a state of its own beside its window's; while that read is not whole (pending, failed, past a bound:
    // then never asked), `R` is its window shrunk to the far grid and the near side draws the rest of its window the far
    // side does not, so nothing is owned and unread, nothing drawn twice
    const o = lay.over;
    o.regions = own.r;
    o.wants = ctx.wants ?? [];
    let whole = true;
    for (let k = 0; k < reads.length; k++) {
      Again[k] = undefined;
      if (lay.readSide[k] !== own.near) continue;
      const j = Jobs[k];
      o.level = reads[k].level;
      Again[k] = read(prep, regions[j]!, sides[j]!, data, lay.rSides[j], o, Ends[j]);
      if (!Again[k]?.data) whole = false;
    }
    if (whole) for (let k = 0; k < reads.length; k++) reads[k] = Again[k] ?? reads[k];
    else {
      own.win = windowsOf(lay, state, own.near);
      rangesInto(own.win, own.step, lay.ownR, true);
    }
    lay.prep = reads[0];
  }
  for (let j = 0; j < 4; j++) {
    const r = regions[j];
    const pending = r && !r.assembly ? state(lay.sides[j]).pending : undefined;
    if (pending?.length) hatchInto(lay, j & 1 ? sb : sa, r!.src.dataSystem, pending);
  }
  return own;
}

/** Each read's region, data end and layout context, as last read (kept: a frame's reads); the reads' jobs; `R`'s. */
const Regions: (PairsRegion | undefined)[] = [undefined, undefined, undefined, undefined];
const Ends: ('a' | 'b')[] = ['a', 'b', 'a', 'b'];
const Sides: (LayoutContext | undefined)[] = [undefined, undefined, undefined, undefined];
const Jobs: number[] = [];
const Again: (PairsPrep | undefined)[] = [];

/** A side's windows by data system (its first read's: a side's reads share them). */
function windowsOf(lay: RibbonsLay, state: NonNullable<(typeof BridgeHooks)['state']>, side: number) {
  return state(lay.sides[Jobs[lay.readSide.indexOf(side)]]).windows;
}

/**
 * The drawn-once rule's ranges (_spec/README.md "Render", Bridges): the near side the finer (its coarsest read; items
 * finest), the bridge's `a` on a tie; `R` its windows (each system's as its read read them: `dataRegions`, sorted and
 * merged) widened outward to the far level's grid; the far side's windows. A few ranges a system (EC-33).
 */
function ownership(lay: RibbonsLay, state: NonNullable<(typeof BridgeHooks)['state']>): Ownership {
  const l = Levels;
  l[0] = l[1] = 0;
  for (let k = 0; k < lay.reads.length; k++) l[lay.readSide[k]] = Math.max(l[lay.readSide[k]], lay.reads[k].level ?? 0);
  const near = l[1] < l[0] ? 1 : 0;
  const step = l[1 - near] || 1;
  rangesInto(windowsOf(lay, state, near), step, lay.ownR);
  rangesInto(windowsOf(lay, state, 1 - near), 1, lay.ownFar);
  return { near, r: lay.ownR, far: lay.ownFar, step };
}
const Levels = [0, 0];

/** Each system's sorted ranges on the grid of `step` (widened, or `inward` shrunk), into `out` (its arrays reused). */
function rangesInto(
  windows: ReadonlyMap<string, readonly number[]>,
  step: number,
  out: Map<string, number[]>,
  inward = false,
) {
  out.forEach(clear);
  Into.out = out;
  Into.step = step;
  Into.inward = inward;
  windows.forEach(into);
}
const Into = { out: undefined as unknown as Map<string, number[]>, step: 1, inward: false };
function clear(r: number[]) {
  r.length = 0;
}
function into(w: readonly number[], system: string) {
  let r = Into.out.get(system);
  if (!r) Into.out.set(system, (r = []));
  widenedToGrid(w, Into.step, r, Into.inward);
}

/** Places of a pending stretch (x pairs) and their blocks. */
const Places = new Float64Array(64);
const Blocks = new Int32Array(32);

/** A side's pending stretches (`[start, end, …]` on `system`) as hatch stretches along its edge, clipped to the band. */
function hatchInto(lay: RibbonsLay, side: BridgeSide, system: string, pending: readonly number[]) {
  const b = lay.band;
  const down = side.edge <= b.y + b.height / 2;
  for (let i = 0; i + 1 < pending.length; i += 2) {
    const n = side.extentInto(system, pending[i], pending[i + 1], Places, Blocks, 32);
    for (let j = 0; j < n; j++) {
      const x0 = Math.max(b.x, Math.min(Places[2 * j], Places[2 * j + 1]));
      const x1 = Math.min(b.x + b.width, Math.max(Places[2 * j], Places[2 * j + 1]));
      if (x1 > x0) lay.hatch.push(x0, x1, side.edge, down ? 1 : 0);
    }
  }
}
