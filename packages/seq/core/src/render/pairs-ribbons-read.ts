/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge's reads of a source's pairs (_spec/README.md "Render", Bridges: what is read, and at which level): each side
// its own window at its own level, on the system it shows; the near side (the finer) its `R`, the drawn-once rule's
// ranges; each side's pending stretches the hatch along its edge (EC-169). Part of the bridge chunk.

import { BridgeHooks } from './bridge-hooks.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import type { BridgeSide } from './layout/bridge.js';
import type { PairsPrep } from './pairs-lay.js';
import type { PairsRegion } from './pairs-region.js';
import type { BridgeLayoutContext, RibbonsLay } from './pairs-ribbons-lay.js';

/** The near side's window and the far side's: what the drawn-once rule reads (`ownership`). */
export interface Ownership {
  /** The near side (0: the bridge's `a`). */
  readonly near: number;
  /** The near side's window widened outward to the far side's bins; the far side's window. */
  readonly r: readonly number[];
  readonly far: readonly number[];
}

/** A side's region: the layer's on the side's system (its name in the source), kept a side while that is the same. */
function sideRegion(lay: RibbonsLay, i: number, region: PairsRegion, name: string, system: string): PairsRegion {
  const kept = lay.regions[i];
  const key = `${name}\u0000${system}`;
  if (kept[0] === key && kept[1]) return kept[1];
  if (region.src.system === name && region.src.dataSystem === system) {
    kept[1] = region;
  } else {
    const length = lay.prep?.model.systems.get(system)?.length;
    const { length: _, ...rest } = region.src;
    kept[1] = {
      ...region,
      src: { ...rest, system: name, dataSystem: system, ...(length !== undefined ? { length } : {}) },
    };
  }
  kept[0] = key;
  return kept[1];
}

/**
 * The frame's reads of a layer in a bridge: its held pairs; a source's pairs by region on each side, each on the system
 * it shows at its own level, with the drawn-once rule's ranges (over an assembly or a side of blocks: one read, on the
 * bridge's `a` side, as the pairs layer reads it); the pending stretches of each side's edge.
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
  lay.hatch.length = 0;
  lay.prep = prep;
  const region = prep.region;
  const read = BridgeHooks.region;
  if (!region || prep.data || !read) {
    if (prep.data) {
      reads.push(prep);
      lay.readSide.push(0);
    }
    return undefined;
  }
  const { a: sa, b: sb } = ctx.bridge.placement;
  const perSide = !region.assembly && !sa.section.blocks && !sb.section.blocks;
  const regions = Regions;
  const sides = Sides;
  for (let i = 0; i < (perSide ? 2 : 1); i++) {
    const s = i ? sb : sa;
    regions[i] = undefined;
    if (!s.ready) continue;
    let r = region;
    if (perSide) {
      const sys = s.section.axis.system;
      const name = sys === prep.aSystem ? region.aName : sys === prep.bSystem ? region.bName : undefined;
      if (name === undefined) continue;
      r = sideRegion(lay, i, region, name, sys);
    }
    regions[i] = r;
    const side = (sides[i] = { ...ctx, section: s.section, axis: s.layout.axis!, place: s.layout.placement });
    const got = read(prep, r, side, data, lay.sides[i]);
    if (got?.data) {
      reads.push(got);
      lay.readSide.push(i);
    }
  }
  if (reads.length) lay.prep = reads[0];
  const own = reads.length === 2 ? ownership(lay, sa, sb) : undefined;
  // the near side reads `R` (its window and the widening to the far bins: a pair anchored there is its alone)
  if (own) {
    const again = read(prep, regions[own.near]!, sides[own.near]!, data, lay.sides[own.near], lay.ownR);
    if (again?.data) reads[own.near] = again;
    lay.prep = reads[0];
  }
  for (let i = 0; i < 2; i++) {
    const r = regions[i];
    const pending = r && !r.assembly ? BridgeHooks.state?.(lay.sides[i]).pending : undefined;
    if (pending?.length) hatchInto(lay, i ? sb : sa, r!.src.dataSystem, pending);
  }
  return own;
}

/** Each side's region and its layout context, as last read (kept: a frame's reads). */
const Regions: (PairsRegion | undefined)[] = [undefined, undefined];
const Sides: (LayoutContext | undefined)[] = [undefined, undefined];

/**
 * The drawn-once rule's ranges (_spec/README.md "Render", Bridges): the near side the finer read (items finest), the
 * bridge's `a` on a tie; `R` its window widened outward to the far level's grid; the far side's window.
 */
function ownership(lay: RibbonsLay, sa: BridgeSide, sb: BridgeSide): Ownership {
  const reads = lay.reads;
  const l0 = reads[0].level ?? 0;
  const l1 = reads[1].level ?? 0;
  const near = l1 < l0 ? 1 : 0;
  const step = reads[1 - near].level ?? 1;
  const r = lay.ownR;
  const far = lay.ownFar;
  r.length = 0;
  far.length = 0;
  (near ? sb : sa).layout.axis!.visible(r);
  for (let j = 0; j + 1 < r.length; j += 2) {
    r[j] = Math.floor(r[j] / step) * step;
    r[j + 1] = Math.ceil(r[j + 1] / step) * step;
  }
  (near ? sa : sb).layout.axis!.visible(far);
  return { near, r, far };
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
