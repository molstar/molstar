/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ribbons' laid-out form (_spec/README.md "Render", Bridges): a bridge's context, the lay a layer keeps, and what
// the layout, the draw and the hit share. Part of the bridge chunk.

import type { Want } from '@molstar/seqx-core/data';
import type { Location } from '@molstar/seqx-core/query';
import type { PairsData } from '@molstar/seqx-core/resolve';
import type { Channel } from '@molstar/seqx-core/theme/index';
import type { RenderConfig } from './config.js';
import type { Box, LayoutContext } from './layer-type.js';
import type { BridgePlacement } from './layout/bridge.js';
import type { PairsLay, PairsPrep, ScratchHit } from './pairs-lay.js';
import type { PairsRegion, RegionOver } from './pairs-region.js';
import type { SeqTokens } from './tokens.js';

/** What a layer laid out in a bridge reads besides a layout's context: its bridge (_spec/README.md "Render", Bridges). */
export interface BridgeContext {
  /** The bridge's band and sides. */
  readonly placement: BridgePlacement;
  /** Its index in `ResolvedView.bridges` (a hit's `bridge`) and its id (a hit's `section`). */
  readonly index: number;
  readonly id: string;
}

/** A layout's context in a bridge: `box` the band, `section` and `axis` its `a` side's, `bridge` the rest. */
export type BridgeLayoutContext = LayoutContext & {
  readonly bridge: BridgeContext;
  /** The layer's wants this frame (an `R` read is asked only when it fits beside them). */
  readonly wants?: readonly Want[];
};

/** A pair's end's system: its row's unit under a complex, else the layer's (the pairs layer's `aOf`, `bOf`). */
export function aOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.aUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.aSystem) : prep.aSystem;
}
export function bOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.bUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.bSystem) : prep.bSystem;
}

/** A ribbon drawn as one curve (its ends both under `minPx`, or `shape: line`), a stub (its far end not shown). */
export const Line = 1;
export const Stub = 2;

/** A ribbon's numbers: its `a` side's corners and y, then its `b` side's (`x0a, x1a, ya, x0b, x1b, yb`). */
export const R = 6;

export class RibbonsLay {
  key = '';
  /** The bridge's id and index (what a hit names). */
  bridgeId = '';
  bridge = -1;
  tokens: SeqTokens | undefined;
  config: RenderConfig | undefined;
  /** The frame's reads: the layer's prep, or a source's read on each side (`reads[src[i]]` ribbon `i`'s). */
  readonly reads: PairsPrep[] = [];
  /** Each read's side (0: the bridge's `a`); the ribbons it draws (1 its data's `a` end on its side, −1 turned, 0 any). */
  readonly readSide: number[] = [];
  readonly readOrient: number[] = [];
  /** The first read (the layer's prep with held pairs). */
  prep: PairsPrep | undefined;
  dpr = 1;
  colorChannel: Channel | undefined;
  /**
   * The ribbons in drawing order (widest first): corners, the pair (its first, merged), its read, its colour, flags, its
   * data's `a` end's side (1: a), whether blocks were merged into it, and its ends' hull on the data's `a` and `b`
   * (`aStart, aEnd, bStart, bEnd`: a merged ribbon's blocks together).
   */
  n = 0;
  rib = new Float64Array(6 * 64);
  pair = new Int32Array(64);
  src = new Uint8Array(64);
  color = new Uint32Array(64);
  flags = new Uint8Array(64);
  orient = new Int32Array(64);
  merged = new Uint8Array(64);
  hull = new Float64Array(4 * 64);
  /** A merged ribbon's extent, made on its first hit (a layout clears them). */
  readonly extents: (Location | undefined)[] = [];
  /** The band (the ribbons' clip) and its sides' edges. */
  band: Box = { x: 0, y: 0, width: 0, height: 0 };
  /** Ribbons past `maxArcs`: their ends' density a device column along each edge, from the band's x. */
  density = new Uint16Array(0);
  densityN = 0;
  densityMax = 0;
  edgeA = 0;
  edgeB = 0;
  /** Pairs with an end on a system neither side shows (`seq.render.bridge.unplaced`). */
  unplaced = 0;
  /** The stretches of an edge whose tiles are pending (EC-169): `x0, x1, y, down` each, clipped to the band. */
  readonly hatch: number[] = [];
  /**
   * Each read's state is kept by its object (the pairs layer's reader's key), its region by its system: a read a side,
   * or (`mirror` over a sided source) a read an end a side, `2 · end + side`; the near side's reads of `R` beside them.
   */
  readonly sides = [0, 1, 2, 3].map(() => ({}) as PairsLay);
  readonly rSides = [0, 1, 2, 3].map(() => ({}) as PairsLay);
  readonly regions = [0, 1, 2, 3].map((): [string, PairsRegion | undefined] => ['', undefined]);
  /** The `R` read's request (reused). */
  readonly over: RegionOver = { regions: new Map(), level: undefined, wants: [] };
  /** The drawn-once rule's ranges by data system: the near side's `R`, the far side's window. */
  readonly ownR = new Map<string, number[]>();
  readonly ownFar = new Map<string, number[]>();
  readonly scratch: ScratchHit = { layer: '', section: '', system: '' };
}

/**
 * Each pair's alignment from its key (RR's identity rule: a row read by region is keyed `r<row>|…`, an alignment's
 * blocks sharing `row`), −1 for none; undefined when no pair has one. Made once a data.
 */
export function groupsOf(d: PairsData): Int32Array | undefined {
  if (Groups.has(d)) return Groups.get(d);
  let out: Int32Array | undefined = new Int32Array(d.n).fill(-1);
  let any = false;
  for (let k = 0; k < d.n; k++) {
    const key = d.keys[k];
    if (typeof key !== 'string' || key.charCodeAt(0) !== 114) continue;
    const bar = key.indexOf('|');
    const row = bar > 1 ? Number(key.slice(1, bar)) : Number.NaN;
    if (!Number.isInteger(row) || row < 0) continue;
    out[k] = row;
    any = true;
  }
  if (!any) out = undefined;
  Groups.set(d, out);
  return out;
}
const Groups = new WeakMap<PairsData, Int32Array | undefined>();
