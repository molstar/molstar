/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs layer in a bridge (_spec/README.md "Render", Bridges): a pair a ribbon between its ends on the two sides,
// its sides cubics with vertical tangents (`x_a + (x_b − x_a)·S(u)`, `S(u) = 3u² − 2u³`), each side joining corresponding
// positions (the a end's first to b's first when the strands' product is 1, to b's last when −1), so a ribbon is
// crossed exactly where GA's rule says, by geometry. Each end at every place of a block showing it; an end whose system
// a side shows on no block a stub (a quarter-height taper from the placed end); an end on a system neither side shows
// not drawn (`seq.render.bridge.unplaced`). Ends at least `ribbons.minPx` device px wide; both under it, one curve
// (`shape: line` forces it); an alignment's consecutive blocks under a device px apart on both sides one ribbon, their
// hull. Filled at `ribbons.fillOpacity` in the layer's colour (by default the `block` of the end on the bridge's `a`
// side, EC-218), widest first; past `pairs.maxArcs` the most aligned length kept, the rest a density strip along each
// edge. A source's pairs read by region on each side, at that side's level, each pair drawn once by ownership (the
// finer side's window widened to the far side's bins); a side's tiles pending the hatch along its edge (EC-169). The
// hit is closed form: at the point's height every ribbon's two sides are one multiply-add each, scanned last drawn
// first. Part of the bridge chunk: the lay (`pairs-ribbons-lay.ts`), the reads (`pairs-ribbons-read.ts`), the draw and
// the hit (`pairs-ribbons-draw.ts`).

import { defineIssues, issueOf } from '@molstar/core/repr/tree';
import { theme as themeForm } from '@molstar/seqx-core/lang';
import type { PairsData } from '@molstar/seqx-core/resolve';
import type { ChannelName, ThemeDatum } from '@molstar/seqx-core/theme/index';
import { channelLegends } from './channels.js';
import type { RenderConfig } from './config.js';
import { grown } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { LayerType, StoreReader } from './layer-type.js';
import type { BridgeSide } from './layout/bridge.js';
import type { PairsPrep } from './pairs-lay.js';
import { drawData, drawEmphasis, hitRibbons } from './pairs-ribbons-draw.js';
import { aOf, type BridgeLayoutContext, bOf, groupsOf, Line, R, RibbonsLay, Stub } from './pairs-ribbons-lay.js';
import { type Ownership, readSides } from './pairs-ribbons-read.js';
import { withAlpha } from './tokens.js';

export type { BridgeContext, BridgeLayoutContext } from './pairs-ribbons-lay.js';
export { hitRibbons, ribbonPath } from './pairs-ribbons-draw.js';
export { RibbonsLay };

/** The ribbons' codes, cataloged with the chunk (a composition that draws no bridge says none). */
export const RibbonIssues = defineIssues('seq.render', {
  'seq.render.bridge.unplaced': {
    severity: 'info',
    description: 'A bridge’s pairs have an end on a system neither side shows: those pairs are not drawn',
    hint: 'Show the pairs’ systems on the bridge’s sides (pairs are not projected through maps)',
  },
});

/** The default colour: the block of the end on the bridge's `a` side (chromosome painting, EC-218). */
const DefaultColor = themeForm.field('block', { list: false });

/**
 * A pair as a ribbon's colour reads it: the pairs layer's datum (its own fields), and `block`, `system` of its end on the
 * bridge's `a` side. Wraps the prep's datum: the chunk imports nothing of the pairs layer's modules at run time (a
 * module shared with the entry would be loaded with it).
 */
class RibbonDatum implements ThemeDatum {
  value = Number.NaN;
  letter = '';
  start = 0;
  end = 0;
  block = '';
  system = '';
  constructor(readonly pair: ThemeDatum & { set(k: number): ThemeDatum }) {}
  set(k: number): this {
    this.pair.set(k);
    this.start = this.pair.start;
    this.end = this.pair.end;
    return this;
  }
  field(name: string): string | number | undefined {
    if (name === 'block') return this.block;
    if (name === 'system') return this.system;
    return this.pair.field(name);
  }
}

/** Whether a field holds numbers: the ends and the numeric extra columns (the pairs layer's `pairEnv`), kept a data. */
const Numeric = new WeakMap<PairsData, (n: string) => boolean>();
function numericOf(d: PairsData): (n: string) => boolean {
  let f = Numeric.get(d);
  if (f) return f;
  const numbers = new Set(
    Object.entries(d.extra ?? {}).flatMap(([k, v]) => (v.some((x) => typeof x === 'number') ? [k] : [])),
  );
  f = (n) => n === 'a_start' || n === 'a_end' || n === 'b_start' || n === 'b_end' || numbers.has(n);
  Numeric.set(d, f);
  return f;
}

// ---- candidates (kept: grown, never shrunk) ----

let Cand = new Float64Array(6 * 256);
let CandPair = new Int32Array(256);
let CandRead = new Uint8Array(256);
let CandFlags = new Uint8Array(256);
let CandOrient = new Int32Array(256);
let CandBlock = new Int32Array(256);
let CandMerged = new Uint8Array(256);
let CandHull = new Float64Array(4 * 256);
let CandWidth = new Float64Array(256);
let CandRank = new Float64Array(256);
let Order = new Int32Array(256);
let nc = 0;

function room(n: number) {
  if (CandPair.length >= n) return;
  const m = Math.max(n, 2 * CandPair.length);
  Cand = grown(Cand, R * m);
  CandPair = grown(CandPair, m);
  CandRead = grown(CandRead, m);
  CandFlags = grown(CandFlags, m);
  CandOrient = grown(CandOrient, m);
  CandBlock = grown(CandBlock, m);
  CandMerged = grown(CandMerged, m);
  CandHull = grown(CandHull, 4 * m);
  CandWidth = grown(CandWidth, m);
  CandRank = grown(CandRank, m);
  Order = grown(Order, m);
}

/** Each end's places (x pairs) and blocks on its side. */
const PlacesA = new Float64Array(64);
const BlocksA = new Int32Array(32);
const PlacesB = new Float64Array(64);
const BlocksB = new Int32Array(32);
const MaxPlaces = 32;
/** Each place's positions (`extentInto`'s `pos`), and a b piece's places within an a piece's image. */
const PosA = new Float64Array(64);
const PosB = new Float64Array(64);
const PosC = new Float64Array(64);
const PlacesC = new Float64Array(64);
const BlocksC = new Int32Array(32);

/** An end narrower than `min` px widened about its middle (into `p[o]`, `p[o + 1]`), its order kept; true when it was. */
function widen(p: Float64Array, o: number, min: number): boolean {
  const w = p[o + 1] - p[o];
  if (Math.abs(w) >= min) return false;
  const m = (p[o] + p[o + 1]) / 2;
  const h = (w < 0 ? -min : min) / 2;
  p[o] = m - h;
  p[o + 1] = m + h;
  return true;
}

/** The read being laid: its data, its index in the lay's reads. */
let Data: PairsData;
let Read = 0;

/** A candidate ribbon from the two ends' places (one side's corners in its order, the other's in the run's). */
function push(
  xa0: number,
  xa1: number,
  ya: number,
  xb0: number,
  xb1: number,
  yb: number,
  k: number,
  flags: number,
  orient: number,
  block: number,
  rank: number,
) {
  room(nc + 1);
  const o = R * nc;
  Cand[o] = xa0;
  Cand[o + 1] = xa1;
  Cand[o + 2] = ya;
  Cand[o + 3] = xb0;
  Cand[o + 4] = xb1;
  Cand[o + 5] = yb;
  CandPair[nc] = k;
  CandRead[nc] = Read;
  CandFlags[nc] = flags;
  CandOrient[nc] = orient;
  CandBlock[nc] = block;
  CandMerged[nc] = 0;
  CandHull[4 * nc] = Data.aStart[k];
  CandHull[4 * nc + 1] = Data.aEnd[k];
  CandHull[4 * nc + 2] = Data.bStart[k];
  CandHull[4 * nc + 3] = Data.bEnd[k];
  CandWidth[nc] = Math.abs(xa1 - xa0) + Math.abs(xb1 - xb0);
  CandRank[nc] = rank;
  nc++;
  Run.c = -1;
}

// ---- the drawn-once rule (_spec/README.md "Render", Bridges: each pair drawn once, by ownership) ----

/** While a read is laid with two reads: whether it is the near side's, the near side, `R` and the far window. */
const Own = { on: false, mine: false, near: 0, r: [] as readonly number[], far: [] as readonly number[] };

/** Whether `p` lies in `[start, end, …]`. */
function inRanges(r: readonly number[], p: number): boolean {
  for (let i = 0; i + 1 < r.length; i += 2) if (p >= r[i] && p < r[i + 1]) return true;
  return false;
}

/**
 * Whether the read being laid draws a ribbon with these anchors on the bridge's sides (each end's first position): the
 * near side every one anchored in `R`; the far side those anchored in its window whose near anchor is outside `R`.
 */
function owns(anchorA: number, anchorB: number): boolean {
  const near = Own.near ? anchorB : anchorA;
  if (Own.mine) return inRanges(Own.r, near);
  return inRanges(Own.far, Own.near ? anchorA : anchorB) && !inRanges(Own.r, near);
}

// ---- an alignment's blocks merged under a pixel ----

/** The open run of merged blocks: its candidate (−1 none), what it continues, its raw corners (before widening). */
const Run = { c: -1, group: -1, orient: 0, turned: false, ba: 0, bb: 0, a0: 0, a1: 0, b0: 0, b1: 0 };
/** The read's alignments (items only), and a device px. */
let Groups: Int32Array | undefined;
let Px = 1;

/** The gap between a hull and an end (negative when they overlap). */
function gap(h0: number, h1: number, x0: number, x1: number): number {
  const hl = h0 < h1 ? h0 : h1;
  const hh = h0 < h1 ? h1 : h0;
  const l = x0 < x1 ? x0 : x1;
  const h = x0 < x1 ? x1 : x0;
  return Math.max(l - hh, hl - h);
}

/** A hull grown over an end, its direction kept (into `Corner[o]`, `Corner[o + 1]`). */
function hullOver(h0: number, h1: number, x0: number, x1: number, o: number) {
  const lo = Math.min(h0, h1, x0, x1);
  const hi = Math.max(h0, h1, x0, x1);
  Corner[o] = h0 <= h1 ? lo : hi;
  Corner[o + 1] = h0 <= h1 ? hi : lo;
}

/**
 * A block's ribbon merged into the open run when it continues it (its alignment, orientation, strand and blocks) under
 * a device px apart on both sides: the run's hull its corners, its rank the sum. True when merged.
 */
function merged(
  k: number,
  orient: number,
  turned: boolean,
  a0: number,
  a1: number,
  b0: number,
  b1: number,
  min: number,
  line: boolean,
  rank: number,
): boolean {
  const g = Groups![k];
  if (g < 0) return false;
  const r = Run;
  if (
    r.c < 0 ||
    r.group !== g ||
    r.orient !== orient ||
    r.turned !== turned ||
    r.ba !== BlocksA[0] ||
    r.bb !== BlocksB[0] ||
    gap(r.a0, r.a1, a0, a1) >= Px ||
    gap(r.b0, r.b1, b0, b1) >= Px
  )
    return false;
  hullOver(r.a0, r.a1, a0, a1, 0);
  hullOver(r.b0, r.b1, b0, b1, 2);
  r.a0 = Corner[0];
  r.a1 = Corner[1];
  r.b0 = Corner[2];
  r.b1 = Corner[3];
  const thinA = widen(Corner, 0, min);
  const thinB = widen(Corner, 2, min);
  const c = r.c;
  const o = R * c;
  Cand[o] = Corner[0];
  Cand[o + 1] = Corner[1];
  Cand[o + 3] = Corner[2];
  Cand[o + 4] = Corner[3];
  CandFlags[c] = line || (thinA && thinB) ? Line : 0;
  CandWidth[c] = Math.abs(Corner[1] - Corner[0]) + Math.abs(Corner[3] - Corner[2]);
  CandRank[c] += rank;
  CandMerged[c] = 1;
  const h = 4 * c;
  CandHull[h] = Math.min(CandHull[h], Data.aStart[k]);
  CandHull[h + 1] = Math.max(CandHull[h + 1], Data.aEnd[k]);
  CandHull[h + 2] = Math.min(CandHull[h + 2], Data.bStart[k]);
  CandHull[h + 3] = Math.max(CandHull[h + 3], Data.bEnd[k]);
  return true;
}

/**
 * One pair's ribbons (`orient` 1: its `a` end on the bridge's `a` side; −1: turned), every place of each end; a stub
 * where one end's system is shown on no block of its side. False when neither end is placed.
 */
function pairRibbons(
  d: PairsData,
  k: number,
  orient: 1 | -1,
  sa: BridgeSide,
  sb: BridgeSide,
  aSys: string,
  bSys: string,
  min: number,
  stub: number,
  line: boolean,
  rank: number,
): boolean {
  // the ends as the sides take them: the bridge's `a` side holds the data's `a` end unless turned
  const as = orient > 0 ? d.aStart[k] : d.bStart[k];
  const ae = orient > 0 ? d.aEnd[k] : d.bEnd[k];
  const bs = orient > 0 ? d.bStart[k] : d.aStart[k];
  const be = orient > 0 ? d.bEnd[k] : d.aEnd[k];
  // another read draws it (or neither: in neither view)
  if (Own.on && !owns(as, bs)) return true;
  const na = sa.extentInto(orient > 0 ? aSys : bSys, as, ae, PlacesA, BlocksA, MaxPlaces, PosA);
  const nb = sb.extentInto(orient > 0 ? bSys : aSys, bs, be, PlacesB, BlocksB, MaxPlaces, PosB);
  if (!na && !nb) return false;
  // corresponding positions: a's first to b's first when the strands' product is 1, to b's last when −1 (GA's twist)
  const product = d.aStrand[k] * d.bStrand[k];
  const ya = sa.edge;
  const yb = sb.edge;
  const yStubB = ya + stub * (yb - ya);
  if (!na || !nb) {
    // a stub: a taper from the placed end into the band, `stubShare` of its height
    const near = na ? sa : sb;
    const n = na || nb;
    const xs = na ? PlacesA : PlacesB;
    const yNear = near.edge;
    const yFar = yNear + stub * ((na ? yb : ya) - yNear);
    for (let i = 0; i < n; i++) {
      widen(xs, 2 * i, min);
      const m = (xs[2 * i] + xs[2 * i + 1]) / 2;
      if (na) push(xs[0 + 2 * i], xs[2 * i + 1], yNear, m, m, yFar, k, Stub, orient, BlocksA[i], rank);
      else push(m, m, yFar, xs[2 * i], xs[2 * i + 1], yNear, k, Stub, orient, -1, rank);
    }
    return true;
  }
  // the b corners in the order corresponding to a's: turned where the run is against
  const turned = product < 0;
  if (na === 1 && nb === 1 && PosA[0] === as && PosA[1] === ae && PosB[0] === bs && PosB[1] === be && Groups) {
    const a0 = PlacesA[0];
    const a1 = PlacesA[1];
    const b0 = turned ? PlacesB[1] : PlacesB[0];
    const b1 = turned ? PlacesB[0] : PlacesB[1];
    if (merged(k, orient, turned, a0, a1, b0, b1, min, line, rank)) return true;
    push1(a0, a1, ya, b0, b1, yb, k, orient, 0, min, line, rank);
    const r = Run;
    r.c = nc - 1;
    r.group = Groups[k];
    r.orient = orient;
    r.turned = turned;
    r.ba = BlocksA[0];
    r.bb = BlocksB[0];
    r.a0 = a0;
    r.a1 = a1;
    r.b0 = b0;
    r.b1 = b1;
    return true;
  }
  // each a piece joined to the b positions corresponding to it through the pair's own linear map (reversed when
  // turned), cut by b's pieces in turn: every ribbon joins positions that correspond; an a piece whose image no block
  // shows is a stub
  const kb = ae > as ? (be - bs) / (ae - as) : 0;
  const bSide = orient > 0 ? bSys : aSys;
  for (let i = 0; i < na; i++) {
    const s = PosA[2 * i];
    const e = PosA[2 * i + 1];
    const xa0 = PlacesA[2 * i];
    const xa1 = PlacesA[2 * i + 1];
    const i0 = turned ? be - (e - as) * kb : bs + (s - as) * kb;
    const i1 = turned ? be - (s - as) * kb : bs + (e - as) * kb;
    const m = sb.extentInto(bSide, Math.floor(i0), Math.ceil(i1), PlacesC, BlocksC, MaxPlaces, PosC);
    if (!m) {
      Corner[0] = xa0;
      Corner[1] = xa1;
      widen(Corner, 0, min);
      const c = (Corner[0] + Corner[1]) / 2;
      push(Corner[0], Corner[1], ya, c, c, yStubB, k, Stub, orient, BlocksA[i], rank);
      continue;
    }
    for (let j = 0; j < m; j++) {
      const p0 = PosC[2 * j];
      const p1 = PosC[2 * j + 1];
      // the a positions corresponding to `[p0, p1)` on b, within the a piece
      let q0 = kb ? (turned ? as + (be - p1) / kb : as + (p0 - bs) / kb) : s;
      let q1 = kb ? (turned ? as + (be - p0) / kb : as + (p1 - bs) / kb) : e;
      q0 = q0 < s ? s : q0 > e ? e : q0;
      q1 = q1 < s ? s : q1 > e ? e : q1;
      push1(
        lerpX(xa0, xa1, s, e, q0),
        lerpX(xa0, xa1, s, e, q1),
        ya,
        turned ? PlacesC[2 * j + 1] : PlacesC[2 * j],
        turned ? PlacesC[2 * j] : PlacesC[2 * j + 1],
        yb,
        k,
        orient,
        i,
        min,
        line,
        rank,
      );
    }
  }
  return true;
}

/** The x of position `p` within a piece `[s, e)` drawn from `x0` to `x1` (a block piece is linear). */
function lerpX(x0: number, x1: number, s: number, e: number, p: number): number {
  return e > s ? x0 + ((p - s) * (x1 - x0)) / (e - s) : x0;
}

/** A ribbon of two placed ends, each widened to `min` (both under it: a line). */
function push1(
  a0: number,
  a1: number,
  ya: number,
  b0: number,
  b1: number,
  yb: number,
  k: number,
  orient: number,
  i: number,
  min: number,
  line: boolean,
  rank: number,
) {
  Corner[0] = a0;
  Corner[1] = a1;
  Corner[2] = b0;
  Corner[3] = b1;
  const thinA = widen(Corner, 0, min);
  const thinB = widen(Corner, 2, min);
  const flags = line || (thinA && thinB) ? Line : 0;
  push(Corner[0], Corner[1], ya, Corner[2], Corner[3], yb, k, flags, orient, BlocksA[i], rank);
}

const Corner = new Float64Array(4);

/** The systems neither side shows, the first few (said once). */
const UnplacedSystems: string[] = [];

/**
 * Lays out a layer's ribbons in a bridge from its reads (one, or a read a side with `own`): each pair placed by the
 * systems its sides show (both: the data's `a` end on the bridge's `a` side; `mirror` the swap too), drawn once by
 * ownership, an alignment's blocks merged under a pixel, ranked past `maxArcs`, the kept ones widest first.
 */
export function layoutRibbons(lay: RibbonsLay, ctx: BridgeLayoutContext, own?: Ownership): void {
  const reads = lay.reads;
  const { placement } = ctx.bridge;
  const sa = placement.a;
  const sb = placement.b;
  lay.band = placement.box;
  lay.edgeA = sa.edge;
  lay.edgeB = sb.edge;
  lay.extents.length = 0;
  if (!reads.length || !sa.ready || !sb.ready) return;
  const config = ctx.config;
  const min = config.ribbons.minPx / ctx.dpr;
  const stub = config.ribbons.stubShare;
  const first = reads[0];
  const line = (first.channels.forms.shape as unknown) === 'line';
  Px = 1 / ctx.dpr;
  nc = 0;
  UnplacedSystems.length = 0;
  let unplaced = 0;
  for (let r = 0; r < reads.length; r++) {
    const prep = reads[r];
    const d = prep.data!;
    Data = d;
    Read = r;
    Run.c = -1;
    Groups = prep.level === undefined ? groupsOf(d) : undefined;
    Own.on = own !== undefined;
    if (own) {
      Own.mine = lay.readSide[r] === own.near;
      Own.near = own.near;
      Own.r = own.r;
      Own.far = own.far;
    }
    for (let k = 0; k < d.n; k++) {
      if (!prep.shown[k]) continue;
      const as = aOf(prep, k);
      const bs = bOf(prep, k);
      // an end on a chromosome the side's `include` leaves out is the side's still: a stub (Defaults, 4)
      const fwd = sa.holds(as) && sb.holds(bs) && (sa.shows(as) || sb.shows(bs));
      const back = sa.holds(bs) && sb.holds(as) && (sa.shows(bs) || sb.shows(as));
      const rank = prep.weight ? prep.weight[k] : d.aEnd[k] - d.aStart[k];
      if (fwd) pairRibbons(d, k, 1, sa, sb, as, bs, min, stub, line, rank);
      // turned when the sides show only that way, or (one genome against itself) the mirror's swap too
      if (back && (!fwd || (prep.mirror && (as !== bs || d.aStart[k] !== d.bStart[k]))))
        pairRibbons(d, k, -1, sa, sb, as, bs, min, stub, line, rank);
      if (fwd || back) continue;
      unplaced++;
      for (const s of [as, bs])
        if (!sa.shows(s) && !sb.shows(s) && UnplacedSystems.length < 3 && !UnplacedSystems.includes(s))
          UnplacedSystems.push(s);
    }
  }
  Own.on = false;
  Groups = undefined;
  lay.unplaced = unplaced;
  if (unplaced)
    ctx.report(
      issueOf(
        RibbonIssues,
        'seq.render.bridge.unplaced',
        `${unplaced} pair(s) have an end on a system neither side shows${UnplacedSystems.length ? ` (${UnplacedSystems.map((s) => s.slice(s.lastIndexOf('/') + 1)).join(', ')})` : ''}: not drawn`,
        { address: ctx.key, hint: RibbonIssues.codes['seq.render.bridge.unplaced'].hint },
      ),
    );
  // past `maxArcs`: the most aligned length kept (ties by span), the rest a density strip along each edge
  const max = config.pairs.maxArcs;
  for (let i = 0; i < nc; i++) Order[i] = i;
  const order = Order.subarray(0, nc);
  let kept = nc;
  lay.densityN = 0;
  lay.densityMax = 0;
  if (nc > max) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.pairs.limit', `More than ${max} ribbons in view: the most aligned kept`, {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.pairs.limit'].hint,
      }),
    );
    order.sort((x, y) => CandRank[y] - CandRank[x] || CandWidth[y] - CandWidth[x]);
    kept = max;
    densityOf(lay, order, max, ctx.dpr);
  }
  // widest first, so narrow ones lie on top and stay hittable
  const top = order.subarray(0, kept);
  top.sort((x, y) => CandWidth[y] - CandWidth[x] || x - y);
  lay.n = 0;
  lay.rib = grown(lay.rib, R * kept);
  lay.pair = grown(lay.pair, kept);
  lay.src = grown(lay.src, kept);
  lay.color = grown(lay.color, kept);
  lay.flags = grown(lay.flags, kept);
  lay.orient = grown(lay.orient, kept);
  lay.merged = grown(lay.merged, kept);
  lay.hull = grown(lay.hull, 4 * kept);
  const a = sa.section;
  const domain = a.domains?.block ?? [a.address.slice(a.address.lastIndexOf('/') + 1)];
  lay.colorChannel = first.channels.get('color' as ChannelName, DefaultColor, ctx, {
    numeric: numericOf(first.data!),
    categories: (f) => (f === 'block' ? domain : f === 'system' ? a.domains?.system : undefined),
  });
  const fill = config.ribbons.fillOpacity;
  // the bins' opacity scale: each read's heaviest bin
  const heaviest = Heaviest;
  for (let r = 0; r < reads.length && r < heaviest.length; r++) heaviest[r] = maxOf(reads[r].weight);
  for (let oi = 0; oi < kept; oi++) {
    const c = top[oi];
    const i = lay.n++;
    for (let r = 0; r < R; r++) lay.rib[R * i + r] = Cand[R * c + r];
    for (let r = 0; r < 4; r++) lay.hull[4 * i + r] = CandHull[4 * c + r];
    const k = CandPair[c];
    const read = CandRead[c];
    const prep = reads[read];
    lay.pair[i] = k;
    lay.src[i] = read;
    lay.flags[i] = CandFlags[c];
    lay.orient[i] = CandOrient[c];
    lay.merged[i] = CandMerged[c];
    const b = CandBlock[c];
    const block = b >= 0 ? a.blocks?.[b] : undefined;
    const datum = ribbonDatum(prep);
    datum.set(k);
    datum.block = block ? block.label : domain[0];
    datum.system = block ? (block.name ?? block.axis.system) : a.axis.system;
    const v = lay.colorChannel.at(datum);
    const color = typeof v === 'number' ? v : lay.tokens!.colors['seq-range'];
    const alpha = ((color & 0xff) / 255) * (prep.weight ? binAlpha(prep.weight[k], heaviest[read], config) : fill);
    lay.color[i] = withAlpha(color, alpha);
  }
}

const Heaviest = new Float64Array(2);
function maxOf(w: Float32Array | undefined): number {
  let max = 0;
  if (w) for (let i = 0; i < w.length; i++) if (w[i] > max) max = w[i];
  return max;
}

/** A bin's opacity by its aligned length (its weight) on a square-root scale over its read's bins. */
function binAlpha(w: number, max: number, config: RenderConfig): number {
  const { min, max: hi } = config.ribbons.binOpacity;
  return max > 0 ? min + (hi - min) * Math.sqrt(w / max) : hi;
}

/** The kept prep's datum for ribbons (one a prep). */
const Datums = new WeakMap<PairsPrep, RibbonDatum>();
function ribbonDatum(prep: PairsPrep): RibbonDatum {
  let d = Datums.get(prep);
  if (!d) Datums.set(prep, (d = new RibbonDatum(prep.datum!)));
  return d;
}

/** The ribbons past `max` as their ends' density a device column along each edge. */
function densityOf(lay: RibbonsLay, order: Int32Array, max: number, dpr: number) {
  const b = lay.band;
  const cols = Math.ceil(b.width * dpr) + 1;
  if (lay.density.length < 2 * cols) lay.density = new Uint16Array(2 * cols);
  lay.density.fill(0, 0, 2 * cols);
  lay.densityN = cols;
  for (let oi = max; oi < order.length; oi++) {
    const o = R * order[oi];
    for (let side = 0; side < 2; side++) {
      const x = (Cand[o + 3 * side] + Cand[o + 3 * side + 1]) / 2;
      const c = Math.floor((x - b.x) * dpr);
      if (c < 0 || c >= cols) continue;
      const k = side * cols + c;
      const v = lay.density[k] < 65535 ? ++lay.density[k] : 65535;
      if (v > lay.densityMax) lay.densityMax = v;
    }
  }
}

/**
 * The pairs layer in a bridge (`scope: 'bridge'`, internal): laid out once a bridge with its `BridgeLayoutContext`,
 * drawn on the page, clipped to its band.
 */
export const RibbonsType: LayerType<PairsPrep, RibbonsLay> = {
  name: 'pairs',
  fallback: 'ranges',
  scope: 'bridge',
  ownHover: true,
  // the pairs layer's prepared form, made by the composition's `pairs` type (`layoutBridges`): never called
  prepare: () => undefined as unknown as PairsPrep,
  layout(prep, ctx, data: StoreReader, _overlay, prev) {
    const lay = prev ?? new RibbonsLay();
    const c = ctx as BridgeLayoutContext;
    lay.key = ctx.key;
    lay.bridgeId = c.bridge.id;
    lay.bridge = c.bridge.index;
    lay.tokens = ctx.tokens;
    lay.config = ctx.config;
    lay.dpr = ctx.dpr;
    lay.n = 0;
    lay.densityN = 0;
    lay.band = c.bridge.placement.box;
    const own = readSides(lay, prep, c, data);
    layoutRibbons(lay, c, own);
    return lay;
  },
  draw(lay, pass, out, emphasis) {
    if (!lay.tokens) return;
    if (pass === 'data') {
      if (lay.prep || lay.hatch.length) drawData(lay, out);
    } else if (lay.prep && emphasis && (pass === 'hover' || pass === 'selection'))
      drawEmphasis(lay, pass, out, emphasis);
  },
  hit: hitRibbons,
  legend(lay, options) {
    return channelLegends([lay.colorChannel], { ...options, layer: lay.key });
  },
  filtered: (lay) => lay.prep?.filtered ?? 0,
};
