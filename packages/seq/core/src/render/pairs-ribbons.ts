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
// (`shape: line` forces it). Filled at `ribbons.fillOpacity` in the layer's colour (by default the `block` of the end on
// the bridge's `a` side, EC-218), widest first; past `pairs.maxArcs` the most aligned length kept, the rest a density
// strip along each edge. The hit is closed form: at the point's height every ribbon's two sides are one multiply-add
// each, scanned last drawn first. Part of the bridge chunk.

import { defineIssues, issueOf } from '@molstar/core/repr/tree';
import { theme as themeForm } from '@molstar/seqx-core/lang';
import { location } from '@molstar/seqx-core/query';
import type { PairsData } from '@molstar/seqx-core/resolve';
import type { Channel, ChannelName, ThemeDatum } from '@molstar/seqx-core/theme/index';
import { channelLegends } from './channels.js';
import type { RenderConfig } from './config.js';
import { grown, NoColor, Op, PathCommand, type PassName, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { Box, Emphasis, Hit, LayerType, LayoutContext, StoreReader } from './layer-type.js';
import type { BridgePlacement, BridgeSide } from './layout/bridge.js';
import type { PairsLay, PairsPrep, ScratchHit } from './pairs-lay.js';
import { BridgeHooks } from './bridge-hooks.js';
import { type SeqTokens, withAlpha } from './tokens.js';

/** The ribbons' codes, cataloged with the chunk (a composition that draws no bridge says none). */
export const RibbonIssues = defineIssues('seq.render', {
  'seq.render.bridge.unplaced': {
    severity: 'info',
    description: 'A bridge’s pairs have an end on a system neither side shows: those pairs are not drawn',
    hint: 'Show the pairs’ systems on the bridge’s sides (pairs are not projected through maps)',
  },
});

/** What a layer laid out in a bridge reads besides a layout's context: its bridge (_spec/README.md "Render", Bridges). */
export interface BridgeContext {
  /** The bridge's band and sides. */
  readonly placement: BridgePlacement;
  /** Its index in `ResolvedView.bridges` (a hit's `bridge`) and its id (a hit's `section`). */
  readonly index: number;
  readonly id: string;
}

/** A layout's context in a bridge: `box` the band, `section` and `axis` its `a` side's, `bridge` the rest. */
export type BridgeLayoutContext = LayoutContext & { readonly bridge: BridgeContext };

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

/** A pair's end's system: its row's unit under a complex, else the layer's (the pairs layer's `aOf`, `bOf`). */
function aOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.aUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.aSystem) : prep.aSystem;
}
function bOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.bUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.bSystem) : prep.bSystem;
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

/** A ribbon drawn as one curve (its ends both under `minPx`, or `shape: line`), and a stub (its far end not shown). */
const Line = 1;
const Stub = 2;

/** A ribbon's numbers: its `a` side's corners and y, then its `b` side's (`x0a, x1a, ya, x0b, x1b, yb`). */
const R = 6;

export class RibbonsLay {
  key = '';
  /** The bridge's id and index (what a hit names). */
  bridgeId = '';
  bridge = -1;
  tokens: SeqTokens | undefined;
  config: RenderConfig | undefined;
  prep: PairsPrep | undefined;
  dpr = 1;
  colorChannel: Channel | undefined;
  /** The ribbons in drawing order (widest first): corners, the pair, its colour, flags, its data's `a` end's side (1: a). */
  n = 0;
  rib = new Float64Array(6 * 64);
  pair = new Int32Array(64);
  color = new Uint32Array(64);
  flags = new Uint8Array(64);
  orient = new Int32Array(64);
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
  /** A region read's state (its level, its items gate) is kept by this object (the pairs layer's reader's key). */
  readonly side = {} as PairsLay;
  readonly scratch: ScratchHit = { layer: '', section: '', system: '' };
}

// ---- candidates (kept: grown, never shrunk) ----

let Cand = new Float64Array(6 * 256);
let CandPair = new Int32Array(256);
let CandFlags = new Uint8Array(256);
let CandOrient = new Int32Array(256);
let CandBlock = new Int32Array(256);
let CandWidth = new Float64Array(256);
let CandRank = new Float64Array(256);
let Order = new Int32Array(256);
let nc = 0;

function room(n: number) {
  if (CandPair.length >= n) return;
  const m = Math.max(n, 2 * CandPair.length);
  Cand = grown(Cand, R * m);
  CandPair = grown(CandPair, m);
  CandFlags = grown(CandFlags, m);
  CandOrient = grown(CandOrient, m);
  CandBlock = grown(CandBlock, m);
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
  CandFlags[nc] = flags;
  CandOrient[nc] = orient;
  CandBlock[nc] = block;
  CandWidth[nc] = Math.abs(xa1 - xa0) + Math.abs(xb1 - xb0);
  CandRank[nc] = rank;
  nc++;
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
  const na = sa.extentInto(orient > 0 ? aSys : bSys, as, ae, PlacesA, BlocksA, MaxPlaces);
  const nb = sb.extentInto(orient > 0 ? bSys : aSys, bs, be, PlacesB, BlocksB, MaxPlaces);
  if (!na && !nb) return false;
  // corresponding positions: a's first to b's first when the strands' product is 1, to b's last when −1 (GA's twist)
  const product = d.aStrand[k] * d.bStrand[k];
  const ya = sa.edge;
  const yb = sb.edge;
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
  for (let i = 0; i < na; i++)
    for (let j = 0; j < nb; j++) {
      const a0 = PlacesA[2 * i];
      const a1 = PlacesA[2 * i + 1];
      // the b corners in the order corresponding to a's: turned where the run is against
      const turned = product < 0;
      const b0 = turned ? PlacesB[2 * j + 1] : PlacesB[2 * j];
      const b1 = turned ? PlacesB[2 * j] : PlacesB[2 * j + 1];
      Corner[0] = a0;
      Corner[1] = a1;
      Corner[2] = b0;
      Corner[3] = b1;
      const thinA = widen(Corner, 0, min);
      const thinB = widen(Corner, 2, min);
      const flags = line || (thinA && thinB) ? Line : 0;
      push(Corner[0], Corner[1], ya, Corner[2], Corner[3], yb, k, flags, orient, BlocksA[i], rank);
    }
  return true;
}

const Corner = new Float64Array(4);

/** The systems neither side shows, the first few (said once). */
const UnplacedSystems: string[] = [];

/**
 * Lays out a layer's ribbons in a bridge: each pair placed by the systems its sides show (both: the data's `a` end on
 * the bridge's `a` side; `mirror` the swap too), its ribbons ranked past `maxArcs`, the kept ones widest first.
 */
export function layoutRibbons(lay: RibbonsLay, prep: PairsPrep, ctx: BridgeLayoutContext): void {
  const d = prep.data;
  const { placement } = ctx.bridge;
  const sa = placement.a;
  const sb = placement.b;
  lay.band = placement.box;
  lay.edgeA = sa.edge;
  lay.edgeB = sb.edge;
  if (!d || !sa.ready || !sb.ready) return;
  const config = ctx.config;
  const min = config.ribbons.minPx / ctx.dpr;
  const stub = config.ribbons.stubShare;
  const line = (prep.channels.forms.shape as unknown) === 'line';
  nc = 0;
  UnplacedSystems.length = 0;
  let unplaced = 0;
  for (let k = 0; k < d.n; k++) {
    if (!prep.shown[k]) continue;
    const as = aOf(prep, k);
    const bs = bOf(prep, k);
    const fwd = sa.shows(as) && sb.shows(bs);
    const back = sa.shows(bs) && sb.shows(as);
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
  lay.color = grown(lay.color, kept);
  lay.flags = grown(lay.flags, kept);
  lay.orient = grown(lay.orient, kept);
  const datum = ribbonDatum(prep);
  const a = sa.section;
  const domain = a.domains?.block ?? [a.address.slice(a.address.lastIndexOf('/') + 1)];
  lay.colorChannel = prep.channels.get('color' as ChannelName, DefaultColor, ctx, {
    numeric: numericOf(d),
    categories: (f) => (f === 'block' ? domain : f === 'system' ? a.domains?.system : undefined),
  });
  const fill = config.ribbons.fillOpacity;
  for (let oi = 0; oi < kept; oi++) {
    const c = top[oi];
    const i = lay.n++;
    for (let r = 0; r < R; r++) lay.rib[R * i + r] = Cand[R * c + r];
    const k = CandPair[c];
    lay.pair[i] = k;
    lay.flags[i] = CandFlags[c];
    lay.orient[i] = CandOrient[c];
    const b = CandBlock[c];
    const block = b >= 0 ? a.blocks?.[b] : undefined;
    datum.set(k);
    datum.block = block ? block.label : domain[0];
    datum.system = block ? (block.name ?? block.axis.system) : a.axis.system;
    const v = lay.colorChannel.at(datum);
    const color = typeof v === 'number' ? v : lay.tokens!.colors['seq-range'];
    const alpha = ((color & 0xff) / 255) * (prep.weight ? binAlpha(prep, k, config) : fill);
    lay.color[i] = withAlpha(color, alpha);
  }
}

/** A bin's opacity by its aligned length (its weight) on a square-root scale over the bins in view. */
function binAlpha(prep: PairsPrep, k: number, config: RenderConfig): number {
  const w = prep.weight!;
  let max = 0;
  for (let i = 0; i < w.length; i++) if (w[i] > max) max = w[i];
  const { min, max: hi } = config.ribbons.binOpacity;
  return max > 0 ? min + (hi - min) * Math.sqrt(w[k] / max) : hi;
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

// ---- draw ----

/** A ribbon's closed path (two cubics with vertical tangents joined by its ends) into a pass, through `args`. */
export function ribbonPath(out: PassWriter, xy: Float64Array, o: number, close: boolean) {
  const a = out.args;
  out.pass.op(Op.path);
  const p = out.pass.paths;
  p.begin();
  p.commands = grown(p.commands, p.commandCount + 5);
  p.points = grown(p.points, 2 * (p.pointCount + 8));
  const c = p.commands;
  const ya = xy[o + 2];
  const yb = xy[o + 5];
  const y1 = ya + (yb - ya) / 3;
  const y2 = ya + (2 * (yb - ya)) / 3;
  const q = p.points;
  let k = 2 * p.pointCount;
  c[p.commandCount++] = PathCommand.move;
  c[p.commandCount++] = PathCommand.cubic;
  q[k++] = xy[o];
  q[k++] = ya;
  q[k++] = xy[o];
  q[k++] = y1;
  q[k++] = xy[o + 3];
  q[k++] = y2;
  q[k++] = xy[o + 3];
  q[k++] = yb;
  p.pointCount += 4;
  if (close) {
    c[p.commandCount++] = PathCommand.line;
    c[p.commandCount++] = PathCommand.cubic;
    c[p.commandCount++] = PathCommand.close;
    q[k++] = xy[o + 4];
    q[k++] = yb;
    q[k++] = xy[o + 4];
    q[k++] = y2;
    q[k++] = xy[o + 1];
    q[k++] = y1;
    q[k++] = xy[o + 1];
    q[k++] = ya;
    p.pointCount += 4;
  }
  const i = p.n++;
  p.cmdCount[i] = p.commandCount - p.cmdStart[i];
  p.fill[i] = a[4];
  p.stroke[i] = a[5];
  p.strokeWidth[i] = a[6];
  p.dash[i] = a[7];
}

/** A line ribbon's middle curve: its corners collapsed onto the ends' middles (into `Mid`). */
const Mid = new Float64Array(R);
function middle(rib: Float64Array, o: number): Float64Array {
  Mid[0] = Mid[1] = (rib[o] + rib[o + 1]) / 2;
  Mid[2] = rib[o + 2];
  Mid[3] = Mid[4] = (rib[o + 3] + rib[o + 4]) / 2;
  Mid[5] = rib[o + 5];
  return Mid;
}

function clipBand(out: PassWriter, b: Box) {
  const g = out.args;
  g[0] = b.x;
  g[1] = b.y;
  g[2] = b.width;
  g[3] = b.height;
  out.clipArgs();
}

function drawData(lay: RibbonsLay, out: PassWriter) {
  const g = out.args;
  clipBand(out, lay.band);
  for (let i = 0; i < lay.n; i++) {
    const line = lay.flags[i] & Line;
    g[4] = line ? NoColor : lay.color[i];
    g[5] = line ? lay.color[i] | 0xff : NoColor;
    g[6] = line ? 1 : 0;
    g[7] = 0;
    if (line) ribbonPath(out, middle(lay.rib, R * i), 0, false);
    else ribbonPath(out, lay.rib, R * i, true);
  }
  if (lay.densityN && lay.densityMax) {
    const dens = lay.tokens!.colors['seq-density-high'];
    const b = lay.band;
    for (let side = 0; side < 2; side++) {
      const y = side === 0 ? lay.edgeA : lay.edgeB;
      const down = y <= b.y + b.height / 2;
      for (let c = 0; c < lay.densityN; c++) {
        const v = lay.density[side * lay.densityN + c];
        if (!v) continue;
        g[0] = b.x + c / lay.dpr;
        g[1] = down ? y : y - 3;
        g[2] = 1 / lay.dpr;
        g[3] = 3;
        g[4] = ((dens & 0xffffff00) | Math.round(64 + (191 * v) / lay.densityMax)) >>> 0;
        g[5] = NoColor;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
    }
  }
  out.unclip();
}

/** Whether a ribbon's pair has an end in a Location's ranges (a hovered position, the selection). */
function endIn(
  prep: PairsPrep,
  k: number,
  parts: readonly { readonly system: string; readonly ranges: ArrayLike<number> }[],
): boolean {
  const d = prep.data!;
  const as = aOf(prep, k);
  const bs = bOf(prep, k);
  for (let i = 0; i < parts.length; i++) {
    const a = parts[i].system === as;
    const b = parts[i].system === bs;
    if (!a && !b) continue;
    const r = parts[i].ranges;
    for (let j = 0; j + 1 < r.length; j += 2)
      if ((a && d.aStart[k] < r[j + 1] && r[j] < d.aEnd[k]) || (b && d.bStart[k] < r[j + 1] && r[j] < d.bEnd[k]))
        return true;
  }
  return false;
}

/** Whether a ribbon's pair is lit: hovered (its element, or an end at the hovered position), or selected. */
function lit(lay: RibbonsLay, pass: PassName, e: Emphasis, k: number, hovered: number): boolean {
  const prep = lay.prep!;
  if (pass === 'hover') {
    if (k === hovered) return true;
    const h = e.hover;
    if (!h || h.element) return false;
    if (h.position !== undefined) {
      const p = Math.floor(h.position);
      const d = prep.data!;
      if (h.system === aOf(prep, k) && p >= d.aStart[k] && p < d.aEnd[k]) return true;
      if (h.system === bOf(prep, k) && p >= d.bStart[k] && p < d.bEnd[k]) return true;
    }
    const projected = e.projected?.hover;
    return projected !== undefined && endIn(prep, k, projected.parts);
  }
  const keys = e.elements?.keys;
  if (keys) {
    const key = prep.data!.keys[k];
    for (let i = 0; i < keys.length; i++) if (keys[i].owner === prep.owner && keys[i].key === key) return true;
  }
  return e.selection !== undefined && endIn(prep, k, e.selection.parts);
}

function drawEmphasis(lay: RibbonsLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const prep = lay.prep!;
  const h = e.hover;
  let hovered = -1;
  if (pass === 'hover') {
    if (!h && !e.projected?.hover) return;
    if (h?.element) {
      if (h.element.owner !== prep.owner) return;
      hovered = prep.index.get(h.element.key) ?? -1;
      if (hovered < 0) return;
    }
  } else if (!e.elements?.keys.length && !e.selection) return;
  const tokens = lay.tokens!;
  const config = lay.config!;
  const ink = pass === 'hover' ? tokens.colors['seq-hover'] : tokens.colors['msui-accent'];
  const width = pass === 'hover' ? config.outline.hover + 1 : config.outline.selection + 1;
  const g = out.args;
  clipBand(out, lay.band);
  for (let i = 0; i < lay.n; i++) {
    const k = lay.pair[i];
    if (!lit(lay, pass, e, k, hovered)) continue;
    const line = lay.flags[i] & Line;
    // the hovered ribbon at full opacity with the hover's ink as its outline; the lit and selected outlined
    g[4] = pass === 'hover' && k === hovered && !line ? lay.color[i] | 0xff : NoColor;
    g[5] = ink;
    g[6] = width;
    g[7] = 0;
    if (line) ribbonPath(out, middle(lay.rib, R * i), 0, false);
    else ribbonPath(out, lay.rib, R * i, true);
  }
  out.unclip();
}

// ---- hit ----

/**
 * The ribbon under a page point, closed form (_spec/README.md "Render", Bridges): at the point's height `u` each
 * ribbon's sides are `x_a + (x_b − x_a)·S(u)`; the point is in it when its x lies between them (widened to
 * `pairs.hitPx` for a thin one), scanned last drawn first (the narrowest on top). Nothing allocated.
 */
export function hitRibbons(lay: RibbonsLay, at: Float64Array): Hit | undefined {
  const x = at[0];
  const y = at[1];
  const b = lay.band;
  if (!lay.prep || y < b.y || y >= b.y + b.height || x < b.x || x >= b.x + b.width) return undefined;
  const r = lay.config!.pairs.hitPx;
  const rib = lay.rib;
  for (let i = lay.n - 1; i >= 0; i--) {
    const o = R * i;
    const ya = rib[o + 2];
    const yb = rib[o + 5];
    const u = (y - ya) / (yb - ya);
    if (!(u >= 0 && u <= 1)) continue;
    const s = u * u * (3 - 2 * u);
    let lo = rib[o] + (rib[o + 3] - rib[o]) * s;
    let hi = rib[o + 1] + (rib[o + 4] - rib[o + 1]) * s;
    if (lay.flags[i] & Line) lo = hi = (lo + hi) / 2;
    if (lo > hi) {
      const t = lo;
      lo = hi;
      hi = t;
    }
    if (hi - lo < r) {
      const m = (lo + hi) / 2;
      lo = m - r / 2;
      hi = m + r / 2;
    }
    if (x >= lo && x <= hi) return ribbonHit(lay, lay.pair[i]);
  }
  return undefined;
}

/** A ribbon's hit: its pair's element and both ends' extent, its bridge (made once a pair: `elements`, `extents`). */
function ribbonHit(lay: RibbonsLay, k: number): Hit {
  const prep = lay.prep!;
  const d = prep.data!;
  let element = prep.elements[k];
  if (!element) element = prep.elements[k] = { owner: prep.owner, key: d.keys[k] };
  let extent = prep.extents[k];
  const as = aOf(prep, k);
  const bs = bOf(prep, k);
  if (!extent) {
    const parts =
      as === bs
        ? [{ system: as, ranges: [d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k]] }]
        : [
            { system: as, ranges: [d.aStart[k], d.aEnd[k]] },
            { system: bs, ranges: [d.bStart[k], d.bEnd[k]] },
          ];
    extent = prep.extents[k] = location({ parts });
  }
  const h = lay.scratch;
  h.layer = lay.key;
  h.section = lay.bridgeId;
  h.system = as;
  h.element = element;
  h.extent = extent;
  h.position = undefined;
  h.ySystem = undefined;
  h.yPosition = undefined;
  h.target = undefined;
  h.bin = prep.level;
  h.bridge = lay.bridge;
  return h;
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
    // a source's pairs read by region: this slice, on the bridge's `a` side's window at its level (bins, or items once
    // they fit), each pair once (_spec/README.md "Render", Bridges: the per-side reads are LA2's)
    let read = prep;
    const sa = c.bridge.placement.a;
    if (prep.region && !prep.data && sa.ready) {
      const side = { ...ctx, section: sa.section, axis: sa.layout.axis!, place: sa.layout.placement };
      const r = BridgeHooks.region?.(prep, prep.region, side, data, lay.side);
      if (r?.data) read = r;
    }
    lay.prep = read;
    layoutRibbons(lay, read, c);
    return lay;
  },
  draw(lay, pass, out, emphasis) {
    if (!lay.tokens || !lay.prep) return;
    if (pass === 'data') drawData(lay, out);
    else if (emphasis && (pass === 'hover' || pass === 'selection')) drawEmphasis(lay, pass, out, emphasis);
  },
  hit: hitRibbons,
  legend(lay, options) {
    return channelLegends([lay.colorChannel], { ...options, layer: lay.key });
  },
  filtered: (lay) => lay.prep?.filtered ?? 0,
};
