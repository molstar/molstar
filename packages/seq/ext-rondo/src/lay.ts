/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type Hit,
  type LayoutContext,
  type PairsRead,
  PairsReader,
  type PairsSource,
  PlacePoint,
  placesOf,
  type ResolvedSection,
  type StaticModel,
  type StoreReader,
} from '@molstar/seqx-core/extension';
import type { RondoConfig } from './config.js';
import { addArc, addCubic, grow, Polygons, RibbonFrame, ribbonCorners } from './geometry.js';
import { hueColors } from './hues.js';

// The ribbons laid out (_spec/README.md, Stems, Ribbons, Stubs): the reader's pairs (bins at a level, or items under
// the gate: `rondo.maxRibbons`) grouped into stems (a coarser level, `2^j` × the read's, no wider than `rondo.stemPx`
// CSS px at the ring's outer circumference), stem pairs as ribbons by strand, an end with no place a stub by partner,
// what is under a device px or past the bound a stem's fan; then each as a polygon, drawn and hit as one.

/** What the layer reads, made once a resolved layer. */
export interface RondoPrep {
  readonly source: PairsSource;
  readonly model: StaticModel;
  readonly stems: 'level' | 'none';
  readonly bundle: number;
  readonly stubs: 'partner' | 'none';
  /** The section a double-click opens the hit's chromosomes in (`section.include`). */
  readonly opens?: string;
  /** The colours: the 24 hues (`rondo:hues` named by the layer's colour), else the categorical tokens. */
  readonly hues: boolean;
  readonly config: RondoConfig;
}

/** A polygon's kind. */
export const Kind = { trumpet: 0, ribbon: 1, stub: 2, fan: 3 } as const;

/** The layer's layout: its stems, ribbons, stubs and fans, their polygons, its hit grid. */
export class RondoLay {
  key = '';
  section = '';
  /** The section's id (a hit's `section`). */
  sectionId = '';
  owner = '';
  prep: RondoPrep | undefined = undefined;
  readonly reader = new PairsReader();
  read: PairsRead | undefined = undefined;
  /** Why nothing is drawn: `pending` (bins not landed), `layout` (text or grid), `none` (no pairs). */
  empty: '' | 'pending' | 'layout' | 'none' = '';
  dpr = 1;
  readonly frame = new RibbonFrame();
  /** The read's level (1: items) and the stems' (positions a stem). */
  level = 1;
  stemLevel = 1;
  items = false;
  /** The px a count (`k`). */
  k = 0;
  fillOpacity = 0.35;
  /** The colours: the hue a system (by its index in the systems' order), the hover outline, the wash, the accent. */
  hues: Uint32Array = new Uint32Array(8);
  hover = 0;
  wash = 0;
  accent = 0;
  ink = 0;
  font = '';
  /** The systems the ends are on, in the assembly's order (the hue domain, taken before `include`), and their names. */
  systems: readonly string[] = [];
  names: string[] = [];

  // ---- stems ----
  ns = 0;
  stemSystem = new Int32Array(64);
  stemBin = new Float64Array(64);
  stemStart = new Float64Array(64);
  stemEnd = new Float64Array(64);
  stemPlaced = new Uint8Array(64);
  stemDir = new Int8Array(64);
  /** Mouth: its places (t) at the band's inner edge, low and high; its middle. */
  stemT0 = new Float64Array(64);
  stemT1 = new Float64Array(64);
  /** The pairs leaving it (ribbons, stubs, fan), within it, the partner systems it reaches. */
  stemCount = new Float64Array(64);
  stemWithin = new Float64Array(64);
  stemFan = new Float64Array(64);
  stemFanTo = new Int32Array(64);
  stemPartners = new Int32Array(64);
  /** Its neck: its place's half width (in place units at the neck) and its first free place (slices are cut from it). */
  stemNeck = new Float64Array(64);
  stemCut = new Float64Array(64);
  stemPoly = new Int32Array(64);
  stemFanPoly = new Int32Array(64);
  stemKey: string[] = [];

  // ---- ribbons (a stem pair a strand) ----
  nr = 0;
  ribA = new Int32Array(64);
  ribB = new Int32Array(64);
  ribStrand = new Int8Array(64);
  ribCount = new Float64Array(64);
  ribDrawn = new Uint8Array(64);
  ribKind = new Int32Array(64);
  /** Its slices: at A (low, high place) and at B. */
  ribA0 = new Float64Array(64);
  ribA1 = new Float64Array(64);
  ribB0 = new Float64Array(64);
  ribB1 = new Float64Array(64);
  ribPoly = new Int32Array(64);
  ribKey: string[] = [];
  /** The kinds the ribbon's rows carry (`type`), the first one met. */
  kinds: readonly string[] = [];

  // ---- stubs (a stem a partner system with no place) ----
  nb = 0;
  stubStem = new Int32Array(16);
  stubPartner = new Int32Array(16);
  stubCount = new Float64Array(16);
  stubDrawn = new Uint8Array(16);
  stubT0 = new Float64Array(16);
  stubT1 = new Float64Array(16);
  stubPoly = new Int32Array(16);
  stubLabel: string[] = [];
  stubLabelX = new Float64Array(16);
  stubLabelY = new Float64Array(16);
  stubKey: string[] = [];

  // ---- polygons, in draw order, and the hit grid over them ----
  readonly polys = new Polygons();
  polyKind = new Int8Array(64);
  polyRef = new Int32Array(64);
  polyColor = new Uint32Array(64);
  order = new Int32Array(64);
  hits: Hit[] = [];
  byKey = new Map<string, number>();
  gridX = 0;
  gridY = 0;
  gridCols = 0;
  gridRows = 0;
  cell = 16;
  gridStart = new Int32Array(1);
  gridItems = new Int32Array(64);
  /** The hole as a polygon (its wash), and the strongest ribbon. */
  hole = -1;
  strongest = -1;
}

const StemKeys = new Map<number, number>();
const RibbonKeys = new Map<number, number>();
const StubKeys = new Map<number, number>();

/** The stem of an end: its index (made on first meeting), `-1` when its system is not one of the layer's. */
function stemOf(lay: RondoLay, sys: number, pos: number): number {
  const bin = Math.floor(pos / lay.stemLevel);
  const key = sys * 4294967296 + bin;
  let s = StemKeys.get(key);
  if (s !== undefined) return s;
  s = lay.ns++;
  if (s >= lay.stemSystem.length) growStems(lay, s + 1);
  lay.stemSystem[s] = sys;
  lay.stemBin[s] = bin;
  lay.stemCount[s] = 0;
  lay.stemWithin[s] = 0;
  lay.stemFan[s] = 0;
  lay.stemFanTo[s] = 0;
  lay.stemPartners[s] = 0;
  lay.stemPlaced[s] = 0;
  lay.stemPoly[s] = -1;
  lay.stemFanPoly[s] = -1;
  StemKeys.set(key, s);
  return s;
}

function growStems(lay: RondoLay, n: number) {
  lay.stemSystem = grow(lay.stemSystem, n);
  lay.stemBin = grow(lay.stemBin, n);
  lay.stemStart = grow(lay.stemStart, n);
  lay.stemEnd = grow(lay.stemEnd, n);
  lay.stemPlaced = grow(lay.stemPlaced, n);
  lay.stemDir = grow(lay.stemDir, n);
  lay.stemT0 = grow(lay.stemT0, n);
  lay.stemT1 = grow(lay.stemT1, n);
  lay.stemCount = grow(lay.stemCount, n);
  lay.stemWithin = grow(lay.stemWithin, n);
  lay.stemFan = grow(lay.stemFan, n);
  lay.stemFanTo = grow(lay.stemFanTo, n);
  lay.stemPartners = grow(lay.stemPartners, n);
  lay.stemNeck = grow(lay.stemNeck, n);
  lay.stemCut = grow(lay.stemCut, n);
  lay.stemPoly = grow(lay.stemPoly, n);
  lay.stemFanPoly = grow(lay.stemFanPoly, n);
}

function growRibbons(lay: RondoLay, n: number) {
  lay.ribA = grow(lay.ribA, n);
  lay.ribB = grow(lay.ribB, n);
  lay.ribStrand = grow(lay.ribStrand, n);
  lay.ribCount = grow(lay.ribCount, n);
  lay.ribDrawn = grow(lay.ribDrawn, n);
  lay.ribKind = grow(lay.ribKind, n);
  lay.ribA0 = grow(lay.ribA0, n);
  lay.ribA1 = grow(lay.ribA1, n);
  lay.ribB0 = grow(lay.ribB0, n);
  lay.ribB1 = grow(lay.ribB1, n);
  lay.ribPoly = grow(lay.ribPoly, n);
}

function growStubs(lay: RondoLay, n: number) {
  lay.stubStem = grow(lay.stubStem, n);
  lay.stubPartner = grow(lay.stubPartner, n);
  lay.stubCount = grow(lay.stubCount, n);
  lay.stubDrawn = grow(lay.stubDrawn, n);
  lay.stubT0 = grow(lay.stubT0, n);
  lay.stubT1 = grow(lay.stubT1, n);
  lay.stubPoly = grow(lay.stubPoly, n);
  lay.stubLabelX = grow(lay.stubLabelX, n);
  lay.stubLabelY = grow(lay.stubLabelY, n);
}

/** The read's systems: the assembly's sequences (the ends' units index them), else the layer's two. */
function systemsOf(read: PairsRead): readonly string[] {
  if (read.endUnits) return read.endUnits;
  return read.aSystem === read.bSystem ? [read.aSystem] : [read.aSystem, read.bSystem];
}

/**
 * The positions a CSS px at the ring's outer circumference (a row's: the band's), the finest strip's: what the stem
 * level is chosen against, so the stems are the same at any device pixel ratio and in the SVG.
 */
function positionsPerCssPx(ctx: LayoutContext): number {
  const place = ctx.place!;
  const s = place.strips;
  let ppp = Number.POSITIVE_INFINITY;
  for (let i = 0; i < s.n; i++) {
    const ax = s.axes[i];
    let pxPerPos = ax.scale;
    if (s.warp[i] === PolarWarpKind && place.warpOf) {
      const w = place.warpOf(i);
      pxPerPos = (ax.scale * Math.abs(w.sweep) * w.r0) / Math.max(1e-9, w.width);
    }
    if (pxPerPos > 0) ppp = Math.min(ppp, 1 / pxPerPos);
  }
  return ppp;
}

/** A strip's polar warp kind (`WarpKind.polar`). */
const PolarWarpKind = 1;

/** The mouth's place of a section unit: an angle on a ring (the strip's frame x through its warp), a page x on a row. */
function placeOfUnit(ctx: LayoutContext, u: number, out: { t: number; strip: number }): boolean {
  const place = ctx.place!;
  if (!place.place(u, Pt)) return false;
  const s = place.strips;
  const i = Pt.strip;
  out.strip = i;
  if (s.warp[i] === PolarWarpKind && place.warpOf) {
    const x = s.axes[i].xOfUnit((u - s.first[i]) / s.scale[i]);
    out.t = place.warpOf(i).angleAt(x);
  } else out.t = Pt.x;
  return true;
}

const Pt = new PlacePoint();
const At = { t: 0, strip: -1 };

/** The frame: a ring's centre and its band's inner edge, or a row's band; false for another layout. */
function frameOf(lay: RondoLay, ctx: LayoutContext): boolean {
  const place = ctx.place;
  const f = lay.frame;
  if (!place || !place.strips.n) return false;
  const s = place.strips;
  const nb = s.bands;
  const band = ctx.band ?? nb - 1;
  if (s.warp[0] === PolarWarpKind && place.warpOf) {
    const w = place.warpOf(0);
    f.polar = true;
    f.cx = w.cx;
    f.cy = w.cy;
    f.r = w.r0 - (s.bandY[band] + s.bandH[band]);
    f.room = Math.max(1, f.r);
    return true;
  }
  if (ctx.section.layout !== 'tracks') return false;
  f.polar = false;
  f.r = s.bandY[band] + s.dy[0];
  f.room = Math.max(1, s.bandH[band]);
  f.x0 = ctx.box.x;
  f.x1 = ctx.box.x + ctx.box.width;
  return true;
}

/** The section's systems' display direction (a block's own: reversed −1). */
function directionOf(section: ResolvedSection, block: number): number {
  const b = section.blocks?.[block];
  return (b ? b.axis.dir : section.axis.dir) === -1 ? -1 : 1;
}

/** The stems' places: each stem's mouth from its first to its last position shown (its block's direction). */
function placeStems(lay: RondoLay, ctx: LayoutContext, model: StaticModel) {
  const section = ctx.section;
  for (let s = 0; s < lay.ns; s++) {
    const sys = lay.systems[lay.stemSystem[s]];
    const len = model.systems.get(sys)?.length ?? Number.POSITIVE_INFINITY;
    const start = lay.stemBin[s] * lay.stemLevel;
    const end = Math.min(start + lay.stemLevel, len);
    lay.stemStart[s] = start;
    lay.stemEnd[s] = end;
    const a = placesOf(section, sys, start);
    if (!a.length || !placeOfUnit(ctx, a[0].unit, At)) continue;
    const ta = At.t;
    const b = placesOf(section, sys, Math.max(start, end - 1));
    const unitB = b.length ? b[0].unit : a[0].unit;
    const dir = directionOf(section, a[0].block);
    // the last position's far edge: one position on (a scale of 1 a position; a block's own scale is the strip's)
    if (!placeOfUnit(ctx, unitB + 0.999, At)) continue;
    lay.stemPlaced[s] = 1;
    lay.stemDir[s] = dir;
    lay.stemT0[s] = Math.min(ta, At.t);
    lay.stemT1[s] = Math.max(ta, At.t);
  }
}

/** The ribbons' read options (the gate the drawn bound: `rondo.maxRibbons`), kept. */
const Options = { gate: 1024 };

/**
 * Lays out the ribbons of a frame (allocation bounded a frame: the polygons and the hit grid are kept buffers, the
 * stems' and ribbons' keys made as they appear).
 */
export function layoutRondo(lay: RondoLay, prep: RondoPrep, ctx: LayoutContext, data: StoreReader): void {
  const cfg = prep.config;
  lay.prep = prep;
  lay.key = ctx.key;
  lay.section = ctx.section.address;
  lay.sectionId = ctx.section.address.slice(ctx.section.address.lastIndexOf('/') + 1);
  lay.owner = prep.source.owner;
  lay.dpr = ctx.dpr;
  lay.ns = lay.nr = lay.nb = 0;
  lay.polys.reset();
  lay.byKey.clear();
  lay.hits.length = 0;
  lay.hole = -1;
  lay.strongest = -1;
  lay.empty = '';
  const tokens = ctx.tokens;
  lay.fillOpacity = ctx.config.ribbons.fillOpacity;
  lay.hover = tokens.colors['seq-hover'];
  lay.wash = tokens.colors['seq-wash'];
  lay.accent = tokens.colors['msui-accent'];
  lay.ink = tokens.colors['msui-text'];
  lay.font = ctx.config.labelFont;
  if (!frameOf(lay, ctx)) {
    lay.empty = 'layout';
    return;
  }
  Options.gate = cfg.maxRibbons;
  const read = lay.reader.read(prep.source, ctx, data, Options);
  lay.read = read;
  const d = read?.data;
  if (!read || !d) {
    lay.empty = ctx.settled ? 'none' : 'pending';
    return;
  }
  lay.systems = systemsOf(read);
  lay.names = lay.systems.map((s) => prep.model.systems.get(s)?.name ?? s);
  lay.kinds = d.kinds;
  lay.items = read.level === undefined;
  lay.level = read.level ?? 1;
  // the stem level: the largest 2^j × the read's level no wider than `stemPx` CSS px; `stems: none`, the read's own
  const target = cfg.stemPx * positionsPerCssPx(ctx);
  let level = lay.level;
  if (prep.stems === 'level') while (level * 2 <= target) level *= 2;
  lay.stemLevel = level;
  // the hues by the systems' order, the dark theme's derived
  const n = lay.systems.length;
  if (prep.hues) lay.hues = hueColors(Math.max(n, 1), tokens.theme === 'dark');
  else {
    lay.hues = new Uint32Array(8);
    for (let i = 0; i < 8; i++) lay.hues[i] = tokens.colors[`seq-cat-${i}` as 'seq-cat-0'];
  }
  gather(lay, read, d.aUnit, d.bUnit);
  placeStems(lay, ctx, prep.model);
  aggregate(lay, read);
  if (!lay.nr && !lay.nb && !lay.ns) {
    // nothing read yet while the view waits for its data (EC-169: pending is said, never an empty hole), else none
    lay.empty = d.n === 0 && !ctx.settled ? 'pending' : 'none';
    return;
  }
  widths(lay, cfg, ctx.config.ribbons.minPx);
  necks(lay, cfg);
  polygons(lay, prep, cfg, ctx);
  hitGrid(lay, cfg.hitCellPx);
}

/** Row ends into stems (each row once), counted: `rows` its stems by row. */
let RowA = new Int32Array(256);
let RowB = new Int32Array(256);

function gather(lay: RondoLay, read: PairsRead, aUnit: Int32Array | undefined, bUnit: Int32Array | undefined) {
  const d = read.data!;
  StemKeys.clear();
  const aSys = aUnit ? -1 : 0;
  const bSys = bUnit ? -1 : lay.systems.length > 1 ? 1 : 0;
  RowA = grow(RowA, d.n);
  RowB = grow(RowB, d.n);
  for (let k = 0; k < d.n; k++) {
    if (!read.shown[k]) {
      RowA[k] = -1;
      continue;
    }
    const sa = aUnit ? aUnit[k] : aSys;
    const sb = bUnit ? bUnit[k] : bSys;
    RowA[k] = stemOf(lay, sa, d.aStart[k]);
    RowB[k] = stemOf(lay, sb, d.bStart[k]);
  }
}

/** Ribbons by stem pair and strand, stubs by stem and partner, within counts; stems' counts the pairs leaving them. */
function aggregate(lay: RondoLay, read: PairsRead) {
  const d = read.data!;
  const w = read.weight;
  RibbonKeys.clear();
  StubKeys.clear();
  for (let k = 0; k < d.n; k++) {
    const a = RowA[k];
    if (a < 0) continue;
    const b = RowB[k];
    const c = w ? w[k] : 1;
    const pa = lay.stemPlaced[a];
    const pb = lay.stemPlaced[b];
    if (!pa && !pb) continue;
    if (a === b) {
      lay.stemWithin[a] += c;
      continue;
    }
    if (pa && pb) {
      // a stem pair a strand: the stem earlier in the systems' order its `a` (its hue), the strand the rows' relative
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      const first = lay.stemSystem[a] < lay.stemSystem[b] || (lay.stemSystem[a] === lay.stemSystem[b] && a < b) ? a : b;
      const strand = d.aStrand[k] * d.bStrand[k] < 0 ? -1 : 1;
      const key = (lo * 1048576 + hi) * 2 + (strand < 0 ? 1 : 0);
      let r = RibbonKeys.get(key);
      if (r === undefined) {
        r = lay.nr++;
        if (r >= lay.ribA.length) growRibbons(lay, r + 1);
        lay.ribA[r] = first;
        lay.ribB[r] = first === a ? b : a;
        lay.ribStrand[r] = strand;
        lay.ribCount[r] = 0;
        lay.ribKind[r] = d.kind[k];
        RibbonKeys.set(key, r);
      }
      lay.ribCount[r] += c;
      lay.stemCount[a] += c;
      lay.stemCount[b] += c;
      continue;
    }
    // a stub: the placed end's stem toward the partner's system
    const near = pa ? a : b;
    const partner = lay.stemSystem[pa ? b : a];
    const key = near * 65536 + partner;
    let s = StubKeys.get(key);
    if (s === undefined) {
      s = lay.nb++;
      if (s >= lay.stubStem.length) growStubs(lay, s + 1);
      lay.stubStem[s] = near;
      lay.stubPartner[s] = partner;
      lay.stubCount[s] = 0;
      StubKeys.set(key, s);
    }
    lay.stubCount[s] += c;
    lay.stemCount[near] += c;
  }
}

/** Ribbons by count, strongest first (ties: the wider span, then the key): the drawn bound's order. */
let Rank = new Int32Array(64);

/**
 * Widths (_spec/README.md, Ribbons, Width): `k` = `maxWidthPx` / the largest stem count in view; a ribbon at least
 * `ribbons.minPx` device px wide drawn, within the strongest `maxRibbons`; the rest in their stems' fans (exact sums).
 */
function widths(lay: RondoLay, cfg: RondoConfig, minPx: number) {
  let max = 0;
  for (let s = 0; s < lay.ns; s++) if (lay.stemPlaced[s] && lay.stemCount[s] > max) max = lay.stemCount[s];
  lay.k = max > 0 ? cfg.maxWidthPx / max : 0;
  const floor = minPx / lay.dpr;
  Rank = grow(Rank, lay.nr);
  for (let r = 0; r < lay.nr; r++) Rank[r] = r;
  const span = (r: number) => Math.abs(lay.stemBin[lay.ribA[r]] - lay.stemBin[lay.ribB[r]]);
  const ranked = Rank.subarray(0, lay.nr).sort(
    (x, y) => lay.ribCount[y] - lay.ribCount[x] || span(y) - span(x) || x - y,
  );
  let drawn = 0;
  for (let i = 0; i < ranked.length; i++) {
    const r = ranked[i];
    const ok = drawn < cfg.maxRibbons && lay.k * lay.ribCount[r] >= floor;
    lay.ribDrawn[r] = ok ? 1 : 0;
    if (ok) {
      drawn++;
      if (lay.strongest < 0) lay.strongest = r;
      continue;
    }
    lay.stemFan[lay.ribA[r]] += lay.ribCount[r];
    lay.stemFan[lay.ribB[r]] += lay.ribCount[r];
    lay.stemFanTo[lay.ribA[r]]++;
    lay.stemFanTo[lay.ribB[r]]++;
  }
  for (let s = 0; s < lay.nb; s++) {
    const ok = lay.k * lay.stubCount[s] >= floor;
    lay.stubDrawn[s] = ok ? 1 : 0;
    if (!ok) {
      lay.stemFan[lay.stubStem[s]] += lay.stubCount[s];
      lay.stemFanTo[lay.stubStem[s]]++;
    }
  }
}

/** Each neck's members: its ribbons by their partners' places clockwise from it, its stubs by partner, its fan last. */
let Members = new Int32Array(64);
let MemberStart = new Int32Array(64);

/**
 * The necks (_spec/README.md, Ribbons, Order in a neck): `k ×` its stem's count wide (at least a device px as a whole),
 * centred under its mouth at `stemDepth` of the hole inward, its slices cut in the fixed order.
 */
function necks(lay: RondoLay, cfg: RondoConfig) {
  const f = lay.frame;
  const dn = cfg.stemDepth * f.room;
  const perPx = f.perPx(dn);
  for (let s = 0; s < lay.ns; s++) {
    if (!lay.stemPlaced[s]) continue;
    const w = Math.max(lay.k * lay.stemCount[s], 1 / lay.dpr);
    lay.stemNeck[s] = (w * perPx) / 2;
    const mid = (lay.stemT0[s] + lay.stemT1[s]) / 2;
    lay.stemCut[s] = mid - lay.stemNeck[s];
  }
  // members by stem (CSR): drawn ribbons twice (each end), drawn stubs once
  MemberStart = grow(MemberStart, lay.ns + 1);
  MemberStart.fill(0, 0, lay.ns + 1);
  for (let r = 0; r < lay.nr; r++)
    if (lay.ribDrawn[r]) {
      MemberStart[lay.ribA[r] + 1]++;
      MemberStart[lay.ribB[r] + 1]++;
    }
  for (let b = 0; b < lay.nb; b++) if (lay.stubDrawn[b]) MemberStart[lay.stubStem[b] + 1]++;
  for (let s = 0; s < lay.ns; s++) MemberStart[s + 1] += MemberStart[s];
  const fill = new Int32Array(lay.ns);
  Members = grow(Members, MemberStart[lay.ns]);
  // a member: a ribbon r as 2r (at a) or 2r + 1 (at b); a stub b as −(b + 1)
  for (let r = 0; r < lay.nr; r++)
    if (lay.ribDrawn[r]) {
      const a = lay.ribA[r];
      const b = lay.ribB[r];
      Members[MemberStart[a] + fill[a]++] = 2 * r;
      Members[MemberStart[b] + fill[b]++] = 2 * r + 1;
    }
  for (let b = 0; b < lay.nb; b++)
    if (lay.stubDrawn[b]) {
      const s = lay.stubStem[b];
      Members[MemberStart[s] + fill[s]++] = -(b + 1);
    }
  const TwoPi = 2 * Math.PI;
  const midOf = (s: number) => (lay.stemT0[s] + lay.stemT1[s]) / 2;
  for (let s = 0; s < lay.ns; s++) {
    const m = Members.subarray(MemberStart[s], MemberStart[s + 1]);
    if (!m.length) continue;
    const at = midOf(s);
    // a partner's place clockwise from the stem (a row: its x); the nearest clockwise last, so ribbons leave uncrossed
    const delta = (x: number) => {
      if (x < 0) return Number.NEGATIVE_INFINITY;
      const r = x >> 1;
      const other = x & 1 ? lay.ribA[r] : lay.ribB[r];
      const t = midOf(other);
      return f.polar ? -(((t - at) % TwoPi) + TwoPi) % TwoPi : t;
    };
    const order = (x: number) => (x < 0 ? lay.stubPartner[-x - 1] : 0);
    m.sort((x, y) => {
      const dx = delta(x);
      const dy = delta(y);
      if (dx !== dy) return dx > dy ? 1 : -1;
      return order(x) - order(y) || x - y;
    });
    // ribbons first (by place), then stubs (by partner): stubs sort last as −∞ … flipped here
    let ribbons = 0;
    for (let i = 0; i < m.length; i++) if (m[i] >= 0) ribbons++;
    const stubs = m.slice(0, m.length - ribbons);
    m.copyWithin(0, m.length - ribbons);
    m.set(stubs, ribbons);
    for (let i = 0; i < m.length; i++) {
      const x = m[i];
      const count = x < 0 ? lay.stubCount[-x - 1] : lay.ribCount[x >> 1];
      const t0 = lay.stemCut[s];
      const t1 = t0 + lay.k * count * perPx;
      lay.stemCut[s] = t1;
      if (x < 0) {
        lay.stubT0[-x - 1] = t0;
        lay.stubT1[-x - 1] = t1;
      } else if (x & 1) {
        lay.ribB0[x >> 1] = t0;
        lay.ribB1[x >> 1] = t1;
      } else {
        lay.ribA0[x >> 1] = t0;
        lay.ribA1[x >> 1] = t1;
      }
    }
  }
}

/** A polygon made: its kind, what it stands for, its colour. */
function made(lay: RondoLay, kind: number, ref: number, color: number): number {
  const i = lay.polys.end();
  if (i >= lay.polyKind.length) {
    lay.polyKind = grow(lay.polyKind, i + 1);
    lay.polyRef = grow(lay.polyRef, i + 1);
    lay.polyColor = grow(lay.polyColor, i + 1);
  }
  lay.polyKind[i] = kind;
  lay.polyRef[i] = ref;
  lay.polyColor[i] = color;
  return i;
}

const alpha = (c: number, a: number) => ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;

/** A system's hue (its index in the systems' order: the domain taken before `include`). */
export function hueOf(lay: RondoLay, sys: number): number {
  return lay.hues[sys % lay.hues.length];
}

const C4 = new Float64Array(4);

/** The polygons in draw order: trumpets, ribbons widest first (the narrow on top, hittable), stubs, fans. */
function polygons(lay: RondoLay, prep: RondoPrep, cfg: RondoConfig, ctx: LayoutContext) {
  const f = lay.frame;
  const P = lay.polys;
  const seg = Math.max(2, cfg.flatten);
  const dn = cfg.stemDepth * f.room;
  const op = lay.fillOpacity;
  // trumpets (a stem's mouth to its neck: two cubics flaring), stems: `none` draws none
  for (let s = 0; s < lay.ns; s++) {
    if (!lay.stemPlaced[s] || lay.stemCount[s] <= 0) continue;
    const t0 = lay.stemT0[s];
    const t1 = lay.stemT1[s];
    const mid = (t0 + t1) / 2;
    const n0 = mid - lay.stemNeck[s];
    const n1 = mid + lay.stemNeck[s];
    if (prep.stems === 'none') {
      lay.stemPoly[s] = -1;
      continue;
    }
    P.begin();
    f.point(t0, 0);
    P.add(f.px, f.py);
    addArc(f, P, t0, t1, 0);
    addCubic(f, P, seg, t1, 0, t1, dn / 2, n1, dn / 2, n1, dn);
    addArc(f, P, n1, n0, dn);
    addCubic(f, P, seg, n0, dn, n0, dn / 2, t0, dn / 2, t0, 0);
    lay.stemPoly[s] = made(lay, Kind.trumpet, s, alpha(hueOf(lay, lay.stemSystem[s]), op));
  }
  // ribbons widest first
  const ribs: number[] = [];
  for (let r = 0; r < lay.nr; r++) if (lay.ribDrawn[r]) ribs.push(r);
  ribs.sort((x, y) => lay.ribCount[y] - lay.ribCount[x] || x - y);
  const bundle = prep.bundle;
  for (const r of ribs) {
    const a = lay.ribA[r];
    const b = lay.ribB[r];
    // uncrossed: A's low edge to B's high, A's high to B's low; crossed: low to low, high to high
    ribbonCorners(
      lay.ribA0[r],
      lay.ribA1[r],
      lay.stemDir[a],
      lay.ribB0[r],
      lay.ribB1[r],
      lay.stemDir[b],
      lay.ribStrand[r],
      C4,
    );
    const dc = controlDepth(f, dn, bundle, C4[0], C4[1]);
    const dc2 = controlDepth(f, dn, bundle, C4[2], C4[3]);
    P.begin();
    f.point(C4[0], dn);
    P.add(f.px, f.py);
    addCubic(f, P, seg, C4[0], dn, C4[0], dc, C4[1], dc, C4[1], dn);
    f.point(C4[2], dn);
    P.add(f.px, f.py);
    addCubic(f, P, seg, C4[2], dn, C4[2], dc2, C4[3], dc2, C4[3], dn);
    lay.ribPoly[r] = made(lay, Kind.ribbon, r, alpha(hueOf(lay, lay.stemSystem[a]), op));
  }
  // stubs: a taper from the neck toward the centre, in the partner's hue, labelled at its tip where it has room
  const tipD = dn + cfg.stubShare * f.room;
  let lastLabel = Number.NEGATIVE_INFINITY;
  for (let s = 0; s < lay.nb; s++) {
    lay.stubPoly[s] = -1;
    lay.stubLabel[s] = '';
    if (!lay.stubDrawn[s] || prep.stubs === 'none') continue;
    const t0 = lay.stubT0[s];
    const t1 = lay.stubT1[s];
    const mid = (t0 + t1) / 2;
    P.begin();
    f.point(t0, dn);
    P.add(f.px, f.py);
    f.point(t1, dn);
    P.add(f.px, f.py);
    f.point(mid, tipD);
    P.add(f.px, f.py);
    lay.stubPoly[s] = made(lay, Kind.stub, s, alpha(hueOf(lay, lay.stubPartner[s]), op));
    // its label: where the tip is `stubLabelPx` from the last labelled one along the tips' circle
    const along = mid / f.perPx(tipD);
    if (along - lastLabel >= cfg.stubLabelPx) {
      lastLabel = along;
      f.point(mid, tipD + 4);
      lay.stubLabel[s] = `${lay.names[lay.stubPartner[s]]} · ${Math.round(lay.stubCount[s]).toLocaleString('en-US')}`;
      lay.stubLabelX[s] = f.px;
      lay.stubLabelY[s] = f.py;
    }
  }
  // fans: a stem's weak links, one taper of their summed width, last in its neck
  for (let s = 0; s < lay.ns; s++) {
    lay.stemFanPoly[s] = -1;
    if (!lay.stemPlaced[s] || lay.stemFan[s] <= 0) continue;
    const t0 = lay.stemCut[s];
    const t1 = t0 + Math.max(lay.k * lay.stemFan[s], 0.5 / lay.dpr) * f.perPx(dn);
    P.begin();
    f.point(t0, dn);
    P.add(f.px, f.py);
    f.point(t1, dn);
    P.add(f.px, f.py);
    f.point((t0 + t1) / 2, dn + (cfg.stubShare / 2) * f.room);
    P.add(f.px, f.py);
    lay.stemFanPoly[s] = made(lay, Kind.fan, s, alpha(hueOf(lay, lay.stemSystem[s]), op * 0.6));
  }
  // the hole (a ring's) or the band (a row's): the wash's polygon, never hit
  P.begin();
  if (f.polar) addArc(f, P, 0, 2 * Math.PI, dn);
  else {
    P.add(f.x0, f.r);
    P.add(f.x1, f.r);
    P.add(f.x1, f.r + f.room);
    P.add(f.x0, f.r + f.room);
  }
  lay.hole = made(lay, -1, -1, 0);
  hitsOf(lay);
}

/** A ribbon's controls' depth: `(1 − bundle)` of the neck's radius from the centre on a ring; on a row by its span. */
function controlDepth(f: RibbonFrame, dn: number, bundle: number, ta: number, tb: number): number {
  if (f.polar) return f.r - (1 - bundle) * (f.r - dn);
  const reach = Math.min(1, Math.abs(tb - ta) / Math.max(1, f.x1 - f.x0));
  return dn + (f.room - dn) * Math.max(0.15, reach) * (0.5 + bundle / 2);
}

/** Each polygon's hit and its key's index (made at layout: a hit allocates nothing). */
function hitsOf(lay: RondoLay) {
  const owner = lay.owner;
  const section = lay.sectionId;
  const P = lay.polys;
  lay.order = grow(lay.order, P.n);
  for (let i = 0; i < P.n; i++) {
    lay.order[i] = i;
    const kind = lay.polyKind[i];
    const ref = lay.polyRef[i];
    let key = '';
    let stem = -1;
    let count = 0;
    if (kind === Kind.trumpet || kind === Kind.fan) {
      stem = ref;
      key = `stem:${lay.systems[lay.stemSystem[ref]]}:${lay.stemBin[ref]}`;
      count = kind === Kind.fan ? lay.stemFan[ref] : lay.stemCount[ref];
      lay.stemKey[ref] = key;
    } else if (kind === Kind.ribbon) {
      stem = lay.ribA[ref];
      const b = lay.ribB[ref];
      key = `ribbon:${lay.systems[lay.stemSystem[stem]]}:${lay.stemBin[stem]}|${lay.systems[lay.stemSystem[b]]}:${lay.stemBin[b]}|${lay.ribStrand[ref]}`;
      count = lay.ribCount[ref];
      lay.ribKey[ref] = key;
    } else if (kind === Kind.stub) {
      stem = lay.stubStem[ref];
      key = `stub:${lay.systems[lay.stemSystem[stem]]}:${lay.stemBin[stem]}|${lay.systems[lay.stubPartner[ref]]}`;
      count = lay.stubCount[ref];
      lay.stubKey[ref] = key;
    } else {
      lay.hits.push(undefined as never);
      continue;
    }
    // a fan is hit as its stem (its key the stem's): the first made (the trumpet) is the key's
    if (!lay.byKey.has(key)) lay.byKey.set(key, i);
    const sys = lay.systems[lay.stemSystem[stem]];
    lay.hits.push({
      layer: lay.key,
      section,
      system: sys,
      element: { owner, key },
      position: (lay.stemStart[stem] + lay.stemEnd[stem]) / 2,
      bin: lay.stemLevel,
      count,
    });
  }
}

/** A uniform grid over the polygons' boxes, each cell listing the polygons over it (`rondo.hitCellPx`). */
function hitGrid(lay: RondoLay, cellPx: number) {
  const P = lay.polys;
  let x0 = Infinity;
  let y0 = Infinity;
  let x1 = -Infinity;
  let y1 = -Infinity;
  for (let i = 0; i < P.n; i++) {
    if (i === lay.hole) continue;
    x0 = Math.min(x0, P.box[4 * i]);
    y0 = Math.min(y0, P.box[4 * i + 1]);
    x1 = Math.max(x1, P.box[4 * i + 2]);
    y1 = Math.max(y1, P.box[4 * i + 3]);
  }
  const c = Math.max(1, cellPx);
  lay.cell = c;
  if (!(x1 >= x0)) {
    lay.gridCols = lay.gridRows = 0;
    return;
  }
  lay.gridX = x0;
  lay.gridY = y0;
  lay.gridCols = Math.max(1, Math.ceil((x1 - x0) / c) + 1);
  lay.gridRows = Math.max(1, Math.ceil((y1 - y0) / c) + 1);
  const cells = lay.gridCols * lay.gridRows;
  lay.gridStart = grow(lay.gridStart, cells + 1);
  lay.gridStart.fill(0, 0, cells + 1);
  const span = (i: number, add: (cell: number) => void) => {
    const cx0 = Math.floor((P.box[4 * i] - x0) / c);
    const cy0 = Math.floor((P.box[4 * i + 1] - y0) / c);
    const cx1 = Math.floor((P.box[4 * i + 2] - x0) / c);
    const cy1 = Math.floor((P.box[4 * i + 3] - y0) / c);
    for (let y = cy0; y <= cy1; y++) for (let x = cx0; x <= cx1; x++) add(y * lay.gridCols + x);
  };
  for (let i = 0; i < P.n; i++) if (i !== lay.hole) span(i, (k) => lay.gridStart[k + 1]++);
  for (let k = 0; k < cells; k++) lay.gridStart[k + 1] += lay.gridStart[k];
  lay.gridItems = grow(lay.gridItems, lay.gridStart[cells]);
  const fill = new Int32Array(cells);
  // topmost first: the last drawn first in each cell
  for (let i = P.n - 1; i >= 0; i--)
    if (i !== lay.hole) span(i, (k) => (lay.gridItems[lay.gridStart[k] + fill[k]++] = i));
}

/** The polygon under a page point, topmost first, or −1 (through the grid: a few candidates; nothing allocated). */
export function polygonAt(lay: RondoLay, x: number, y: number): number {
  if (!lay.gridCols) return -1;
  const cx = Math.floor((x - lay.gridX) / lay.cell);
  const cy = Math.floor((y - lay.gridY) / lay.cell);
  if (cx < 0 || cy < 0 || cx >= lay.gridCols || cy >= lay.gridRows) return -1;
  const k = cy * lay.gridCols + cx;
  for (let j = lay.gridStart[k]; j < lay.gridStart[k + 1]; j++) {
    const i = lay.gridItems[j];
    if (lay.polys.contains(i, x, y)) return i;
  }
  return -1;
}

/** The stem showing a position of a system, or −1. */
export function stemAt(lay: RondoLay, system: string, position: number): number {
  for (let s = 0; s < lay.ns; s++)
    if (
      lay.stemPlaced[s] &&
      lay.systems[lay.stemSystem[s]] === system &&
      position >= lay.stemStart[s] &&
      position < lay.stemEnd[s]
    )
      return s;
  return -1;
}
