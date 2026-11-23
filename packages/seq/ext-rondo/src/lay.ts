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
  type PairsReadOptions,
  type PairsSource,
  PlacePoint,
  placesOf,
  type ResolvedSection,
  type StaticModel,
  type StoreReader,
} from '@molstar/seqx-core/extension';
import { ColorMode, classOf, clearClasses, colorsOf, peersOf, rowWeights, weightsOf } from './colour.js';
import { composedConfig, DefaultRondoConfig, type RondoConfig } from './config.js';
import { grow, Polygons, RibbonFrame } from './geometry.js';
import { hueColors } from './hues.js';
import { holeWords } from './hole.js';
import { hitGrid, polygons } from './shapes.js';

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
  /** Whether the layer writes its colour (its `color` channel read: a fixed colour, a field's); else the block hue. */
  readonly colored: boolean;
  /** What a ribbon is wide by: `count` (the links), or a numeric field summed (`q`). */
  readonly width: string;
  /** The layer's name in its card and legend (its kind's: `rondoRing`'s `kinds`); empty when none. */
  readonly label: string;
  /**
   * Its widths' scale: `layer` (its own `k`, the default) or `section` (one `k` over the section's layers of
   * `section` and its `width`, so widths compare across kinds).
   */
  readonly scale: 'layer' | 'section';
  /** Whether its data names a filter and its links are read by region (their bins cannot be filtered). */
  readonly regionFilter: boolean;
  readonly config: RondoConfig;
}

export { ColorMode } from './colour.js';

export { Kind } from './shapes.js';

/** A compiled colour channel (the layer's kept channels' `get`). */
type Channel = ReturnType<PairsSource['channels']['get']>;

/** The layer's layout: its stems, ribbons, stubs and fans, their polygons, its hit grid. */
export class RondoLay {
  key = '';
  section = '';
  /** The section's id (a hit's `section`). */
  sectionId = '';
  owner = '';
  prep: RondoPrep | undefined = undefined;
  /** The settings it was laid out with (the descriptor's, the composition's over them). */
  cfg: RondoConfig = DefaultRondoConfig;
  readonly reader = new PairsReader();
  read: PairsRead | undefined = undefined;
  /** Why nothing is drawn: `pending` (bins not landed), `layout` (text or grid), `none` (no pairs), `failed` (its source). */
  empty: '' | 'pending' | 'layout' | 'none' | 'failed' | 'shared' = '';
  /** What the hole says while it is empty for a reason (`Loading links…`, a failure's issue), its measured width. */
  holeText = '';
  holeW = 0;
  /**
   * Its bins are read by region under a filter they cannot be read against (`rondo.ribbons.filter`); `shared`: another
   * layer over the same links draws them (the bins drawn once, never once a layer).
   */
  regionBins = false;
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
  /** A folded stub's colour (no one partner's hue). */
  muted = 0;
  /** How it is coloured (`ColorMode`), its colour channel when it reads one. */
  mode: number = ColorMode.system;
  channel: Channel | undefined = undefined;
  /**
   * The layer's one colour (`0xRRGGBBAA`), stems', stubs' and fans' and every ribbon's but a field's at the item level:
   * the fixed colour; a field colour's value when every link read shares it, else the ranges' ink. −1: the blocks' hues.
   */
  one = -1;
  /** A field colour read link by link (items, its values mixed): each ribbon its links' colour. */
  byRow = false;
  /** Its bins carry not the field it is coloured by: one colour (the ink), said, and the legend says so. */
  inkBins = false;
  /** More colours than a stem pair's key holds (`MaxColorClasses`): some ribbons merged, said. */
  classOverflow = false;
  /** Each row's colour, and the read and channel it was scanned for (kept: a warm layout scans none). */
  rowColor = new Uint32Array(64);
  colorRead: unknown = undefined;
  colorKept: unknown = undefined;
  colorFirst = -1;
  colorMixed = false;
  colorReadable = false;
  /** The layer's `opacity` (its channel's uniform value), a factor of `ribbons.fillOpacity`. */
  opacity = 1;
  /** Its name (`label`) and what it is wide by (`count`, or a field summed). */
  label = '';
  width = 'count';
  /** Its own widest stem's measure, and whether it shares a scale with its section (`scale: section`). */
  ownMax = 0;
  shared = false;
  /** The lays of its section (itself among them): a bin's card counts each layer's links there. */
  peers: readonly RondoLay[] = [];
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
  /**
   * Links, apart from the width measure (`stemCount`, `stubCount`, `stemWithin`, `stemFan`: a field's sum with `width`
   * a field): what folding and the cards count.
   */
  stemLinks = new Float64Array(64);
  withinLinks = new Float64Array(64);
  fanLinks = new Float64Array(64);
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
  /** Its links (its count when wide by `count`), its colour (a field colour by link), its first row (its card's fields). */
  ribLinks = new Float64Array(64);
  ribColor = new Uint32Array(64);
  ribRow = new Int32Array(64);
  /** Each row's ribbon (−1: none), so a card lists the fields a ribbon's links share. */
  rowRib = new Int32Array(64);
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
  stubLinks = new Float64Array(16);
  stubDrawn = new Uint8Array(16);
  stubT0 = new Float64Array(16);
  stubT1 = new Float64Array(16);
  stubPoly = new Int32Array(16);
  stubLabel: string[] = [];
  stubLabelX = new Float64Array(16);
  stubLabelY = new Float64Array(16);
  /** A label's measured width, CSS px. */
  stubLabelW = new Float64Array(16);
  stubKey: string[] = [];
  /**
   * Folded stubs (`rondo.stubMin`): a stub of partner −1 is a stem's small stubs as one; its partners and counts are the
   * entries `e < nf` with `foldOwner[e]` the stub.
   */
  nf = 0;
  foldOwner = new Int32Array(16);
  foldPartner = new Int32Array(16);
  foldCount = new Float64Array(16);
  foldLinks = new Float64Array(16);
  /** The ring shows one chromosome or a pair: the stems drawn lighter (`rondo.focusStems`). */
  focus = false;

  // ---- polygons, in draw order, and the hit grid over them ----
  readonly polys = new Polygons();
  polyKind = new Int8Array(64);
  polyRef = new Int32Array(64);
  polyColor = new Uint32Array(64);
  /** A thin polygon's hit slop past its outline, CSS px (a ribbon or stub under `2 × rondo.hitSlopPx` wide); else 0. */
  polySlop = new Float64Array(64);
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
  lay.stemLinks[s] = 0;
  lay.withinLinks[s] = 0;
  lay.fanLinks[s] = 0;
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
  lay.stemLinks = grow(lay.stemLinks, n);
  lay.withinLinks = grow(lay.withinLinks, n);
  lay.fanLinks = grow(lay.fanLinks, n);
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
  lay.ribLinks = grow(lay.ribLinks, n);
  lay.ribColor = grow(lay.ribColor, n);
  lay.ribRow = grow(lay.ribRow, n);
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
  lay.stubLinks = grow(lay.stubLinks, n);
  lay.stubDrawn = grow(lay.stubDrawn, n);
  lay.stubT0 = grow(lay.stubT0, n);
  lay.stubT1 = grow(lay.stubT1, n);
  lay.stubPoly = grow(lay.stubPoly, n);
  lay.stubLabelX = grow(lay.stubLabelX, n);
  lay.stubLabelY = grow(lay.stubLabelY, n);
  lay.stubLabelW = grow(lay.stubLabelW, n);
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
    // the hole's rim: the innermost band's inner edge, whichever band holds the layer (kinds of links in tracks of
    // their own share one rim)
    f.r = w.r0 - (s.bandY[nb - 1] + s.bandH[nb - 1]);
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

/** The ribbons' read options (the gate the drawn bound: `rondo.maxRibbons`; a field summed when wide by one), kept. */
const Options = { gate: 1024 };
const OptionsSum = { gate: 1024, sum: '' };

/**
 * Lays out the ribbons of a frame (allocation bounded a frame: the polygons and the hit grid are kept buffers, the
 * stems' and ribbons' keys made as they appear).
 */
export function layoutRondo(lay: RondoLay, prep: RondoPrep, ctx: LayoutContext, data: StoreReader): void {
  // the descriptor's numbers, the composition's `rondo.*` settings over them
  const cfg = composedConfig(prep.config, ctx.config);
  lay.cfg = cfg;
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
  lay.holeText = '';
  lay.regionBins = false;
  const tokens = ctx.tokens;
  lay.fillOpacity = ctx.config.ribbons.fillOpacity;
  lay.hover = tokens.colors['seq-hover'];
  lay.wash = tokens.colors['seq-wash'];
  lay.accent = tokens.colors['msui-accent'];
  lay.ink = tokens.colors['msui-text'];
  lay.muted = tokens.colors['msui-text-muted'];
  lay.font = ctx.config.labelFont;
  if (!frameOf(lay, ctx)) {
    lay.empty = 'layout';
    return;
  }
  lay.label = prep.label;
  lay.width = prep.width;
  peersOf(lay, ctx);
  let options: PairsReadOptions = Options;
  if (prep.width !== 'count') {
    OptionsSum.gate = cfg.maxRibbons;
    OptionsSum.sum = prep.width;
    options = OptionsSum;
  } else Options.gate = cfg.maxRibbons;
  const read = lay.reader.read(prep.source, ctx, data, options);
  lay.read = read;
  const d = read?.data;
  if (!read || !d) {
    // none only once settled with no bins' tile awaited (a frame between a source's open and its tiles' wants is not)
    lay.empty = lay.reader.failed ? 'failed' : ctx.settled && !lay.reader.waiting ? 'none' : 'pending';
    holeWords(lay, prep, ctx);
    return;
  }
  lay.systems = systemsOf(read);
  lay.names = lay.systems.map((s) => prep.model.systems.get(s)?.name ?? s);
  lay.kinds = d.kinds;
  lay.items = read.level === undefined;
  // bins read by region carry no fields: a filter is not read on them, so the first layer over the links draws them
  lay.regionBins = prep.regionFilter && !lay.items;
  if (lay.regionBins) {
    const first = lay.peers.find((p) => p.owner === lay.owner && p.regionBins);
    if (first && first !== lay) {
      lay.empty = 'shared';
      return;
    }
  }
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
  colorsOf(lay, prep, ctx, read);
  weightsOf(lay, prep, ctx, read);
  gather(lay, read, d.aUnit, d.bUnit);
  placeStems(lay, ctx, prep.model);
  aggregate(lay, read);
  foldStubs(lay, cfg.stubMin);
  lay.focus = (ctx.section.blocks?.length ?? 1) <= 2;
  if (!lay.nr && !lay.nb && !lay.ns) {
    // nothing read yet while the view waits for its data (EC-169: pending is said, never an empty hole), else none
    lay.empty =
      d.n === 0 && lay.reader.failed
        ? 'failed'
        : d.n === 0 && (!ctx.settled || lay.reader.waiting)
          ? 'pending'
          : 'none';
    holeWords(lay, prep, ctx);
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
  const W = rowWeights();
  const RowColor = lay.rowColor;
  lay.classOverflow = false;
  RibbonKeys.clear();
  clearClasses();
  StubKeys.clear();
  lay.rowRib = grow(lay.rowRib, d.n);
  lay.rowRib.fill(-1, 0, d.n);
  for (let k = 0; k < d.n; k++) {
    const a = RowA[k];
    if (a < 0) continue;
    const b = RowB[k];
    const c = W[k];
    const n = w ? w[k] : 1;
    const pa = lay.stemPlaced[a];
    const pb = lay.stemPlaced[b];
    if (!pa && !pb) continue;
    if (a === b) {
      lay.stemWithin[a] += c;
      lay.withinLinks[a] += n;
      continue;
    }
    if (pa && pb) {
      // a stem pair a strand: the stem earlier in the systems' order its `a` (its hue), the strand the rows' relative
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      const first = lay.stemSystem[a] < lay.stemSystem[b] || (lay.stemSystem[a] === lay.stemSystem[b] && a < b) ? a : b;
      const strand = d.aStrand[k] * d.bStrand[k] < 0 ? -1 : 1;
      // a field colour by link: one ribbon a stem pair, a strand and a colour (kinds are layers, never mixed in a ribbon)
      const key = ((lo * 1048576 + hi) * 2 + (strand < 0 ? 1 : 0)) * 64 + (lay.byRow ? classOf(lay, RowColor[k]) : 0);
      let r = RibbonKeys.get(key);
      if (r === undefined) {
        r = lay.nr++;
        if (r >= lay.ribA.length) growRibbons(lay, r + 1);
        lay.ribA[r] = first;
        lay.ribB[r] = first === a ? b : a;
        lay.ribStrand[r] = strand;
        lay.ribCount[r] = 0;
        lay.ribKind[r] = d.kind[k];
        lay.ribLinks[r] = 0;
        lay.ribColor[r] = lay.byRow ? RowColor[k] : 0;
        lay.ribRow[r] = k;
        RibbonKeys.set(key, r);
      }
      lay.ribCount[r] += c;
      lay.ribLinks[r] += n;
      lay.rowRib[k] = r;
      lay.stemCount[a] += c;
      lay.stemCount[b] += c;
      lay.stemLinks[a] += n;
      lay.stemLinks[b] += n;
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
      lay.stubLinks[s] = 0;
      StubKeys.set(key, s);
    }
    lay.stubCount[s] += c;
    lay.stubLinks[s] += n;
    lay.stemCount[near] += c;
    lay.stemLinks[near] += n;
  }
}

/** A stem's small stubs (scratch, grown): how many, its folded stub's index. */
let SmallN = new Int32Array(64);
let FoldAt = new Int32Array(64);
let TmpStem = new Int32Array(16);
let TmpPartner = new Int32Array(16);
let TmpCount = new Float64Array(16);
let TmpLinks = new Float64Array(16);

/**
 * Folds a stem's stubs of fewer than `min` links into one stub of partner −1 when it has two or more (its count their
 * sum, so a neck is still the sum of its parts), their partners kept for its card (`foldOwner`, `foldPartner`).
 */
function foldStubs(lay: RondoLay, min: number) {
  lay.nf = 0;
  const n0 = lay.nb;
  if (!(min > 1) || n0 < 2) return;
  SmallN = grow(SmallN, lay.ns);
  FoldAt = grow(FoldAt, lay.ns);
  SmallN.fill(0, 0, lay.ns);
  FoldAt.fill(-1, 0, lay.ns);
  let folds = false;
  // by links, never by a width's field (`q` sums are not counts)
  for (let b = 0; b < n0; b++) if (lay.stubLinks[b] < min && ++SmallN[lay.stubStem[b]] >= 2) folds = true;
  if (!folds) return;
  TmpStem = grow(TmpStem, n0);
  TmpPartner = grow(TmpPartner, n0);
  TmpCount = grow(TmpCount, n0);
  TmpLinks = grow(TmpLinks, n0);
  for (let b = 0; b < n0; b++) {
    TmpLinks[b] = lay.stubLinks[b];
    TmpStem[b] = lay.stubStem[b];
    TmpPartner[b] = lay.stubPartner[b];
    TmpCount[b] = lay.stubCount[b];
  }
  lay.foldOwner = grow(lay.foldOwner, n0);
  lay.foldPartner = grow(lay.foldPartner, n0);
  lay.foldCount = grow(lay.foldCount, n0);
  lay.foldLinks = grow(lay.foldLinks, n0);
  lay.nb = 0;
  for (let b = 0; b < n0; b++) {
    const stem = TmpStem[b];
    const c = TmpCount[b];
    const n = TmpLinks[b];
    if (n < min && SmallN[stem] >= 2) {
      let f = FoldAt[stem];
      if (f < 0) {
        f = FoldAt[stem] = lay.nb++;
        lay.stubStem[f] = stem;
        lay.stubPartner[f] = -1;
        lay.stubCount[f] = 0;
        lay.stubLinks[f] = 0;
      }
      lay.stubCount[f] += c;
      lay.stubLinks[f] += n;
      const e = lay.nf++;
      lay.foldOwner[e] = f;
      lay.foldPartner[e] = TmpPartner[b];
      lay.foldCount[e] = c;
      lay.foldLinks[e] = n;
      continue;
    }
    const k = lay.nb++;
    lay.stubStem[k] = stem;
    lay.stubPartner[k] = TmpPartner[b];
    lay.stubCount[k] = c;
    lay.stubLinks[k] = n;
  }
}

/** A folded stub's partners and counts, strongest first. */
export function foldedOf(lay: RondoLay, stub: number): { partner: number; count: number }[] {
  const out: { partner: number; count: number }[] = [];
  for (let e = 0; e < lay.nf; e++)
    if (lay.foldOwner[e] === stub) out.push({ partner: lay.foldPartner[e], count: lay.foldLinks[e] });
  return out.sort((a, b) => b.count - a.count || a.partner - b.partner);
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
  lay.ownMax = max;
  lay.shared = lay.prep?.scale === 'section';
  // one scale over the section (`scale: section`): the widest stem of its layers of that scale and width rule (a layer
  // laid out after this one is read as last laid out: a frame late when its widest stem grows)
  if (lay.shared)
    for (const p of lay.peers)
      if (p !== lay && p.shared && p.width === lay.width && p.read && !p.empty && p.ownMax > max) max = p.ownMax;
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
    lay.fanLinks[lay.ribA[r]] += lay.ribLinks[r];
    lay.fanLinks[lay.ribB[r]] += lay.ribLinks[r];
    lay.stemFanTo[lay.ribA[r]]++;
    lay.stemFanTo[lay.ribB[r]]++;
  }
  for (let s = 0; s < lay.nb; s++) {
    const ok = lay.k * lay.stubCount[s] >= floor;
    lay.stubDrawn[s] = ok ? 1 : 0;
    if (!ok) {
      lay.stemFan[lay.stubStem[s]] += lay.stubCount[s];
      lay.fanLinks[lay.stubStem[s]] += lay.stubLinks[s];
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
