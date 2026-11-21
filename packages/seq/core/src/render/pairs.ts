/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs layer (_spec/README.md "Render", The pairs layer): two locations that belong together, drawn by where its
// section puts them. Laid out once a section (`scope: 'section'`): each end on its block (`pairs-ends.ts`: a pair on
// one system again on each copy) placed with the placement's `place`, its arc's control the mid offset at the arc's
// height in its strip's frame through the strip's warp, so one curve is an arc over a track or a block, an arc between
// two lines, a chord through a ring's hole, with no layout named here. An end not placed leaves a caret at the placed
// end, never a dangling arc; on one straight strip an end past the window is taken `reach` (a share of the window,
// 0.25) past its edge, so an arc moves continuously on pan. An arc's height is by its span against the widest drawn
// (nested arcs nest, the widest fills the band); `side` above or below (a ring's below reaches the hole: chords).
// Crossing pairs are drawn crossing (`same`), or `split`: the largest nested set on `side`, the rest on the other side
// (exact to 2,000 pairs, greedy past it). At most `pairs.maxArcs` arcs: past it one a bin of the ends' device columns
// (widest first), the rest a density strip of ends. A hit is the arc or caret within `pairs.hitPx` of the point: the
// pair's element, its extent both ends. In a grid (a strip with `ctx.axisY`) the layer is laid out a strip at a time
// instead: a dot a pair of single positions, a diagonal for a stretch, at (a, b) and with `mirror` at (b, a) too; kinds
// at one cell side by side.

import { issueOf } from '@molstar/core/repr/tree';
import { type Key, type Location, location } from '@molstar/seqx-core/query';
import { theme as themeForm } from '@molstar/seqx-core/lang';
import type { PairsData, StaticModel } from '@molstar/seqx-core/resolve';
import type { Channel, ChannelName, ThemeDatum } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer, ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, channelLegends, Extents, type LayerChannels } from './channels.js';
import type { RenderConfig } from './config.js';
import { grown, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { compileRowFilter } from './filter.js';
import { RenderIssues } from './issues.js';
import { type PairEnds, pairEnds } from './pairs-ends.js';
import {
  type Box,
  type Emphasis,
  type EmphasisHighlight,
  type Hit,
  type ItemDescription,
  type LayerType,
  type LayoutContext,
  sectionId,
} from './layer-type.js';
import { WarpKind } from './layout/strips.js';
import { PlacePoint, type Placement } from './layout/types.js';
import type { SeqTokens } from './tokens.js';

/** Pairs a `split` places exactly (the largest nested set: quadratic); past it, greedily (fixed: no figure asks more). */
export const SplitExact = 2000;

/** An arc's least height, px: never flat. */
const MinArc = 3;

/** A dot's least side, CSS px. */
const DotMin = 4;

/** The feet's inset from the band's edge, px: the stroke's half and a hair, so an arc is drawn within its band. */
const Inset = 1;

/** The points a curve is sampled at for a hit. */
const HitSamples = 16;

/** The caret's half width at a placed end whose partner is not placed, px. */
const Caret = 3;

type Side = 'above' | 'below';

/** The default color: by kind (the `type` column), a categorical palette whose colors are fixed by value. */
const DefaultColor = themeForm.field('type', { list: false });

/** A pair as the themes and the filter read it: its kind, label, key, strands and ends. */
class PairDatum implements ThemeDatum {
  value = Number.NaN;
  letter = '';
  start = 0;
  end = 0;
  k = 0;
  key: string | number = 0;
  constructor(readonly data: PairsData) {}
  set(k: number): this {
    const d = this.data;
    this.k = k;
    this.key = d.keys[k];
    this.start = d.aStart[k];
    this.end = d.aEnd[k];
    return this;
  }
  field(name: string): string | number | undefined {
    const d = this.data;
    const k = this.k;
    switch (name) {
      case 'type':
        return d.kind[k] < 0 ? undefined : d.kinds[d.kind[k]];
      case 'label':
        return d.labels?.[k] ?? undefined;
      case 'id':
        return d.keys[k];
      case 'a_strand':
        return d.aStrand[k] < 0 ? '-' : '+';
      case 'b_strand':
        return d.bStrand[k] < 0 ? '-' : '+';
      case 'a_start':
        return d.aStart[k];
      case 'a_end':
        return d.aEnd[k];
      case 'b_start':
        return d.bStart[k];
      case 'b_end':
        return d.bEnd[k];
      default:
        return d.extra?.[name]?.[k] ?? undefined;
    }
  }
}

const PairEnds = new Set(['a_start', 'a_end', 'b_start', 'b_end']);
const PairEnvs = new WeakMap<PairsData, { numeric(n: string): boolean }>();
/** The color channel's environment: the ends and the numeric extra columns are numbers, the rest text. */
function pairEnv(d: PairsData | undefined): { numeric(n: string): boolean } {
  let env = d && PairEnvs.get(d);
  if (env) return env;
  const numbers = new Set(
    Object.entries(d?.extra ?? {}).flatMap(([k, v]) => (v.some((x) => typeof x === 'number') ? [k] : [])),
  );
  env = { numeric: (n: string) => PairEnds.has(n) || numbers.has(n) };
  if (d) PairEnvs.set(d, env);
  return env;
}

export interface PairsPrep {
  readonly owner: string;
  readonly data?: PairsData;
  readonly aSystem: string;
  readonly bSystem: string;
  /** A complex's links: the units their ends are on (`PairsData.aUnit`, `bUnit` index them), each row's own. */
  readonly endUnits?: readonly string[];
  readonly side: Side;
  readonly reach: number;
  readonly acrossLines: boolean;
  readonly mirror: boolean;
  readonly triangle: 'both' | 'upper' | 'lower';
  readonly channels: LayerChannels;
  /** 1 where the layer's data filter keeps a pair (all when none). */
  readonly shown: Uint8Array;
  /** The pairs the filter left out (the track's "n filtered"); whether its form is not read row by row (said). */
  readonly filtered: number;
  readonly unreadFilter: boolean;
  readonly datum?: PairDatum;
  /** With `split`: 1 where a pair is drawn on the other side. */
  readonly other?: Uint8Array;
  /** Each pair's index by its key. */
  readonly index: ReadonlyMap<string | number, number>;
  /** Each pair's element and extent (both ends), made when first hit. */
  readonly elements: (Key | undefined)[];
  readonly extents: (Location | undefined)[];
  /** The resolver's steps: a block of a unit shows what it instantiates (an end on a molecule reaches its units). */
  readonly model: StaticModel;
  /** The ends drawn in a section (`pairEnds`), by section. */
  readonly units: WeakMap<ResolvedSection, PairEnds>;
}

type ScratchHit = { -readonly [K in keyof Hit]: Hit[K] };

export class PairsLay {
  key = '';
  section = '';
  owner = '';
  tokens: SeqTokens | undefined;
  config: RenderConfig | undefined;
  prep: PairsPrep | undefined;
  dpr = 1;
  /** `arcs` (a section's) or `dots` (a grid's strip). */
  mode: 'arcs' | 'dots' | 'none' = 'none';
  /** The color channel as last laid out, its open domains' extents. */
  colorChannel: Channel | undefined;
  readonly extents = new Extents();
  /** Arcs: the pair's index, and its curve on the page (start, two controls, end: 8 numbers), its color. */
  n = 0;
  pair = new Int32Array(64);
  curve = new Float64Array(512);
  color = new Uint32Array(64);
  /** Carets at ends whose partner is not placed: x, y, 1 above or −1 below, the pair. */
  carets: number[] = [];
  /** The density strip of ends past `maxArcs`: counts a device column from `densityX`, at `densityY`. */
  density = new Uint16Array(0);
  densityN = 0;
  densityX = 0;
  densityY = 0;
  densityH = 0;
  densityMax = 0;
  /** The page's box the arcs are clipped to (the data's, for translated strips). */
  clip: Box | undefined;
  /** Dots: x, y, w, h each, the pair, its color; segments: x0, y0, x1, y1, width each, the pair, its color. */
  dots = new Float64Array(64);
  dotPair = new Int32Array(16);
  dotColor = new Uint32Array(16);
  nd = 0;
  segs = new Float64Array(80);
  segPair = new Int32Array(16);
  segColor = new Uint32Array(16);
  ns = 0;
  readonly scratch: ScratchHit = { layer: '', section: '', system: '' };
  /** A path's points for a segment (two). */
  readonly seg = new Float64Array(4);
}

/** A lay's pairs lit by the highlight channels (module-private): what they were made for, and each pair's channel. */
interface Lit {
  /** The highlight channels and the prep `lit` was made for (the viewer keeps the list until a channel changes). */
  for: readonly EmphasisHighlight[] | undefined;
  prep: PairsPrep | undefined;
  /** Each pair's highlight channel (its index + 1; 0 none): the last channel listing its key. */
  lit: Int32Array;
}

const LitOf = new WeakMap<PairsLay, Lit>();

/** The annotation a layer's data names: its address, else (a query) the annotation it starts at. */
function annotationOf(model: StaticModel, address: string) {
  const a = model.annotations.get(address);
  if (a) return { address, entry: a };
  for (const [at, entry] of model.annotations) if (address.startsWith(`${at}/`)) return { address: at, entry };
  return undefined;
}

function prepare(layer: ResolvedLayer, model: StaticModel): PairsPrep {
  const named = annotationOf(model, layer.data?.address ?? '');
  const owner = named?.address ?? '';
  const a = named?.entry;
  const data = a?.shape === 'pairs' ? a.pairs : undefined;
  const p = layer.params;
  const index = new Map<string | number, number>();
  if (data) for (let k = 0; k < data.n; k++) if (!index.has(data.keys[k])) index.set(data.keys[k], k);
  const t = p.triangle;
  // the layer's data filter, read pair by pair once (_spec/README.md "Themes", Filters)
  const form = layer.data?.filter;
  const filter = form ? compileRowFilter(form) : undefined;
  const baseForm = layer.data?.base;
  const base = baseForm ? compileRowFilter(baseForm) : undefined;
  const datum = data ? new PairDatum(data) : undefined;
  const shown = new Uint8Array(data?.n ?? 0).fill(1);
  let filtered = 0;
  if ((filter || base) && datum)
    for (let k = 0; k < shown.length; k++)
      if (base && !base(datum.set(k))) shown[k] = 0;
      else if (filter && !filter(datum.set(k))) {
        shown[k] = 0;
        filtered++;
      }
  return {
    owner,
    ...(data ? { data } : {}),
    ...(datum ? { datum } : {}),
    channels: channelForms(layer.theme),
    shown,
    filtered,
    unreadFilter: (form !== undefined && !filter) || (baseForm !== undefined && !base) || layer.data?.unread === true,
    aSystem: a?.target ?? '',
    bSystem: a?.ySystem ?? '',
    ...(a?.units ? { endUnits: a.units } : {}),
    side: p.side === 'below' ? 'below' : 'above',
    reach: typeof p.reach === 'number' && p.reach >= 0 ? p.reach : 0.25,
    acrossLines: p.across_lines === true,
    mirror: p.mirror === true,
    triangle: t === 'upper' || t === 'lower' ? t : 'both',
    ...(data && p.crossing === 'split' && a?.target === a?.ySystem ? { other: splitCrossing(data, shown) } : {}),
    index,
    elements: [],
    extents: [],
    model,
    units: new WeakMap(),
  };
}

/**
 * The pairs drawn on the other side under `split`: the complement of the largest nested (non-crossing) set of the
 * pairs as intervals between their ends' centres; exact by the interval recurrence to `SplitExact` pairs, greedy past
 * it (narrowest first, kept when it crosses none kept).
 */
export function splitCrossing(d: PairsData, shown?: Uint8Array): Uint8Array {
  // the pairs drawn (a filter's): those left out cross nothing
  const at: number[] = [];
  for (let k = 0; k < d.n; k++) if (!shown || shown[k]) at.push(k);
  const n = at.length;
  const lo = new Float64Array(n);
  const hi = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const k = at[i];
    const ca = (d.aStart[k] + d.aEnd[k]) / 2;
    const cb = (d.bStart[k] + d.bEnd[k]) / 2;
    lo[i] = Math.min(ca, cb);
    hi[i] = Math.max(ca, cb);
  }
  const keep = new Uint8Array(n);
  if (n <= SplitExact) nestedExact(n, lo, hi, keep);
  else nestedGreedy(n, lo, hi, keep);
  const other = new Uint8Array(d.n);
  for (let i = 0; i < n; i++) other[at[i]] = keep[i] ? 0 : 1;
  return other;
}

/** Whether two intervals cross: one's ends either side of one end of the other. */
export function crosses(alo: number, ahi: number, blo: number, bhi: number): boolean {
  return (alo < blo && blo < ahi && ahi < bhi) || (blo < alo && alo < bhi && bhi < ahi);
}

/**
 * The largest set of non-crossing intervals (nested or apart): over the 2n ends in order, F(i, j) the most held within
 * ends i..j, F(i, j) = max(F(i + 1, j), 1 + F(i + 1, k − 1) + F(k + 1, j)) where k is i's partner in (i, j]. At one
 * place closing ends come first, then each self-pair's two ends, then opening ends; a wider interval outside a narrower
 * one and equal ones nested, as `crosses` reads them (touching ends do not cross). O(n²) time, a triangle of int16.
 */
function nestedExact(n: number, lo: Float64Array, hi: Float64Array, keep: Uint8Array) {
  const m = 2 * n;
  const ends = new Int32Array(m);
  for (let i = 0; i < m; i++) ends[i] = i;
  const at = (e: number) => (e < n ? lo[e] : hi[e - n]);
  // at one place: the others' closing ends, each self-pair's two ends together (it crosses nothing), the others'
  // opening ends; equal intervals nest (opened in order, closed in reverse)
  const rank = (e: number) => (lo[e % n] === hi[e % n] ? 1 : e < n ? 2 : 0);
  ends.sort((x, y) => {
    const d = at(x) - at(y);
    if (d) return d;
    const r = rank(x) - rank(y);
    if (r) return r;
    const px = x % n;
    const py = y % n;
    if (r === 0 && rank(x) === 1) return px !== py ? px - py : x < n ? -1 : 1;
    const w = hi[px] - lo[px] - (hi[py] - lo[py]);
    if (w) return x < n ? -w : w;
    return x < n ? px - py : py - px;
  });
  const pos = new Int32Array(m);
  for (let i = 0; i < m; i++) pos[ends[i]] = i;
  const partner = new Int32Array(m);
  for (let i = 0; i < m; i++) {
    const e = ends[i];
    partner[i] = pos[e < n ? e + n : e - n];
  }
  // F over the triangle i ≤ j, row i from `row(i)`
  const row = (i: number) => i * m - (i * (i - 1)) / 2 - i;
  const F = new Int16Array((m * (m + 1)) / 2);
  const get = (i: number, j: number) => (i > j ? 0 : F[row(i) + j]);
  for (let len = 2; len <= m; len++)
    for (let i = 0; i + len - 1 < m; i++) {
      const j = i + len - 1;
      let best = get(i + 1, j);
      const k = partner[i];
      if (k > i && k <= j) {
        const v = 1 + get(i + 1, k - 1) + get(k + 1, j);
        if (v > best) best = v;
      }
      F[row(i) + j] = best;
    }
  // back through the choices
  const stack: number[] = [0, m - 1];
  while (stack.length) {
    const j = stack.pop()!;
    let i = stack.pop()!;
    while (i < j) {
      const k = partner[i];
      if (k > i && k <= j && get(i, j) === 1 + get(i + 1, k - 1) + get(k + 1, j)) {
        keep[ends[i] % n] = 1;
        stack.push(i + 1, k - 1);
        i = k + 1;
      } else i++;
    }
  }
}

/** Narrowest first, a pair kept when it crosses no pair kept. */
function nestedGreedy(n: number, lo: Float64Array, hi: Float64Array, keep: Uint8Array) {
  const order = Array.from({ length: n }, (_, k) => k).sort((x, y) => hi[x] - lo[x] - (hi[y] - lo[y]));
  const kept: number[] = [];
  for (const k of order) {
    let ok = true;
    for (let i = 0; i < kept.length && ok; i++) if (crosses(lo[k], hi[k], lo[kept[i]], hi[kept[i]])) ok = false;
    if (ok) {
      kept.push(k);
      keep[k] = 1;
    }
  }
}

/** Each pair's ends in a section, kept by the section and the pieces they were placed on. */
function endsIn(prep: PairsPrep, section: ResolvedSection, layout: AxisTransform): PairEnds {
  const kept = prep.units.get(section);
  if (kept && kept.pieces === (section.axis.ranges ?? layout.pieces)) return kept;
  const ends = pairEnds({ ...prep, data: prep.data! }, section, layout);
  prep.units.set(section, ends);
  return ends;
}

/** A pair's color from the layer's channel (its kind's by default); the ranges' ink where it gives none. */
function pairColor(lay: PairsLay, prep: PairsPrep, k: number): number {
  const c = lay.colorChannel?.at(prep.datum!.set(k));
  return typeof c === 'number' ? c : lay.tokens!.colors['seq-range'];
}

const PA = new PlacePoint();
const PB = new PlacePoint();
const PM = new PlacePoint();

export const PairsType: LayerType<PairsPrep, PairsLay> = {
  name: 'pairs',
  fallback: 'ranges',
  scope: 'section',
  prepare,

  layout(prep, ctx: LayoutContext, _data, _overlay, prev?: PairsLay): PairsLay {
    const lay = prev ?? new PairsLay();
    lay.key = ctx.key;
    lay.section = sectionId(ctx.section);
    lay.owner = prep.owner;
    lay.tokens = ctx.tokens;
    lay.config = ctx.config;
    lay.prep = prep;
    lay.dpr = ctx.dpr;
    lay.n = 0;
    lay.nd = 0;
    lay.ns = 0;
    lay.carets.length = 0;
    lay.densityN = 0;
    lay.clip = undefined;
    lay.mode = 'none';
    lay.colorChannel = undefined;
    if (!prep.data) return lay;
    if (prep.unreadFilter)
      ctx.report(
        issueOf(
          RenderIssues,
          'seq.render.filter',
          'The layer’s data query is not read pair by pair: no pair is left out',
          {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.filter'].hint,
          },
        ),
      );
    lay.colorChannel = prep.channels.get('color' as ChannelName, DefaultColor, ctx, pairEnv(prep.data));
    lay.extents.reset([lay.colorChannel]);
    if (ctx.axisY) layoutDots(lay, prep, ctx, ctx.axis, ctx.axisY);
    else if (ctx.place) layoutArcs(lay, prep, ctx, ctx.place);
    if (lay.extents.any) lay.extents.flush(ctx.theme);
    return lay;
  },

  filtered: (lay) => lay.prep?.filtered ?? 0,

  legend(lay, options) {
    return channelLegends([lay.colorChannel], { ...options, layer: lay.key });
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    if (!lay.tokens || !lay.config || !lay.prep) return;
    if (pass === 'data') {
      if (lay.mode === 'arcs') drawArcs(lay, out);
      else if (lay.mode === 'dots') drawDots(lay, out);
      return;
    }
    if (emphasis && (pass === 'hover' || pass === 'selection')) drawEmphasis(lay, pass, out, emphasis);
    else if (pass === 'highlight' && emphasis?.highlights) drawHighlights(lay, out, emphasis.highlights);
  },

  items: pairItems,

  hit(lay, at) {
    if (lay.mode === 'arcs') return hitArcs(lay, at);
    if (lay.mode === 'dots') return hitDots(lay, at);
    return undefined;
  },
};

// ---- arcs ----

/** The page point of a strip-frame point: through its translation, or its polar warp. */
function toPage(place: Placement, strip: number, x: number, y: number, out: Float64Array, o: number) {
  const s = place.strips;
  if (s.warp[strip] === WarpKind.polar) {
    const w = place.warpOf!(strip);
    w.map(x, y);
    out[o] = w.px;
    out[o + 1] = w.py;
  } else {
    out[o] = x + s.dx[strip];
    out[o + 1] = y + s.dy[strip];
  }
}

/** A display unit's strip (into `pt`) and its frame x, NaN when no strip holds it. */
function frameX(place: Placement, u: number, pt: PlacePoint): number {
  if (!place.place(u, pt)) return Number.NaN;
  const s = place.strips;
  const i = pt.strip;
  // a polar strip's point is on the page: its frame x from the strip's first unit (a ring's arcs)
  if (s.warp[i] === WarpKind.polar) return s.axes[i].xOfUnit((u - s.first[i]) / s.scale[i]);
  return pt.x - s.dx[i];
}

/** What a layout of arcs reads a pair at a time (kept: grown, never shrunk). */
let Spans = new Float64Array(64);
let StripOf = new Int32Array(128);
let Raw = new Float64Array(256);

function layoutArcs(lay: PairsLay, prep: PairsPrep, ctx: LayoutContext, place: Placement) {
  const s = place.strips;
  if (!s.n) return;
  lay.mode = 'arcs';
  const band = ctx.band ?? 0;
  const ends = endsIn(prep, ctx.section, ctx.axis);
  const units = ends.u;
  const nb = s.bands;
  const one = s.n === 1 && s.warp[0] === WarpKind.translate;
  const a0 = s.axes[0];
  const w = a0.u1 - a0.u0;
  const reachLo = a0.u0 - prep.reach * w;
  const reachHi = a0.u1 + prep.reach * w;
  if (Spans.length < ends.n) {
    Spans = new Float64Array(2 * ends.n);
    StripOf = new Int32Array(4 * ends.n);
    Raw = new Float64Array(8 * ends.n);
  }
  let m = 0;
  let widest = 0;
  let unplaced = false;
  let copied = 0;
  // each entry a pair, or a pair again on a copy (`pairEnds`)
  for (let e = 0; e < ends.n; e++) {
    const k = ends.pair[e];
    if (!prep.shown[k]) continue;
    let oa = units[2 * e];
    let ob = units[2 * e + 1];
    const unnamed = ends.unnamed[e];
    if (unnamed) {
      // an end on split copies the data names no unit of: on no copy (`pairEnds`)
      copied++;
      if (unnamed & 4) unplaced = true;
      continue;
    }
    if (Number.isNaN(oa) && Number.isNaN(ob)) {
      unplaced = true;
      continue;
    }
    let xa = !Number.isNaN(oa) ? frameX(place, oa, PA) : Number.NaN;
    let xb = !Number.isNaN(ob) ? frameX(place, ob, PB) : Number.NaN;
    let sa = PA.strip;
    let sb = PB.strip;
    if (one) {
      // one straight strip: an end past the window taken `reach` past its edge
      const outA = !Number.isNaN(oa) && Number.isNaN(xa);
      const outB = !Number.isNaN(ob) && Number.isNaN(xb);
      if (outA && outB && oa < a0.u0 === ob < a0.u0) continue;
      if (outA) {
        oa = Math.max(reachLo, Math.min(reachHi, oa));
        xa = a0.xOfUnit(oa);
        sa = 0;
      }
      if (outB) {
        ob = Math.max(reachLo, Math.min(reachHi, ob));
        xb = a0.xOfUnit(ob);
        sb = 0;
      }
    }
    if (Number.isNaN(xa) && Number.isNaN(xb)) continue;
    // an end not placed (a chain not shown, another line without `across_lines`): a caret at the placed one
    const across = !Number.isNaN(xa) && !Number.isNaN(xb) && s.top[sa] !== s.top[sb] && !prep.acrossLines;
    if (Number.isNaN(xa) || Number.isNaN(xb) || across) {
      if (!Number.isNaN(xa)) caret(lay, prep, place, sa, band, xa, k);
      if (!Number.isNaN(xb)) caret(lay, prep, place, sb, band, xb, k);
      continue;
    }
    // the span: px on a straight strip, display units round a ring
    const span = s.warp[sa] === WarpKind.polar ? Math.abs(ob - oa) : Math.abs(xb - xa);
    Spans[m] = span;
    StripOf[2 * m] = sa;
    StripOf[2 * m + 1] = sb;
    Raw[4 * m] = xa;
    Raw[4 * m + 1] = xb;
    Raw[4 * m + 2] = (oa + ob) / 2;
    Raw[4 * m + 3] = k;
    if (span > widest) widest = span;
    m++;
  }
  if (copied) reportCopies(ctx, prep, ends, copied);
  if (unplaced)
    ctx.report(
      issueOf(RenderIssues, 'seq.render.pairs.systems', 'Some pairs’ ends are on no system the section shows', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.pairs.systems'].hint,
      }),
    );
  // past `maxArcs`: one arc a bin of the ends' device columns, the widest first; the rest a density strip of ends
  const max = ctx.config.pairs.maxArcs;
  const binned = m > max;
  let order: Int32Array | undefined;
  if (binned) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.pairs.limit', `More than ${max} pairs in view: binned`, {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.pairs.limit'].hint,
      }),
    );
    order = new Int32Array(m);
    for (let i = 0; i < m; i++) order[i] = i;
    order.sort((x, y) => Spans[y] - Spans[x]);
    densityStrip(lay, place, band, ctx.dpr);
  }
  const seen = binned ? new Set<number>() : undefined;
  const dpr = ctx.dpr;
  for (let oi = 0; oi < m; oi++) {
    const i = order ? order[oi] : oi;
    const k = Raw[4 * i + 3];
    const sa = StripOf[2 * i];
    const sb = StripOf[2 * i + 1];
    let xa = Raw[4 * i];
    let xb = Raw[4 * i + 1];
    if (seen) {
      const ca = Math.floor((xa + s.dx[sa]) * dpr);
      const cb = Math.floor((xb + s.dx[sb]) * dpr);
      const bin = Math.min(ca, cb) * 65536 + Math.max(ca, cb);
      if (lay.n >= max || seen.has(bin)) {
        addDensity(lay, xa + s.dx[sa], dpr);
        addDensity(lay, xb + s.dx[sb], dpr);
        continue;
      }
      seen.add(bin);
    }
    const other = prep.other?.[k] === 1;
    const side: Side = other ? (prep.side === 'above' ? 'below' : 'above') : prep.side;
    const kb = sa * nb + band;
    const top = s.bandY[kb];
    const h = s.bandH[kb];
    const polar = s.warp[sa] === WarpKind.polar;
    // the feet on the band's edge toward the sequence (`above`: its bottom, the arcs rising; `below`: its top), inset
    // so the stroke stays in the band; split, on its middle, each side half the band (a ring: its edges, as drawn)
    const split = prep.other !== undefined && !polar;
    const base = split
      ? top + h / 2
      : polar
        ? side === 'above'
          ? top + h
          : top
        : side === 'above'
          ? top + h - Inset
          : top + Inset;
    const room = split ? h / 2 - Inset : h - 2 * Inset;
    const frac = widest > 0 ? Spans[i] / widest : 1;
    // the control's offset from the base: twice the apex's height (a curve's apex is halfway to its control); on a
    // ring below, toward the centre by the span (a chord through the hole)
    const depth =
      polar && side === 'below'
        ? (place.warpOf!(sa).r0 - base) * Math.min(1, frac) * 1.5
        : 2 * Math.min(Math.max(0, room), Math.max(MinArc, room * frac));
    const cy = side === 'above' ? base - depth : base + depth;
    lay.pair = grown(lay.pair, lay.n + 1);
    lay.color = grown(lay.color, lay.n + 1);
    lay.curve = grown(lay.curve, 8 * lay.n + 8);
    const o = 8 * lay.n;
    const c = lay.curve;
    if (xa > xb && !polar && sa === sb) {
      const t = xa;
      xa = xb;
      xb = t;
    }
    toPage(place, sa, xa, base, c, o);
    toPage(place, sb, xb, base, c, o + 6);
    // the control: the mid offset at the arc's height in its strip's frame (one strip: between the ends)
    if (sa === sb) {
      const mid = polar ? s.axes[sa].xOfUnit((Raw[4 * i + 2] - s.first[sa]) / s.scale[sa]) : (xa + xb) / 2;
      toPage(place, sa, mid, cy, c, o + 2);
    } else {
      const mx = frameX(place, Raw[4 * i + 2], PM);
      if (!Number.isNaN(mx)) toPage(place, PM.strip, mx, cy, c, o + 2);
      else {
        c[o + 2] = (c[o] + c[o + 6]) / 2;
        c[o + 3] = Math.min(c[o + 1], c[o + 7]) - depth;
      }
    }
    // a quadratic's control as a cubic's two
    const qx = c[o + 2];
    const qy = c[o + 3];
    c[o + 2] = c[o] + (2 / 3) * (qx - c[o]);
    c[o + 3] = c[o + 1] + (2 / 3) * (qy - c[o + 1]);
    c[o + 4] = c[o + 6] + (2 / 3) * (qx - c[o + 6]);
    c[o + 5] = c[o + 7] + (2 / 3) * (qy - c[o + 7]);
    lay.pair[lay.n] = k;
    lay.color[lay.n] = pairColor(lay, prep, k);
    if (lay.extents.any) lay.extents.add(prep.datum!.set(k));
    lay.n++;
  }
  // translated strips clip to the data's box (an end taken past the window); a ring needs none
  if (s.warp[0] === WarpKind.translate) {
    let x0 = Number.POSITIVE_INFINITY;
    let x1 = Number.NEGATIVE_INFINITY;
    for (let i = 0; i < s.n; i++) {
      x0 = Math.min(x0, s.axes[i].x0 + s.dx[i]);
      x1 = Math.max(x1, s.axes[i].x0 + s.axes[i].width + s.dx[i]);
    }
    const b = place.box;
    lay.clip = { x: x0, y: b.y, width: Math.max(0, x1 - x0), height: b.height };
  }
}

/** The pairs on split copies whose data names no units: the count, each system with its copies. */
function reportCopies(ctx: LayoutContext, prep: PairsPrep, ends: PairEnds, n: number) {
  const name = (system: string, copies: number) => `${system.slice(system.lastIndexOf('/') + 1)} ×${copies}`;
  const shown = [
    ...(ends.aCopies ? [name(prep.aSystem, ends.aCopies)] : []),
    ...(ends.bCopies && prep.bSystem !== prep.aSystem ? [name(prep.bSystem, ends.bCopies)] : []),
  ];
  const pairs = n === 1 ? 'pair' : 'pairs';
  const message = `${n} ${pairs} on split copies (${shown.join(', ')}) not drawn: the data names no units`;
  ctx.report(
    issueOf(RenderIssues, 'seq.render.pairs.copy', message, {
      address: ctx.key,
      hint: RenderIssues.codes['seq.render.pairs.copy'].hint,
    }),
  );
}

/** A caret at a placed end whose partner is not placed (a tick into the band), its pair's. */
function caret(lay: PairsLay, prep: PairsPrep, place: Placement, strip: number, band: number, x: number, k: number) {
  const s = place.strips;
  const kb = strip * s.bands + band;
  const above = prep.side === 'above';
  const y = above ? s.bandY[kb] + s.bandH[kb] : s.bandY[kb];
  toPage(place, strip, x, y, CaretAt, 0);
  lay.carets.push(CaretAt[0], CaretAt[1], above ? 1 : -1, k);
}

const CaretAt = new Float64Array(2);

function densityStrip(lay: PairsLay, place: Placement, band: number, dpr: number) {
  const s = place.strips;
  const a = s.axes[0];
  const cols = Math.ceil(a.width * dpr) + 1;
  if (lay.density.length < cols) lay.density = new Uint16Array(cols);
  lay.density.fill(0, 0, cols);
  lay.densityN = cols;
  lay.densityX = a.x0 + s.dx[0];
  lay.densityH = Math.min(3, s.bandH[band]);
  lay.densityY = s.bandY[band] + s.bandH[band] - lay.densityH + s.dy[0];
  lay.densityMax = 0;
}

function addDensity(lay: PairsLay, x: number, dpr: number) {
  const c = Math.floor((x - lay.densityX) * dpr);
  if (c < 0 || c >= lay.densityN) return;
  const v = lay.density[c] < 65535 ? ++lay.density[c] : 65535;
  if (v > lay.densityMax) lay.densityMax = v;
}

function clipTo(out: PassWriter, b: Box) {
  const g = out.args;
  g[0] = b.x;
  g[1] = b.y;
  g[2] = b.width;
  g[3] = b.height;
  out.clipArgs();
}

function drawArcs(lay: PairsLay, out: PassWriter) {
  const g = out.args;
  if (lay.clip) clipTo(out, lay.clip);
  for (let i = 0; i < lay.n; i++) {
    g[4] = NoColor;
    g[5] = lay.color[i];
    g[6] = 1;
    g[7] = 0;
    out.curveArgs(lay.curve, 8 * i);
  }
  const ink = lay.tokens!.colors['seq-range'];
  for (let i = 0; i + 3 < lay.carets.length; i += 4) {
    const x = lay.carets[i];
    const y = lay.carets[i + 1];
    g[0] = x - Caret / 2;
    g[1] = lay.carets[i + 2] > 0 ? y - 2 * Caret : y;
    g[2] = Caret;
    g[3] = 2 * Caret;
    g[4] = ink;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  if (lay.densityN && lay.densityMax) {
    const dpr = lay.dpr;
    const dens = lay.tokens!.colors['seq-density-high'];
    for (let c = 0; c < lay.densityN; c++) {
      const v = lay.density[c];
      if (!v) continue;
      g[0] = lay.densityX + c / dpr;
      g[1] = lay.densityY;
      g[2] = 1 / dpr;
      g[3] = lay.densityH;
      g[4] = ((dens & 0xffffff00) | Math.round(64 + (191 * v) / lay.densityMax)) >>> 0;
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
  if (lay.clip) out.unclip();
}

/**
 * The squared distance from a point to a curve, sampled, into `Dist[0]`: no number returned (a double returned from a
 * call that is not inlined is boxed: a hover over many arcs would allocate a number an arc).
 */
const Dist = new Float64Array(1);

function curveDistance(c: Float64Array, o: number, at: Float64Array): void {
  const x = at[0];
  const y = at[1];
  let best = Number.POSITIVE_INFINITY;
  let px = c[o];
  let py = c[o + 1];
  for (let s = 1; s <= HitSamples; s++) {
    const t = s / HitSamples;
    const u = 1 - t;
    const qx = u * u * u * c[o] + 3 * u * u * t * c[o + 2] + 3 * u * t * t * c[o + 4] + t * t * t * c[o + 6];
    const qy = u * u * u * c[o + 1] + 3 * u * u * t * c[o + 3] + 3 * u * t * t * c[o + 5] + t * t * t * c[o + 7];
    const d2 = segmentDistance2(px, py, qx, qy, x, y);
    if (d2 < best) best = d2;
    px = qx;
    py = qy;
  }
  Dist[0] = best;
}

/** The squared distance from (x, y) to the segment (x0, y0)–(x1, y1) (small: inlined in the hit loops). */
function segmentDistance2(x0: number, y0: number, x1: number, y1: number, x: number, y: number): number {
  const dx = x1 - x0;
  const dy = y1 - y0;
  const l2 = dx * dx + dy * dy;
  let f = l2 > 0 ? ((x - x0) * dx + (y - y0) * dy) / l2 : 0;
  f = f < 0 ? 0 : f > 1 ? 1 : f;
  const ex = x0 + f * dx - x;
  const ey = y0 + f * dy - y;
  return ex * ex + ey * ey;
}

/** The arc or caret nearest a page point within `pairs.hitPx` (a tie: the later, drawn on top). */
function hitArcs(lay: PairsLay, at: Float64Array): Hit | undefined {
  const x = at[0];
  const y = at[1];
  const r = lay.config!.pairs.hitPx;
  const c = lay.curve;
  if (lay.clip) {
    const b = lay.clip;
    if (x < b.x || x >= b.x + b.width || y < b.y || y >= b.y + b.height) return undefined;
  }
  let best = -1;
  let bestD = r * r;
  for (let i = 0; i < lay.n; i++) {
    const o = 8 * i;
    // the hull of its four points (the curve lies in it): x, then y
    if (outside(c, o, x, r) || outside(c, o + 1, y, r)) continue;
    curveDistance(c, o, at);
    const d = Dist[0];
    if (d <= bestD) {
      bestD = d;
      best = i;
    }
  }
  // a caret (drawn after the arcs: a tie is its), by the distance to its box
  let pair = best < 0 ? -1 : lay.pair[best];
  const cs = lay.carets;
  for (let i = 0; i + 3 < cs.length; i += 4) {
    const top = cs[i + 2] > 0 ? cs[i + 1] - 2 * Caret : cs[i + 1];
    const dx = Math.max(0, Math.abs(x - cs[i]) - Caret / 2);
    const dy = y < top ? top - y : y > top + 2 * Caret ? y - top - 2 * Caret : 0;
    if (dx * dx + dy * dy <= bestD) {
      bestD = dx * dx + dy * dy;
      pair = cs[i + 3];
    }
  }
  return pair < 0 ? undefined : pairHit(lay, pair);
}

/** Whether `v` is more than `r` outside the range of a curve's four points' coordinates at `o`, `o + 2`, … */
function outside(c: Float64Array, o: number, v: number, r: number): boolean {
  let lo = c[o];
  let hi = lo;
  for (let k = 2; k < 8; k += 2) {
    const w = c[o + k];
    if (w < lo) lo = w;
    if (w > hi) hi = w;
  }
  return v < lo - r || v > hi + r;
}

/** A pair's `a` end's system: its row's unit under a complex, else the layer's. */
function aOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.aUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.aSystem) : prep.aSystem;
}

/** A pair's `b` end's system: its row's unit under a complex, else the layer's. */
function bOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.bUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.bSystem) : prep.bSystem;
}

/** A pair's hit: its element and its extent, both ends (made once a pair). */
function pairHit(lay: PairsLay, k: number): Hit {
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
  h.section = lay.section;
  h.system = as;
  h.element = element;
  h.extent = extent;
  h.position = undefined;
  h.ySystem = undefined;
  h.yPosition = undefined;
  h.target = undefined;
  return h;
}

// ---- dots (a grid's strip) ----

function layoutDots(lay: PairsLay, prep: PairsPrep, ctx: LayoutContext, a: AxisTransform, ay: AxisTransform) {
  const d = prep.data!;
  lay.mode = 'dots';
  const forward = prep.aSystem === a.system && prep.bSystem === ay.system;
  const back = prep.aSystem === ay.system && prep.bSystem === a.system;
  if (!forward && !(prep.mirror && back)) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.pairs.systems', 'The pairs are not on the grid’s two systems', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.pairs.systems'].hint,
      }),
    );
    return;
  }
  for (let k = 0; k < d.n; k++) {
    if (!prep.shown[k]) continue;
    // a stretch runs with its partner when both ends are on the same strand, against it otherwise (EC-206)
    const parallel = d.aStrand[k] === d.bStrand[k];
    if (forward) dot(lay, prep, a, ay, k, d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k], parallel);
    // a self-pair on a self grid is its own mirror: one dot
    const self = forward && d.aStart[k] === d.bStart[k] && d.aEnd[k] === d.bEnd[k];
    if (prep.mirror && back && !self)
      dot(lay, prep, a, ay, k, d.bStart[k], d.bEnd[k], d.aStart[k], d.aEnd[k], parallel);
  }
  sideBySide(lay);
}

/**
 * A pair's dot (single positions) or diagonal (a stretch: x's first with y's first when `parallel`, else with y's last),
 * kept by the triangle.
 */
function dot(
  lay: PairsLay,
  prep: PairsPrep,
  a: AxisTransform,
  ay: AxisTransform,
  k: number,
  xs: number,
  xe: number,
  ys: number,
  ye: number,
  parallel: boolean,
) {
  const t = prep.triangle;
  const cx = (xs + xe - 1) / 2;
  const cy = (ys + ye - 1) / 2;
  if ((t === 'upper' && !(cy < cx)) || (t === 'lower' && !(cy > cx))) return;
  const color = pairColor(lay, prep, k);
  if (lay.extents.any) lay.extents.add(prep.datum!.set(k));
  if (xe - xs === 1 && ye - ys === 1) {
    if (!a.cell(xs)) return;
    const ax0 = a.sx0;
    const ax1 = a.sx1;
    if (!ay.cell(ys)) return;
    const i = lay.nd++;
    lay.dots = grown(lay.dots, 4 * i + 4);
    lay.dotPair = grown(lay.dotPair, i + 1);
    lay.dotColor = grown(lay.dotColor, i + 1);
    // at least `DotMin` CSS px a side: a dot reads at any zoom
    const m = DotMin;
    const w = Math.max(Math.abs(ax1 - ax0), m);
    const h = Math.max(Math.abs(ay.sx1 - ay.sx0), m);
    lay.dots[4 * i] = (ax0 + ax1) / 2 - w / 2;
    lay.dots[4 * i + 1] = (ay.sx0 + ay.sx1) / 2 - h / 2;
    lay.dots[4 * i + 2] = w;
    lay.dots[4 * i + 3] = h;
    lay.dotPair[i] = k;
    lay.dotColor[i] = color;
    return;
  }
  // a stretch: from x's first to x's last, y's first to last (parallel) or last to first, through the cells' centres
  if (!a.cell(xs)) return;
  const x0 = (a.sx0 + a.sx1) / 2;
  if (!a.cell(xe - 1)) return;
  const x1 = (a.sx0 + a.sx1) / 2;
  if (!ay.cell(parallel ? ys : ye - 1)) return;
  const y0 = (ay.sx0 + ay.sx1) / 2;
  if (!ay.cell(parallel ? ye - 1 : ys)) return;
  const y1 = (ay.sx0 + ay.sx1) / 2;
  const i = lay.ns++;
  lay.segs = grown(lay.segs, 5 * i + 5);
  lay.segPair = grown(lay.segPair, i + 1);
  lay.segColor = grown(lay.segColor, i + 1);
  lay.segs[5 * i] = x0;
  lay.segs[5 * i + 1] = y0;
  lay.segs[5 * i + 2] = x1;
  lay.segs[5 * i + 3] = y1;
  lay.segs[5 * i + 4] = Math.max(1.5, Math.min(a.scale, ay.scale) * 0.6);
  lay.segPair[i] = k;
  lay.segColor[i] = color;
}

/** Dots at one cell share it: each a slice of its width, in order (EC-63). */
function sideBySide(lay: PairsLay) {
  if (lay.nd < 2) return;
  const d = lay.dots;
  const at = new Map<string, number[]>();
  for (let i = 0; i < lay.nd; i++) {
    const key = `${d[4 * i]}|${d[4 * i + 1]}`;
    let l = at.get(key);
    if (!l) at.set(key, (l = []));
    l.push(i);
  }
  for (const l of at.values()) {
    if (l.length < 2) continue;
    const w = d[4 * l[0] + 2] / l.length;
    for (let j = 0; j < l.length; j++) {
      d[4 * l[j]] += j * w;
      d[4 * l[j] + 2] = w;
    }
  }
}

function drawDots(lay: PairsLay, out: PassWriter) {
  const g = out.args;
  // outlined in the text's ink (dark on light, light on dark): a dot reads over any color of the cells under it
  const line = lay.tokens!.colors['msui-text'];
  for (let i = 0; i < lay.nd; i++) {
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = lay.dotColor[i];
    g[5] = line;
    g[6] = 1;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  const p = lay.seg;
  for (let i = 0; i < lay.ns; i++) {
    p[0] = lay.segs[5 * i];
    p[1] = lay.segs[5 * i + 1];
    p[2] = lay.segs[5 * i + 2];
    p[3] = lay.segs[5 * i + 3];
    g[4] = NoColor;
    g[5] = lay.segColor[i];
    g[6] = lay.segs[5 * i + 4];
    g[7] = 0;
    out.pathArgs(p, 2, false);
  }
}

/** The dot under a point (its box), else the diagonal within `pairs.hitPx` of it. */
function hitDots(lay: PairsLay, at: Float64Array): Hit | undefined {
  const x = at[0];
  const y = at[1];
  const d = lay.dots;
  for (let i = lay.nd - 1; i >= 0; i--)
    if (x >= d[4 * i] && x < d[4 * i] + d[4 * i + 2] && y >= d[4 * i + 1] && y < d[4 * i + 1] + d[4 * i + 3])
      return pairHit(lay, lay.dotPair[i]);
  const r = lay.config!.pairs.hitPx;
  const s = lay.segs;
  for (let i = lay.ns - 1; i >= 0; i--) {
    const w = Math.max(r, s[5 * i + 4] / 2);
    if (segmentDistance2(s[5 * i], s[5 * i + 1], s[5 * i + 2], s[5 * i + 3], x, y) <= w * w)
      return pairHit(lay, lay.segPair[i]);
  }
  return undefined;
}

// ---- emphasis ----

/** Whether a pair has an end holding position `p` on `system` (a hovered base: its partner lights). */
function endHolds(prep: PairsPrep, k: number, system: string, p: number): boolean {
  const d = prep.data!;
  return (
    (system === aOf(prep, k) && p >= d.aStart[k] && p < d.aEnd[k]) ||
    (system === bOf(prep, k) && p >= d.bStart[k] && p < d.bEnd[k])
  );
}

/** Whether the selection's elements hold a pair (few keys: read in order). */
function selected(prep: PairsPrep, e: Emphasis, k: number): boolean {
  const keys = e.elements?.keys;
  if (!keys) return false;
  const key = prep.data!.keys[k];
  for (let i = 0; i < keys.length; i++) if (keys[i].owner === prep.owner && keys[i].key === key) return true;
  return false;
}

/** Whether a pair's `a` end holds `x` and its `b` end `y`. */
function endsAt(d: PairsData, k: number, x: number, y: number): boolean {
  return x >= d.aStart[k] && x < d.aEnd[k] && y >= d.bStart[k] && y < d.bEnd[k];
}

/** Whether a pair has an end in a Location's ranges (a projected hover: a cell's two positions on linked sections). */
function endIn(prep: PairsPrep, k: number, loc: Location): boolean {
  const d = prep.data!;
  for (let i = 0; i < loc.parts.length; i++) {
    const part = loc.parts[i];
    const a = part.system === aOf(prep, k);
    const b = part.system === bOf(prep, k);
    if (!a && !b) continue;
    const r = part.ranges;
    for (let j = 0; j + 1 < r.length; j += 2)
      if ((a && d.aStart[k] < r[j + 1] && r[j] < d.aEnd[k]) || (b && d.bStart[k] < r[j + 1] && r[j] < d.bEnd[k]))
        return true;
  }
  return false;
}

/**
 * Whether a pair is lit in a pass. The hover: the hovered pair; as arcs, every pair with an end at a hovered position
 * (a base, or either position of a cell: its partners light) or in the hover projected here; as dots, the pair at the
 * hovered cell. The selection: a selected pair.
 */
function lit(lay: PairsLay, pass: PassName, e: Emphasis, k: number, hovered: number): boolean {
  const prep = lay.prep!;
  if (pass !== 'hover') return selected(prep, e, k);
  if (k === hovered) return true;
  const h = e.hover;
  if (h?.element) return false;
  const x = h?.position === undefined ? -1 : Math.floor(h.position);
  const y = h?.yPosition === undefined ? -1 : Math.floor(h.yPosition);
  if (lay.mode === 'dots') {
    if (!h || x < 0 || y < 0 || h.ySystem === undefined) return false;
    return (
      (h.system === prep.aSystem && h.ySystem === prep.bSystem && endsAt(prep.data!, k, x, y)) ||
      (prep.mirror && h.system === prep.bSystem && h.ySystem === prep.aSystem && endsAt(prep.data!, k, y, x))
    );
  }
  if (h && x >= 0 && endHolds(prep, k, h.system, x)) return true;
  if (h && y >= 0 && h.ySystem !== undefined && endHolds(prep, k, h.ySystem, y)) return true;
  const p = e.projected?.hover;
  return p !== undefined && endIn(prep, k, p);
}

/**
 * The hover: the hovered pair (its element), or the pairs at the hovered positions (`lit`), drawn over in the hover's
 * ink; the selection: each selected pair in the accent.
 */
function drawEmphasis(lay: PairsLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const prep = lay.prep!;
  const tokens = lay.tokens!;
  const h = e.hover;
  let hovered = -1;
  if (pass === 'hover') {
    if (!h && !e.projected?.hover) return;
    if (h?.element) {
      if (h.element.owner !== prep.owner) return;
      hovered = prep.index.get(h.element.key) ?? -1;
      if (hovered < 0) return;
    }
  } else if (!e.elements?.keys.length) return;
  const ink = pass === 'hover' ? tokens.colors['seq-hover'] : tokens.colors['msui-accent'];
  const width = pass === 'hover' ? lay.config!.outline.hover + 1 : lay.config!.outline.selection + 1;
  const g = out.args;
  if (lay.mode === 'arcs') {
    if (lay.clip) clipTo(out, lay.clip);
    for (let i = 0; i < lay.n; i++) {
      if (!lit(lay, pass, e, lay.pair[i], hovered)) continue;
      g[4] = NoColor;
      g[5] = ink;
      g[6] = width;
      g[7] = 0;
      out.curveArgs(lay.curve, 8 * i);
    }
    if (lay.clip) out.unclip();
    return;
  }
  for (let i = 0; i < lay.nd; i++) {
    if (!lit(lay, pass, e, lay.dotPair[i], hovered)) continue;
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = NoColor;
    g[5] = ink;
    g[6] = width - 1;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  const p = lay.seg;
  for (let i = 0; i < lay.ns; i++) {
    if (!lit(lay, pass, e, lay.segPair[i], hovered)) continue;
    p[0] = lay.segs[5 * i];
    p[1] = lay.segs[5 * i + 1];
    p[2] = lay.segs[5 * i + 2];
    p[3] = lay.segs[5 * i + 3];
    g[4] = NoColor;
    g[5] = ink;
    g[6] = lay.segs[5 * i + 4] + 2;
    g[7] = 0;
    out.pathArgs(p, 2, false);
  }
}

/**
 * Each pair's highlight channel, made again when the channels change (by identity: the viewer keeps the list): the
 * last channel whose element keys name the pair (a key names its first pair, as `index` does).
 */
function litOf(lay: PairsLay, highlights: readonly EmphasisHighlight[]): Int32Array {
  const prep = lay.prep!;
  let l = LitOf.get(lay);
  if (!l) LitOf.set(lay, (l = { for: undefined, prep: undefined, lit: new Int32Array(0) }));
  if (l.for === highlights && l.prep === prep) return l.lit;
  const n = prep.data?.n ?? 0;
  if (l.lit.length < n) l.lit = new Int32Array(n);
  const lit = l.lit;
  lit.fill(0);
  for (let c = 0; c < highlights.length; c++) {
    const keys = highlights[c].elements?.keys;
    if (!keys) continue;
    for (let i = 0; i < keys.length; i++) {
      if (keys[i].owner !== prep.owner) continue;
      const k = prep.index.get(keys[i].key);
      if (k !== undefined) lit[k] = c + 1;
    }
  }
  l.for = highlights;
  l.prep = prep;
  return lit;
}

/**
 * The highlight channels' pairs (_spec/README.md "Render", The pairs layer, Highlight): each pair a channel lists by its
 * key drawn over in the channel's color, as the selection's are in the accent (arcs, chords, dots and stretches).
 */
function drawHighlights(lay: PairsLay, out: PassWriter, highlights: readonly EmphasisHighlight[]) {
  if (!lay.prep?.data) return;
  let any = false;
  for (let c = 0; c < highlights.length && !any; c++) any = !!highlights[c].elements?.keys.length;
  if (!any) return;
  const lit = litOf(lay, highlights);
  const width = lay.config!.outline.selection + 1;
  const g = out.args;
  if (lay.mode === 'arcs') {
    if (lay.clip) clipTo(out, lay.clip);
    for (let i = 0; i < lay.n; i++) {
      const c = lit[lay.pair[i]];
      if (!c) continue;
      g[4] = NoColor;
      g[5] = highlights[c - 1].color;
      g[6] = width;
      g[7] = 0;
      out.curveArgs(lay.curve, 8 * i);
    }
    if (lay.clip) out.unclip();
    return;
  }
  for (let i = 0; i < lay.nd; i++) {
    const c = lit[lay.dotPair[i]];
    if (!c) continue;
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = NoColor;
    g[5] = highlights[c - 1].color;
    g[6] = width - 1;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  const p = lay.seg;
  for (let i = 0; i < lay.ns; i++) {
    const c = lit[lay.segPair[i]];
    if (!c) continue;
    p[0] = lay.segs[5 * i];
    p[1] = lay.segs[5 * i + 1];
    p[2] = lay.segs[5 * i + 2];
    p[3] = lay.segs[5 * i + 3];
    g[4] = NoColor;
    g[5] = highlights[c - 1].color;
    g[6] = lay.segs[5 * i + 4] + 2;
    g[7] = 0;
    out.pathArgs(p, 2, false);
  }
}

/** The pairs drawn (the arcs, the dots) as items: a pair's label (else its key) and its first end's extent. */
function pairItems(lay: PairsLay): ItemDescription[] {
  const d = lay.prep?.data;
  if (!d) return [];
  const out: ItemDescription[] = [];
  const add = (k: number, x0: number, x1: number) => {
    const key = d.keys[k];
    out.push({
      key,
      label: d.labels?.[k] || String(key),
      start: d.aStart[k],
      end: d.aEnd[k],
      x0: Math.min(x0, x1),
      x1: Math.max(x0, x1),
    });
  };
  if (lay.mode === 'arcs') for (let i = 0; i < lay.n; i++) add(lay.pair[i], lay.curve[8 * i], lay.curve[8 * i + 6]);
  else if (lay.mode === 'dots')
    for (let i = 0; i < lay.nd; i++) add(lay.dotPair[i], lay.dots[4 * i], lay.dots[4 * i] + lay.dots[4 * i + 2]);
  return out;
}
