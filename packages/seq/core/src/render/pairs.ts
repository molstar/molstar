/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs layer (_spec/README.md "Render", The pairs layer): two locations that belong together, drawn by where its
// section puts them. Laid out once a section (`scope: 'section'`): each end placed with the placement's `place`, its
// arc's control the mid offset at the arc's height in its strip's frame through the strip's warp, so one curve is an
// arc over a track or a block, an arc between two lines, a chord through a ring's hole, with no layout named here. An
// end not placed leaves a caret at the placed end, never a dangling arc; on one straight strip an end past the window is
// taken `reach` (a share of the window, 0.25) past its edge, so an arc moves continuously on pan. An arc's height is by
// its span against the widest drawn (nested arcs nest, the widest fills the band); `side` above or below (a ring's
// below reaches the hole: chords). Crossing pairs are drawn crossing (`same`), or `split`: the largest nested set on
// `side`, the rest on the other side (exact to 2,000 pairs, greedy past it). At most `pairs.maxArcs` arcs: past it one a
// bin of the ends' device columns (widest first), the rest a density strip of ends. A hit is the arc within
// `pairs.hitPx` of the point: the pair's element, its extent both ends. In a grid (a strip with `ctx.axisY`) the layer
// is laid out a strip at a time instead: a dot a pair of single positions, a diagonal for a stretch, at (a, b) and with
// `mirror` at (b, a) too; kinds at one cell side by side.

import { issueOf } from '@molstar/core/repr/tree';
import { type ElementKey, type Location, location } from '@molstar/seqx-core/query';
import type { PairsData, StaticModel } from '@molstar/seqx-core/resolve';
import { cellOffset, type ResolvedLayer, type ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import type { RenderConfig } from './config.js';
import { NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { type Box, type Emphasis, type Hit, type LayerType, type LayoutContext, sectionId } from './layer-type.js';
import { WarpKind } from './layout/strips.js';
import { PlacePoint, type Placement } from './layout/types.js';
import type { SeqTokens } from './tokens.js';

/** Pairs a `split` places exactly (the largest nested set: quadratic); past it, greedily (fixed: no figure asks more). */
export const SplitExact = 2000;

/** An arc's least height, px: never flat. */
const MinArc = 3;

/** The points a curve is sampled at for a hit. */
const HitSamples = 16;

/** The caret's half width at a placed end whose partner is not placed, px. */
const Caret = 3;

type Side = 'above' | 'below';

export interface PairsPrep {
  readonly owner: string;
  readonly data?: PairsData;
  readonly aSystem: string;
  readonly bSystem: string;
  readonly side: Side;
  readonly reach: number;
  readonly acrossLines: boolean;
  readonly mirror: boolean;
  readonly triangle: 'both' | 'upper' | 'lower';
  /** The `color` param when uniform (a token or a color), else each kind its category's. */
  readonly color?: string | number;
  /** With `split`: 1 where a pair is drawn on the other side. */
  readonly other?: Uint8Array;
  /** Each pair's index by its key. */
  readonly index: ReadonlyMap<string | number, number>;
  /** Each pair's element and extent (both ends), made when first hit. */
  readonly elements: (ElementKey | undefined)[];
  readonly extents: (Location | undefined)[];
  /** Each pair's two ends' centres as display units of a section (NaN: not on it), by section. */
  readonly units: WeakMap<ResolvedSection, Float64Array>;
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

function prepare(layer: ResolvedLayer, model: StaticModel): PairsPrep {
  const owner = layer.data?.address ?? '';
  const a = model.annotations.get(owner);
  const data = a?.shape === 'pairs' ? a.pairs : undefined;
  const p = layer.params;
  const index = new Map<string | number, number>();
  if (data) for (let k = 0; k < data.n; k++) if (!index.has(data.keys[k])) index.set(data.keys[k], k);
  const t = p.triangle;
  const uniform = layer.theme.color;
  return {
    owner,
    ...(data ? { data } : {}),
    aSystem: a?.target ?? '',
    bSystem: a?.ySystem ?? '',
    side: p.side === 'below' ? 'below' : 'above',
    reach: typeof p.reach === 'number' && p.reach >= 0 ? p.reach : 0.25,
    acrossLines: p.across_lines === true,
    mirror: p.mirror === true,
    triangle: t === 'upper' || t === 'lower' ? t : 'both',
    ...(typeof uniform === 'string' || typeof uniform === 'number' ? { color: uniform } : {}),
    ...(data && p.crossing === 'split' && a?.target === a?.ySystem ? { other: splitCrossing(data) } : {}),
    index,
    elements: [],
    extents: [],
    units: new WeakMap(),
  };
}

/**
 * The pairs drawn on the other side under `split`: the complement of the largest nested (non-crossing) set of the
 * pairs as intervals between their ends' centres; exact by the interval recurrence to `SplitExact` pairs, greedy past
 * it (narrowest first, kept when it crosses none kept).
 */
export function splitCrossing(d: PairsData): Uint8Array {
  const n = d.n;
  const lo = new Float64Array(n);
  const hi = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const ca = (d.aStart[k] + d.aEnd[k]) / 2;
    const cb = (d.bStart[k] + d.bEnd[k]) / 2;
    lo[k] = Math.min(ca, cb);
    hi[k] = Math.max(ca, cb);
  }
  const keep = new Uint8Array(n);
  if (n <= SplitExact) nestedExact(n, lo, hi, keep);
  else nestedGreedy(n, lo, hi, keep);
  const other = new Uint8Array(n);
  for (let k = 0; k < n; k++) other[k] = keep[k] ? 0 : 1;
  return other;
}

/** Whether two intervals cross: one's ends either side of one end of the other. */
export function crosses(alo: number, ahi: number, blo: number, bhi: number): boolean {
  return (alo < blo && blo < ahi && ahi < bhi) || (blo < alo && alo < bhi && bhi < ahi);
}

/**
 * The largest set of non-crossing intervals (nested or apart): over the 2n ends in order, F(i, j) the most held within
 * ends i..j, F(i, j) = max(F(i + 1, j), 1 + F(i + 1, k − 1) + F(k + 1, j)) where k is i's partner in (i, j]. At one
 * place closing ends come first and a wider interval outside a narrower one, as `crosses` reads them. O(n²) time, a
 * triangle of int16.
 */
function nestedExact(n: number, lo: Float64Array, hi: Float64Array, keep: Uint8Array) {
  const m = 2 * n;
  const ends = new Int32Array(m);
  for (let i = 0; i < m; i++) ends[i] = i;
  const at = (e: number) => (e < n ? lo[e] : hi[e - n]);
  ends.sort((x, y) => {
    const d = at(x) - at(y);
    if (d) return d;
    const xo = x < n;
    const yo = y < n;
    if (xo !== yo) return xo ? 1 : -1;
    const wx = hi[x % n] - lo[x % n];
    const wy = hi[y % n] - lo[y % n];
    return xo ? wy - wx : wx - wy;
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

/** The centre of `[start, end)` on an axis as a display unit, NaN when not shown. */
function centreOf(axis: Parameters<typeof cellOffset>[0], start: number, end: number): number {
  const a = cellOffset(axis, start);
  const b = cellOffset(axis, end - 1);
  if (a === undefined || b === undefined) return Number.NaN;
  return (Math.min(a, b) + Math.max(a, b) + 1) / 2;
}

/**
 * Each pair's ends' centres as the section's display units (NaN: not on it): on the section's axis, or on the first
 * block on the end's system (merged copies are one block: one arc a copy group).
 */
function unitsOn(prep: PairsPrep, section: ResolvedSection): Float64Array {
  let u = prep.units.get(section);
  if (u) return u;
  const d = prep.data!;
  u = new Float64Array(2 * d.n).fill(Number.NaN);
  const on = (system: string, s: number, e: number): number => {
    const blocks = section.blocks;
    if (!blocks) return section.axis.system === system ? centreOf(section.axis, s, e) : Number.NaN;
    for (const b of blocks)
      if (b.axis.system === system && b.offset !== undefined) {
        const c = centreOf(b.axis, s, e);
        if (c === c) return b.offset + c;
      }
    return Number.NaN;
  };
  for (let k = 0; k < d.n; k++) {
    u[2 * k] = on(prep.aSystem, d.aStart[k], d.aEnd[k]);
    u[2 * k + 1] = on(prep.bSystem, d.bStart[k], d.bEnd[k]);
  }
  prep.units.set(section, u);
  return u;
}

const CategoryTokens = [
  'seq-cat-0',
  'seq-cat-1',
  'seq-cat-2',
  'seq-cat-3',
  'seq-cat-4',
  'seq-cat-5',
  'seq-cat-6',
  'seq-cat-7',
] as const;

/** A pair's color: the uniform param's, else its kind's category when there are several (else the ranges' ink). */
function pairColor(prep: PairsPrep, tokens: SeqTokens, k: number): number {
  const c = prep.color;
  if (typeof c === 'number') return c;
  if (typeof c === 'string' && c.startsWith('token:')) {
    const t = tokens.colors[c.slice(6) as keyof SeqTokens['colors']];
    if (t !== undefined) return t;
  }
  const kind = prep.data!.kind[k];
  if (kind < 0 || prep.data!.kinds.length < 2) return tokens.colors['seq-range'];
  return tokens.colors[CategoryTokens[kind % CategoryTokens.length]];
}

/** The point a hit tests (x, y). */
const Pt = new Float64Array(2);

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
    if (!prep.data) return lay;
    if (ctx.axisY) layoutDots(lay, prep, ctx, ctx.axis, ctx.axisY);
    else if (ctx.place) layoutArcs(lay, prep, ctx, ctx.place);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    if (!lay.tokens || !lay.config || !lay.prep) return;
    if (pass === 'data') {
      if (lay.mode === 'arcs') drawArcs(lay, out);
      else if (lay.mode === 'dots') drawDots(lay, out);
      return;
    }
    if (emphasis && (pass === 'hover' || pass === 'selection')) drawEmphasis(lay, pass, out, emphasis);
  },

  hit(lay, x, y) {
    // the point through a scratch: a fractional number passed to a call that is not inlined is boxed
    Pt[0] = x;
    Pt[1] = y;
    if (lay.mode === 'arcs') return hitArcs(lay);
    if (lay.mode === 'dots') return hitDots(lay);
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
  if (s.warp[i] === WarpKind.polar) return s.axes[i].xOfUnit(u - s.first[i]);
  return pt.x - s.dx[i];
}

function grow(lay: PairsLay, n: number) {
  if (n <= lay.pair.length) return;
  const m = Math.max(n, 2 * lay.pair.length);
  const pair = new Int32Array(m);
  pair.set(lay.pair);
  lay.pair = pair;
  const color = new Uint32Array(m);
  color.set(lay.color);
  lay.color = color;
  const curve = new Float64Array(8 * m);
  curve.set(lay.curve);
  lay.curve = curve;
}

/** What a layout of arcs reads a pair at a time (kept: grown, never shrunk). */
let Spans = new Float64Array(64);
let StripOf = new Int32Array(128);
let Raw = new Float64Array(256);

function layoutArcs(lay: PairsLay, prep: PairsPrep, ctx: LayoutContext, place: Placement) {
  const d = prep.data!;
  const s = place.strips;
  if (!s.n) return;
  lay.mode = 'arcs';
  const band = ctx.band ?? 0;
  const units = unitsOn(prep, ctx.section);
  const nb = s.bands;
  const one = s.n === 1 && s.warp[0] === WarpKind.translate;
  const a0 = s.axes[0];
  const w = a0.u1 - a0.u0;
  const reachLo = a0.u0 - prep.reach * w;
  const reachHi = a0.u1 + prep.reach * w;
  if (Spans.length < d.n) {
    Spans = new Float64Array(2 * d.n);
    StripOf = new Int32Array(4 * d.n);
    Raw = new Float64Array(8 * d.n);
  }
  let m = 0;
  let widest = 0;
  let unplaced = false;
  for (let k = 0; k < d.n; k++) {
    let oa = units[2 * k];
    let ob = units[2 * k + 1];
    if (oa !== oa && ob !== ob) {
      unplaced = true;
      continue;
    }
    let xa = oa === oa ? frameX(place, oa, PA) : Number.NaN;
    let xb = ob === ob ? frameX(place, ob, PB) : Number.NaN;
    let sa = PA.strip;
    let sb = PB.strip;
    if (one) {
      // one straight strip: an end past the window taken `reach` past its edge
      const outA = oa === oa && xa !== xa;
      const outB = ob === ob && xb !== xb;
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
    if (xa !== xa && xb !== xb) continue;
    // an end not placed (a chain not shown, another line without `across_lines`): a caret at the placed one
    const across = xa === xa && xb === xb && s.top[sa] !== s.top[sb] && !prep.acrossLines;
    if (xa !== xa || xb !== xb || across) {
      if (xa === xa) caret(lay, prep, place, sa, band, xa, k);
      if (xb === xb) caret(lay, prep, place, sb, band, xb, k);
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
    const base = side === 'above' ? top + h : top;
    const frac = widest > 0 ? Spans[i] / widest : 1;
    // the control's offset from the base: twice the apex's height (a curve's apex is halfway to its control); on a
    // ring below, toward the centre by the span (a chord through the hole)
    const depth =
      polar && side === 'below'
        ? (place.warpOf!(sa).r0 - base) * Math.min(1, frac) * 1.5
        : 2 * Math.max(MinArc, (h - 2) * frac);
    const cy = side === 'above' ? base - depth : base + depth;
    grow(lay, lay.n + 1);
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
      const mid = polar ? s.axes[sa].xOfUnit(Raw[4 * i + 2] - s.first[sa]) : (xa + xb) / 2;
      toPage(place, sa, mid, cy, c, o + 2);
    } else {
      const mx = frameX(place, Raw[4 * i + 2], PM);
      if (mx === mx) toPage(place, PM.strip, mx, cy, c, o + 2);
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
    lay.color[lay.n] = pairColor(prep, ctx.tokens, k);
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

function curveDistance(c: Float64Array, o: number): void {
  const x = Pt[0];
  const y = Pt[1];
  let best = Number.POSITIVE_INFINITY;
  let px = c[o];
  let py = c[o + 1];
  for (let s = 1; s <= HitSamples; s++) {
    const t = s / HitSamples;
    const u = 1 - t;
    const qx = u * u * u * c[o] + 3 * u * u * t * c[o + 2] + 3 * u * t * t * c[o + 4] + t * t * t * c[o + 6];
    const qy = u * u * u * c[o + 1] + 3 * u * u * t * c[o + 3] + 3 * u * t * t * c[o + 5] + t * t * t * c[o + 7];
    // the segment's nearest point
    const dx = qx - px;
    const dy = qy - py;
    const l2 = dx * dx + dy * dy;
    let f = l2 > 0 ? ((x - px) * dx + (y - py) * dy) / l2 : 0;
    f = f < 0 ? 0 : f > 1 ? 1 : f;
    const ex = px + f * dx - x;
    const ey = py + f * dy - y;
    const d2 = ex * ex + ey * ey;
    if (d2 < best) best = d2;
    px = qx;
    py = qy;
  }
  Dist[0] = best;
}

/** The arc nearest a page point within `pairs.hitPx` (a tie: the later, drawn on top). */
function hitArcs(lay: PairsLay): Hit | undefined {
  const x = Pt[0];
  const y = Pt[1];
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
    curveDistance(c, o);
    const d = Dist[0];
    if (d <= bestD) {
      bestD = d;
      best = i;
    }
  }
  return best < 0 ? undefined : pairHit(lay, lay.pair[best]);
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

/** A pair's hit: its element and its extent, both ends (made once a pair). */
function pairHit(lay: PairsLay, k: number): Hit {
  const prep = lay.prep!;
  const d = prep.data!;
  let element = prep.elements[k];
  if (!element) element = prep.elements[k] = { owner: prep.owner, key: d.keys[k] };
  let extent = prep.extents[k];
  if (!extent) {
    const parts =
      prep.aSystem === prep.bSystem
        ? [{ system: prep.aSystem, ranges: [d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k]] }]
        : [
            { system: prep.aSystem, ranges: [d.aStart[k], d.aEnd[k]] },
            { system: prep.bSystem, ranges: [d.bStart[k], d.bEnd[k]] },
          ];
    extent = prep.extents[k] = location({ parts });
  }
  const h = lay.scratch;
  h.layer = lay.key;
  h.section = lay.section;
  h.system = prep.aSystem;
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
    if (forward) dot(lay, prep, a, ay, k, d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k]);
    if (prep.mirror && back) dot(lay, prep, a, ay, k, d.bStart[k], d.bEnd[k], d.aStart[k], d.aEnd[k]);
  }
  sideBySide(lay);
}

function growDots(lay: PairsLay) {
  const dots = new Float64Array(2 * lay.dots.length);
  dots.set(lay.dots);
  lay.dots = dots;
  const p = new Int32Array(2 * lay.dotPair.length);
  p.set(lay.dotPair);
  lay.dotPair = p;
  const c = new Uint32Array(2 * lay.dotColor.length);
  c.set(lay.dotColor);
  lay.dotColor = c;
}

function growSegs(lay: PairsLay) {
  const segs = new Float64Array(2 * lay.segs.length);
  segs.set(lay.segs);
  lay.segs = segs;
  const p = new Int32Array(2 * lay.segPair.length);
  p.set(lay.segPair);
  lay.segPair = p;
  const c = new Uint32Array(2 * lay.segColor.length);
  c.set(lay.segColor);
  lay.segColor = c;
}

/** A pair's dot (single positions) or diagonal (a stretch, antiparallel: x's first with y's last), kept by the triangle. */
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
) {
  const t = prep.triangle;
  const cx = (xs + xe - 1) / 2;
  const cy = (ys + ye - 1) / 2;
  if ((t === 'upper' && !(cy < cx)) || (t === 'lower' && !(cy > cx))) return;
  const color = pairColor(prep, lay.tokens!, k);
  if (xe - xs === 1 && ye - ys === 1) {
    if (!a.cell(xs)) return;
    const ax0 = a.sx0;
    const ax1 = a.sx1;
    if (!ay.cell(ys)) return;
    const i = lay.nd++;
    if (4 * i + 4 > lay.dots.length || i >= lay.dotPair.length) growDots(lay);
    // at least two device pixels a side
    const m = 2 / a.dpr;
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
  // a stretch: from (x's first, y's last) to (x's last, y's first), through the cells' centres
  if (!a.cell(xs)) return;
  const x0 = (a.sx0 + a.sx1) / 2;
  if (!a.cell(xe - 1)) return;
  const x1 = (a.sx0 + a.sx1) / 2;
  if (!ay.cell(ye - 1)) return;
  const y0 = (ay.sx0 + ay.sx1) / 2;
  if (!ay.cell(ys)) return;
  const y1 = (ay.sx0 + ay.sx1) / 2;
  const i = lay.ns++;
  if (5 * i + 5 > lay.segs.length || i >= lay.segPair.length) growSegs(lay);
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
  for (let i = 0; i < lay.nd; i++) {
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = lay.dotColor[i];
    g[5] = NoColor;
    g[6] = 0;
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
function hitDots(lay: PairsLay): Hit | undefined {
  const x = Pt[0];
  const y = Pt[1];
  const d = lay.dots;
  for (let i = lay.nd - 1; i >= 0; i--)
    if (x >= d[4 * i] && x < d[4 * i] + d[4 * i + 2] && y >= d[4 * i + 1] && y < d[4 * i + 1] + d[4 * i + 3])
      return pairHit(lay, lay.dotPair[i]);
  const r = lay.config!.pairs.hitPx;
  const s = lay.segs;
  for (let i = lay.ns - 1; i >= 0; i--) {
    const x0 = s[5 * i];
    const y0 = s[5 * i + 1];
    const dx = s[5 * i + 2] - x0;
    const dy = s[5 * i + 3] - y0;
    const l2 = dx * dx + dy * dy;
    let f = l2 > 0 ? ((x - x0) * dx + (y - y0) * dy) / l2 : 0;
    f = f < 0 ? 0 : f > 1 ? 1 : f;
    const ex = x0 + f * dx - x;
    const ey = y0 + f * dy - y;
    const w = Math.max(r, s[5 * i + 4] / 2);
    if (ex * ex + ey * ey <= w * w) return pairHit(lay, lay.segPair[i]);
  }
  return undefined;
}

// ---- emphasis ----

/** Whether a pair has an end holding position `p` on `system` (a hovered base: its partner lights). */
function endHolds(prep: PairsPrep, k: number, system: string, p: number): boolean {
  const d = prep.data!;
  return (
    (system === prep.aSystem && p >= d.aStart[k] && p < d.aEnd[k]) ||
    (system === prep.bSystem && p >= d.bStart[k] && p < d.bEnd[k])
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

/** Whether a pair is lit in a pass: the hovered one or one with an end at the hovered position; a selected one. */
function lit(prep: PairsPrep, pass: PassName, e: Emphasis, k: number, hovered: number, at: number): boolean {
  if (pass !== 'hover') return selected(prep, e, k);
  return k === hovered || (at >= 0 && endHolds(prep, k, e.hover!.system, at));
}

/**
 * The hover: the hovered pair (its element), or every pair with an end at the hovered position (a base: its partner
 * lights), drawn over in the hover's ink; the selection: each selected pair in the accent.
 */
function drawEmphasis(lay: PairsLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const prep = lay.prep!;
  const tokens = lay.tokens!;
  const h = e.hover;
  let hovered = -1;
  let at = -1;
  if (pass === 'hover') {
    if (!h) return;
    if (h.element) {
      if (h.element.owner !== prep.owner) return;
      hovered = prep.index.get(h.element.key) ?? -1;
    } else if (h.position !== undefined && h.yPosition === undefined) at = Math.floor(h.position);
    if (hovered < 0 && at < 0) return;
  } else if (!e.elements?.keys.length) return;
  const ink = pass === 'hover' ? tokens.colors['seq-hover'] : tokens.colors['msui-accent'];
  const width = pass === 'hover' ? lay.config!.outline.hover + 1 : lay.config!.outline.selection + 1;
  const g = out.args;
  if (lay.mode === 'arcs') {
    if (lay.clip) clipTo(out, lay.clip);
    for (let i = 0; i < lay.n; i++) {
      if (!lit(prep, pass, e, lay.pair[i], hovered, at)) continue;
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
    if (!lit(prep, pass, e, lay.dotPair[i], hovered, at)) continue;
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
    if (!lit(prep, pass, e, lay.segPair[i], hovered, at)) continue;
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
