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
// end, never a dangling arc; on one line of straight strips (a track, a row of blocks) an end past the window is taken
// at most `reach` (a share of the window, 0.25) past the edge it lies beyond, so an arc moves continuously on pan. An
// arc's height is by its span against the widest within reach (nested arcs nest, the widest fills the band); `side`
// above or below (a ring's below reaches the hole: chords). Crossing pairs are drawn crossing (`same`), or `split`: the
// largest nested set on `side`, the rest on the other side (exact to 2,000 pairs, greedy past it). At most
// `pairs.maxArcs` arcs: past it one a bin of the ends' device columns (widest first), the rest a density strip of ends.
// A hit is the arc or caret within `pairs.hitPx` of the point: the pair's element, its extent both ends. In a grid (a
// strip with `ctx.axisY`) the layer is laid out a strip at a time instead: a dot a pair of single positions, a diagonal
// for a stretch, at (a, b) and with `mirror` at (b, a) too; kinds at one cell side by side.

import { issueOf } from '@molstar/core/repr/tree';
import type { ChannelName } from '@molstar/seqx-core/theme/index';
import { channelLegends } from './channels.js';
import { grown, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { type PairEnds } from './pairs-ends.js';
import { type Emphasis, type Hit, type LayerType, type LayoutContext, sectionId } from './layer-type.js';
import type { AxisTransform } from './axis.js';
import { WarpKind } from './layout/strips.js';
import {
  Caret,
  clipTo,
  DefaultColor,
  endsIn,
  HitSamples,
  Inset,
  MinArc,
  pairColor,
  pairEnv,
  pairHit,
  PairsLay,
  type PairsPrep,
  prepare,
  segmentDistance2,
  type Side,
} from './pairs-lay.js';
import { drawDots, hitDots, layoutDots } from './pairs-dots.js';
import { heldBins } from './pairs-held.js';
import { regionPairs } from './pairs-region.js';
import { drawEmphasis, drawHighlights, pairItems } from './pairs-emphasis.js';
import { PlacePoint, type Placement } from './layout/types.js';
import { withAlpha } from './tokens.js';

export { PairsLay, type PairsPrep } from './pairs-lay.js';

export { crosses, SplitExact, splitCrossing } from './pairs-split.js';

const PA = new PlacePoint();
const PB = new PlacePoint();
const PM = new PlacePoint();

export const PairsType: LayerType<PairsPrep, PairsLay> = {
  name: 'pairs',
  fallback: 'ranges',
  scope: 'section',
  prepare,

  layout(prep, ctx: LayoutContext, data, _overlay, prev?: PairsLay): PairsLay {
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
    // a source's pairs: this frame's read (bins, or items once they fit) as the pairs laid out
    if (prep.region && !prep.data) {
      const read = ctx.place ? regionPairs(prep, prep.region, ctx, data, lay) : undefined;
      if (!read?.data) return lay;
      prep = read;
      lay.prep = prep;
    }
    if (!prep.data) return lay;
    // held pairs past the gate: bin pairs at the pairs level a source's are read at (`pairs-held.ts`)
    if (!prep.region && ctx.place) {
      const binned = heldBins(prep, ctx, data);
      if (binned) {
        prep = binned;
        lay.prep = prep;
      }
    }
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

/** `virtualX`'s result (no number returned: a double returned from a call not inlined is boxed). */
const VX = new Float64Array(1);

/**
 * The page x of an end past the window on one line (unit `u` outside `[u0, u1]`) into `VX[0]`: past the edge it lies
 * beyond by its distance at the window's px a unit, at most `reach` of the window's width. Continuous as a pan moves
 * the end across the edge (it meets the end's placed x there) and as it goes past reach.
 */
function virtualX(ax: AxisTransform, u: number, reach: number) {
  const ppu = ax.width / (ax.u1 - ax.u0);
  const m = reach * ax.width;
  VX[0] = u < ax.u0 ? ax.x0 - Math.min(m, (ax.u0 - u) * ppu) : ax.x0 + ax.width + Math.min(m, (u - ax.u1) * ppu);
}

function layoutArcs(lay: PairsLay, prep: PairsPrep, ctx: LayoutContext, place: Placement) {
  const s = place.strips;
  if (!s.n) return;
  lay.mode = 'arcs';
  const band = ctx.band ?? 0;
  const ends = endsIn(prep, ctx.section, ctx.axis);
  const units = ends.u;
  const nb = s.bands;
  // one line of translated strips (a track, a row of blocks): an end past the window is taken at most `reach` of the
  // window past the edge its unit lies beyond (`virtualX`), on its partner's strip
  let line = true;
  for (let i = 0; i < s.n; i++) if (s.warp[i] !== WarpKind.translate || s.top[i] !== s.top[0]) line = false;
  const ax = ctx.axis;
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
    const oa = units[2 * e];
    const ob = units[2 * e + 1];
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
    if (line) {
      // past the window: its unit outside it (an end inside it on no strip, a block too narrow to read, stays a caret)
      const outA = !Number.isNaN(oa) && Number.isNaN(xa) && (oa < ax.u0 || oa > ax.u1);
      const outB = !Number.isNaN(ob) && Number.isNaN(xb) && (ob < ax.u0 || ob > ax.u1);
      if (outA && outB && oa < ax.u0 === ob < ax.u0) {
        // both past one edge: not drawn, but its span counts in the scale while within reach, so a pair entering the
        // view brings no snap
        virtualX(ax, oa, prep.reach);
        const va = VX[0];
        virtualX(ax, ob, prep.reach);
        const span = Math.abs(VX[0] - va);
        if (span > widest) widest = span;
        continue;
      }
      // the ends on one strip's frame (one line: the strips' frames differ by x alone), the arc symmetric between them
      const home = !outA && !Number.isNaN(xa) ? sa : !outB && !Number.isNaN(xb) ? sb : 0;
      if (outA) {
        virtualX(ax, oa, prep.reach);
        xa = VX[0] - s.dx[home];
      } else if (!Number.isNaN(xa)) xa += s.dx[sa] - s.dx[home];
      if (outB) {
        virtualX(ax, ob, prep.reach);
        xb = VX[0] - s.dx[home];
      } else if (!Number.isNaN(xb)) xb += s.dx[sb] - s.dx[home];
      if (!Number.isNaN(xa)) sa = home;
      if (!Number.isNaN(xb)) sb = home;
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
  // a ring says it after its chords (the count of those it drops); a line now, its ends a density strip
  const ring = s.warp[0] === WarpKind.polar;
  if (binned) {
    if (!ring)
      ctx.report(
        issueOf(RenderIssues, 'seq.render.pairs.limit', `More than ${max} pairs in view: binned`, {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.pairs.limit'].hint,
        }),
      );
    order = new Int32Array(m);
    for (let i = 0; i < m; i++) order[i] = i;
    // bins: the strongest first (by count, ties by span); items: the widest first
    const w = prep.weight;
    if (w) order.sort((x, y) => w[Raw[4 * y + 3]] - w[Raw[4 * x + 3]] || Spans[y] - Spans[x]);
    else order.sort((x, y) => Spans[y] - Spans[x]);
    densityStrip(lay, place, band, ctx.dpr);
  }
  const seen = binned && !prep.weight ? new Set<number>() : undefined;
  // a bin pair's stroke by its count: `pairs.binStroke` CSS px on a square-root scale over the counts in view
  let maxWeight = 0;
  if (prep.weight) for (let i = 0; i < m; i++) maxWeight = Math.max(maxWeight, prep.weight[Raw[4 * i + 3]]);
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
    } else if (binned && lay.n >= max) {
      // bin pairs past `maxArcs`: the strongest kept (they come first), the rest the density strip's
      addDensity(lay, xa + s.dx[sa], dpr);
      addDensity(lay, xb + s.dx[sb], dpr);
      continue;
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
    lay.width = grown(lay.width, lay.n + 1);
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
    const stroke = ctx.config.pairs.binStroke;
    // a bin pair by its count: its stroke within `binStroke`, its colour's alpha from a quarter to whole, so the many
    // weak chords of a genome ring are a haze and the strong ones read (region-reads Default 3)
    const share = prep.weight && maxWeight > 0 ? Math.sqrt(prep.weight[k] / maxWeight) : 1;
    const color = pairColor(lay, prep, k);
    lay.color[lay.n] = prep.weight ? withAlpha(color, ((color & 0xff) / 255) * (0.25 + 0.75 * share)) : color;
    lay.width[lay.n] = prep.weight ? stroke.min + (stroke.max - stroke.min) * share : 1;
    if (lay.extents.any) lay.extents.add(prep.datum!.set(k));
    lay.n++;
  }
  if (binned && ring)
    ctx.report(
      issueOf(
        RenderIssues,
        'seq.render.pairs.limit',
        `More than ${max} pairs in view: the ${lay.n} strongest are drawn on the ring, ${m - lay.n} are not`,
        { address: ctx.key, hint: RenderIssues.codes['seq.render.pairs.limit'].hint },
      ),
    );
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
  // a ring has no straight strip: the chords past the cap are dropped, never drawn as ticks along a line above it
  if (s.warp[0] === WarpKind.polar) {
    lay.densityN = 0;
    lay.densityMax = 0;
    return;
  }
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

function drawArcs(lay: PairsLay, out: PassWriter) {
  const g = out.args;
  if (lay.clip) clipTo(out, lay.clip);
  for (let i = 0; i < lay.n; i++) {
    g[4] = NoColor;
    g[5] = lay.color[i];
    g[6] = lay.width[i];
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
