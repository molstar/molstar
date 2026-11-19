/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rings layout (_spec/README.md "Render", Rings layout): a section's tracks as concentric rings, the first the
// outermost, round its whole axis; a strip a block, unrolled along the outer ring (the frame the outer circumference
// wide, so a layer reduces per angular column there) and bent onto the page by one polar warp (x to an angle clockwise
// from `rotation`, y to a radius inward). Blocks are arcs with `rings.blockGap` between them (a whole cyclic block the
// full circle, its seam no break). Items' labels go to a sink and are placed outside when the ring is laid out (never
// on a window move or a hover): sorted by angle, spread a line apart per side, each with a leader out along its angle;
// those a side cannot hold dropped by the layer's order, then the shorter item. The hole holds the title. The window
// (what a follower shows, `render/window.ts`) is a sector over the rings with its grips on the rim, out of the ring's
// key: a pan of its follower draws the ring from its kept lays.

import { piecesTotal, type ResolvedSection } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from '../axis.js';
import type { RenderConfig } from '../config.js';
import { fontSize } from '../config.js';
import type { ElementKey } from '@molstar/seqx-core/query';
import { Anchor, grown, NoColor, type PassWriter } from '../draw-lists.js';
import type { Box, Emphasis, Hit, LabelSink } from '../layer-type.js';
import { sectionId } from '../layer-type.js';
import type { TextMetrics } from '../text.js';
import type { SeqTokens } from '../tokens.js';
import type { SectionLayout } from '../tracks.js';
import { drawWindowArgs, WindowMetrics, windowPart } from '../window.js';
import { PlacePoint, Strips, WarpKind } from './strips.js';
import { clampWindow, TracksNav } from './tracks.js';
import type { NavContext, Placement, PlaceContext, SectionLayoutType, SectionNav } from './types.js';
import { PolarWarp } from './warp.js';

// anchors read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;
const AnchorEnd: number = Anchor.end;
const AnchorMiddle: number = Anchor.middle;

const TwoPi = 2 * Math.PI;

/**
 * Fixed sizes (CSS px, as other hit boxes are): the room a side keeps for outside labels (they fade past it), the
 * window's rim outside the outer ring (its body and grips are hit there, so the data under the window stays hit), the
 * margin a ring keeps when it has no labels.
 */
export const RingMetrics = { labelRoom: 96, rim: 12, margin: 4 } as const;

/** A rings section's params as the spec writes them (_spec/README.md "Resolved view"), with their defaults. */
export interface RingParams {
  /** Where the origin is, radians clockwise from the top. */
  readonly rotation: number;
  readonly labels: 'outside' | 'none';
}

export function ringParams(params: ResolvedSection['params']): RingParams {
  const r = params.rotation;
  return {
    rotation: typeof r === 'number' && Number.isFinite(r) ? (r * Math.PI) / 180 : 0,
    labels: params.labels === 'none' ? 'none' : 'outside',
  };
}

/** A ring's outside labels: given while its strips are laid out, placed after; their arrays kept frame to frame. */
class RingLabels {
  /** Given: n of them. */
  n = 0;
  text: string[] = [];
  /** The label's angle (radians, clockwise from the top, in [0, 2π)), its item's length, its order (band, then layer). */
  angle = new Float64Array(16);
  length = new Float64Array(16);
  order = new Float64Array(16);
  hits: (Hit | undefined)[] = [];
  /** Placed: m of them, by index into the given; each its text's anchor point, side, width and its leader's points. */
  m = 0;
  placed = new Int32Array(16);
  x = new Float64Array(16);
  y = new Float64Array(16);
  right = new Uint8Array(16);
  width = new Float64Array(16);
  /** A leader: from the ring out along its angle, to its bend, to the label (x, y a point). */
  leader = new Float64Array(96);

  reset() {
    this.n = 0;
    this.m = 0;
  }

  add(text: string, angle: number, length: number, order: number, hit: Hit | undefined) {
    const i = this.n++;
    this.angle = grown(this.angle, i + 1);
    this.length = grown(this.length, i + 1);
    this.order = grown(this.order, i + 1);
    this.text[i] = text;
    this.angle[i] = angle;
    this.length[i] = length;
    this.order[i] = order;
    this.hits[i] = hit;
  }

  ensurePlaced(n: number) {
    this.placed = grown(this.placed, n);
    this.x = grown(this.x, n);
    this.y = grown(this.y, n);
    this.right = grown(this.right, n);
    this.width = grown(this.width, n);
    this.leader = grown(this.leader, 6 * n);
  }
}

/**
 * A ring laid out: its center and outer radius, its bands' radii, one polar warp for all its strips (their frames
 * share the unrolled circle: x is the arc length along the outer ring from the origin), its labels, its title and its
 * window. `at` brings a page point back through the warp (`atan2` and the radius), `place` sends a unit out.
 */
export class RingPlacement implements Placement {
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  contentHeight = 0;
  readonly strips = new Strips();
  readonly warp = new PolarWarp();
  readonly labels = new RingLabels();
  cx = 0;
  cy = 0;
  /** The outer radius, and the rings' depth (their bands and the gaps between them). */
  r = 0;
  depth = 0;
  /** The unrolled circle's length (the outer circumference). */
  circumference = 1;
  /** The section's display units, and the title in the hole. */
  total = 0;
  title = '';
  subtitle = '';
  labelFont = '';
  /** The window: drawn, active (something follows), and its ends in the section's display units (end past start). */
  windowShown = false;
  windowActive = false;
  windowStart = 0;
  windowEnd = 0;
  /** A leader's points, reused. */
  readonly path = new Float64Array(6);

  warpOf(): PolarWarp {
    return this.warp;
  }

  /** The unit `unitFrame` reads and the frame x it writes: fields, not an argument and a return (boxed). */
  fu = 0;
  fx = 0;
  /**
   * Frame x of the section unit `fu` into `fx` (past the circle for a unit past the total): what the window's draw
   * calls, so nothing is passed (a number passed to a call that is not inlined is boxed).
   */
  private unitFrame(): void {
    const s = this.strips;
    this.fx = 0;
    if (!s.n) return;
    const t = this.total;
    let turns = 0;
    let v = this.fu;
    if (t > 0 && (v >= t || v < 0)) {
      turns = Math.floor(v / t);
      v -= turns * t;
    }
    let i = s.n - 1;
    while (i > 0 && s.first[i] > v) i--;
    const a = s.axes[i];
    // `xOfUnit` inline (a ring has no groups)
    this.fx = a.x0 + (Math.min(v - s.first[i], a.total) - a.u0) * a.scale + turns * this.circumference;
  }

  /** Frame x of a section unit (past the circle for a unit past the total). */
  frameX(u: number): number {
    this.fu = u;
    this.unitFrame();
    return this.fx;
  }

  /** The window's ends in the frame into `wx0`, `wx1` (`wx1` past the circle over the origin). */
  wx0 = 0;
  wx1 = 0;
  windowFrame(): void {
    this.fu = this.windowStart;
    this.unitFrame();
    this.wx0 = this.fx;
    this.fu = this.windowEnd;
    this.unitFrame();
    this.wx1 = this.fx <= this.wx0 ? this.fx + this.circumference : this.fx;
  }

  /** The item labels' font size, px. */
  labelSize = 10;
  /** The title's font size, px. */
  titleSize = 12;

  /**
   * Whether a page point is in placed label `k`'s hit box: its line (`labelLine` tall, centred where it was placed), so
   * two labels' boxes never overlap (their lines are a line apart at least).
   */
  inLabel(k: number, x: number, y: number): boolean {
    const l = this.labels;
    const x0 = l.right[k] ? l.x[k] : l.x[k] - l.width[k];
    const cy = l.y[k] - 0.35 * this.labelSize;
    const h = labelLine(this.labelSize) / 2;
    return x >= x0 - 2 && x < x0 + l.width[k] + 2 && y >= cy - h && y < cy + h;
  }

  /** The label under a page point, or −1. */
  labelAt(x: number, y: number): number {
    for (let k = 0; k < this.labels.m; k++) if (this.inLabel(k, x, y)) return k;
    return -1;
  }

  at(x: number, y: number, out: PlacePoint): boolean {
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height || x < b.x || x >= b.x + b.width) return false;
    if (this.labelAt(x, y) >= 0) {
      out.kind = 'label';
      return true;
    }
    const w = this.warp;
    if (!w.unmap(x, y)) return true;
    const fx = w.px;
    const fy = w.py;
    if (fy < -RingMetrics.rim || fy >= this.depth) return true;
    if (fy < 0) {
      if (this.windowActive) out.kind = 'window';
      return true;
    }
    const s = this.strips;
    const i = s.atX(fx);
    out.x = fx;
    out.y = fy;
    if (i < 0) {
      out.kind = 'gap';
      return true;
    }
    out.strip = i;
    for (let j = 0; j < s.bands; j++) {
      const k = i * s.bands + j;
      if (s.bandH[k] > 0 && fy >= s.bandY[k] && fy < s.bandY[k] + s.bandH[k]) {
        out.band = j;
        break;
      }
    }
    if (out.band < 0) {
      out.kind = 'gap';
      return true;
    }
    out.kind = 'data';
    out.offset = s.first[i] + s.axes[i].unitAt(fx);
    return true;
  }

  offsetAt(x: number, y: number): number {
    const s = this.strips;
    if (!s.n || !this.warp.unmap(x, y)) return -1;
    const fx = this.warp.px;
    // in a gap between blocks: the nearer block's edge
    const i = s.atX(fx);
    return s.offsetIn(i < 0 ? s.nearestX(fx) : i, fx);
  }

  place(offset: number, out: PlacePoint): boolean {
    if (!this.strips.placeAcross(offset, out)) return false;
    // through the warp: a strip with no band at the outer circle
    this.warp.map(out.x, out.band < 0 ? 0 : out.y);
    out.x = this.warp.px;
    out.y = this.warp.py;
    return true;
  }
}

/** The rings layout's arithmetic: its viewport is the window round the ring (wrapping over a cyclic origin). */
export const RingsNav: SectionNav = {
  pan(c, v, dx) {
    const d = (dx * (v.end - v.start)) / c.width;
    return clampWindow(c, { ...v, start: v.start - d, end: v.end - d });
  },
  zoom(c, v, factor, at) {
    // the unit under the pointer's angle, the turn nearest the window's middle (a window over the origin runs past
    // the total)
    const t = c.total;
    const a = c.cyclic && t > 0 ? at + Math.round(((v.start + v.end) / 2 - at) / t) * t : at;
    return TracksNav.zoom(c, v, factor, a);
  },
  fit: TracksNav.fit,
  focus(c, v, start, end) {
    // two ends on the circle: the shorter arc between them (an item over the origin, not the rest of the ring)
    const t = c.total;
    if (c.cyclic && t > 0 && end - start > t / 2 && end <= t)
      return clampWindow(c, { ...v, start: end, end: start + t });
    return clampWindow(c, { ...v, start, end });
  },
  reveal(c, v, start, end) {
    const t = c.total;
    const w = v.end - v.start;
    if (c.cyclic && t > 0) {
      // within the window, round the origin too
      const s = (((start - v.start) % t) + t) % t;
      if (s + (end - start) <= w) return v;
    } else if (start >= v.start && end <= v.end) return v;
    if (end - start > w) return clampWindow(c, { ...v, start, end });
    const m = (start + end) / 2 - w / 2;
    return clampWindow(c, { ...v, start: m, end: m + w });
  },
  step: TracksNav.step,
  clamp: clampWindow,
};

/** A ring's title: its system's name (a section of blocks: its id) and its length. */
function ringTitle(section: ResolvedSection, ctx: PlaceContext, total: number): [string, string] {
  const system = ctx.model?.systems.get(section.axis.system);
  const name = section.blocks ? sectionId(section) : (system?.name ?? sectionId(section));
  const unit = section.blocks ? '' : system?.alphabet === 'protein' ? ' aa' : ' bp';
  return [name, `${Math.round(total).toLocaleString('en-US')}${unit}`];
}

/** Rings outside in: the section's tracks as rings round its whole axis, a block an arc. */
export const RingsLayout: SectionLayoutType<RingPlacement> = {
  name: 'rings',
  fillBands: false,
  // its window passes the origin of a cyclic axis and draws there: the viewer's navigation wraps it
  wraps: true,
  place(section, ctx, _viewport, prev) {
    const p = prev ?? new RingPlacement();
    const { config, axis } = ctx;
    const rc = config.rings;
    const params = ringParams(section.params);
    const outside = params.labels === 'outside';
    const strips = p.strips;
    const nb = ctx.bands.length;
    strips.reset(nb);
    p.labels.reset();
    p.labelFont = config.itemFont;
    const size = fontSize(config.itemFont);
    p.labelSize = size;
    p.titleSize = fontSize(config.labelFont);
    // the diameter fits the box's width (its labels' room aside), capped
    const side = outside ? rc.labelGap + rc.leaderBend + RingMetrics.labelRoom : RingMetrics.margin + RingMetrics.rim;
    const written = section.params.height;
    const avail = Math.max(48, ctx.width - 2 * config.inset - 2 * side);
    const d = Math.min(rc.maxDiameter, typeof written === 'number' && written > 0 ? written : avail, avail);
    const r = d / 2;
    const margin = outside ? rc.labelGap + rc.leaderBend + size : RingMetrics.margin + RingMetrics.rim;
    const top = ctx.box.y;
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: d + 2 * margin };
    p.contentHeight = p.box.height;
    p.cx = config.inset + p.box.width / 2;
    p.cy = top + margin + r;
    p.r = r;
    // the rings' depth: their heights and gaps, scaled to keep the hole
    let depth = 0;
    for (let j = 0; j < nb; j++) depth += (j ? config.trackGap : 0) + ctx.heights[j];
    const room = r * (1 - Math.min(Math.max(rc.hole, 0), 0.95));
    const k = depth > room ? room / depth : 1;
    p.depth = depth * k;
    // the blocks: arcs of the circle less their gaps (a whole cyclic block alone has none)
    const blocks = section.blocks;
    const nBlocks = blocks?.length ?? 1;
    const whole =
      !blocks && axis.cyclic && axis.pieces.length === 2 && axis.pieces[0] === 0 && axis.total === axis.length;
    const gaps = whole ? 0 : nBlocks;
    const gap = (rc.blockGap * Math.PI) / 180;
    const arc = Math.max(TwoPi - gaps * gap, 1e-6);
    const circ = TwoPi * r;
    p.circumference = circ;
    const w = p.warp;
    w.cx = p.cx;
    w.cy = p.cy;
    w.r0 = r;
    w.x0 = 0;
    w.width = circ;
    w.a0 = params.rotation;
    w.sweep = TwoPi;
    w.tangent = rc.letters === 'tangent';
    w.step = 4;
    // each block's units over its share of the arc, from half a gap past the origin
    let total = 0;
    if (blocks) for (const b of blocks) total += b.axis.ranges ? piecesTotal(b.axis.ranges) : 0;
    else total = axis.total;
    p.total = total;
    let x = (gaps ? gap / 2 : 0) * r;
    let first = 0;
    for (let bi = 0; bi < nBlocks; bi++) {
      const b = blocks?.[bi];
      const pieces = b ? b.axis.ranges : axis.pieces;
      const len = pieces ? piecesTotal(pieces) : 0;
      if (!pieces || len <= 0) continue;
      const width = (len / Math.max(total, 1e-9)) * arc * r;
      const a = AxisTransform.on(b ? b.axis : axis, pieces, x, width, ctx.dpr);
      const i = strips.add(a, b ? bi : 0, 0, 0, first);
      strips.warp[i] = WarpKind.polar;
      strips.top[i] = p.box.y;
      strips.height[i] = p.box.height;
      let y = 0;
      for (let j = 0; j < nb; j++) {
        if (j) y += config.trackGap * k;
        strips.bandY[i * nb + j] = y;
        strips.bandH[i * nb + j] = ctx.heights[j] * k;
        y += ctx.heights[j] * k;
      }
      x += width + gap * r;
      first += len;
    }
    // the strips' layers, their labels into the sink (placed below, once the ring is laid out)
    for (let i = 0; i < strips.n; i++) {
      Sink.at(p, i);
      ctx.layStrip(p, i, outside ? Sink : undefined);
      // a block's own label leads its labels
      const b = blocks?.[strips.block[i]];
      if (b && outside) {
        const a = strips.axes[i];
        p.labels.add(b.label, angleOf(p, a.x0 + a.width / 2), Infinity, -1, undefined);
      }
    }
    if (outside) placeLabels(p, config, size, ctx.text);
    const [title, sub] = rc.center === 'title' ? ringTitle(section, ctx, total) : ['', ''];
    p.title = title;
    p.subtitle = sub;
    return p;
  },
  // rotation and size are the section's and the view's: the window is never part of the ring's key
  sameKey: () => true,
  nav: RingsNav,
  followWindow(_section, axis, viewport) {
    return viewport ? { start: viewport.start, end: viewport.end } : { start: 0, end: axis.total };
  },
  window(placement, viewport, active, shown) {
    // a pending section's placement is not a ring's
    if (!(placement instanceof RingPlacement)) return;
    const p = placement;
    p.windowShown = shown;
    p.windowStart = viewport ? viewport.start : 0;
    p.windowEnd = viewport ? viewport.end : p.total;
    // a window with nowhere to go (the whole ring) is idle: quiet, no wash, grips or hover
    p.windowActive = active && p.windowEnd - p.windowStart < p.total;
  },
  chrome(s, pass, out, emphasis, strip, band) {
    // pending (not a ring's placement), or no strips (an empty range): a track's chrome is called for strip 0
    const p = s.placement;
    if (!(p instanceof RingPlacement)) return;
    const { tokens, config } = s;
    if (band >= 0) {
      // each ring's own ground under its layers, a sector a strip (through `args`: nothing boxed)
      if (pass !== 'data' || strip >= p.strips.n) return;
      const a = p.strips.axes[strip];
      const k = strip * p.strips.bands + band;
      const g = out.args;
      g[0] = a.x0;
      g[1] = p.strips.bandY[k];
      g[2] = a.width;
      g[3] = p.strips.bandH[k];
      g[4] = tokens.colors['seq-track'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
      return;
    }
    if (pass === 'data') p.warp.refused = 0;
    else if (pass === 'window') drawRingWindow(p, s, out, emphasis, tokens);
    else if (pass === 'labels') drawRingLabels(p, out, tokens, config);
  },
  hitChrome(s, x, y) {
    const p = s.placement;
    if (!(p instanceof RingPlacement)) return undefined;
    const k = p.labelAt(x, y);
    if (k >= 0) return p.labels.hits[p.labels.placed[k]];
    if (!p.strips.n || !p.windowActive || !p.windowShown || !p.warp.unmap(x, y)) return undefined;
    const fx = p.warp.px;
    const fy = p.warp.py;
    p.windowFrame();
    const c = p.circumference;
    // the pointer's frame x a turn either way too: a grip reaches across the origin from either side
    const part =
      windowPart(fx, fy, p.wx0, p.wx1, -RingMetrics.rim, RingMetrics.rim) ??
      windowPart(fx + c, fy, p.wx0, p.wx1, -RingMetrics.rim, RingMetrics.rim) ??
      windowPart(fx - c, fy, p.wx0, p.wx1, -RingMetrics.rim, RingMetrics.rim);
    if (!part) return undefined;
    return { layer: s.section.address, section: s.id, system: s.section.axis.system, target: 'window', part };
  },
};

/** A ring label's line: what two labels on a side are spaced by, and each one's hit box's height. */
function labelLine(size: number): number {
  return size + 2;
}

/** The angle of frame x, in [0, 2π), clockwise from the top. */
function angleOf(p: RingPlacement, x: number): number {
  const t = p.warp.angleAt(x);
  return ((t % TwoPi) + TwoPi) % TwoPi;
}

/** The label sink a strip's layers give their labels to: each at its item's middle angle. */
class RingSink implements LabelSink {
  private p!: RingPlacement;
  private strip = 0;
  at(p: RingPlacement, strip: number) {
    this.p = p;
    this.strip = strip;
  }
  label(text: string, start: number, end: number, _element: ElementKey, hit?: Hit, order = 0): void {
    const p = this.p;
    const a = p.strips.axes[this.strip];
    const xs = a.x(start);
    let xe = a.x(end);
    if (!Number.isFinite(xs) || !Number.isFinite(xe)) return;
    // an item over the origin runs past the circle
    if (xe < xs) xe += p.circumference;
    const len = end > start ? end - start : end + (a.length ?? a.total) - start;
    p.labels.add(text, angleOf(p, (xs + xe) / 2), len, order, hit);
  }
}
const Sink = new RingSink();

/**
 * Places a ring's labels: each side (right: angles in [0, π)) sorted by height, spread a line apart within the box
 * (a forward pass then a backward one), past what a side holds dropped by order then length; each label hugging the
 * circle at its leader's outer radius, its leader out along its angle to the bend, then to the label.
 */
function placeLabels(p: RingPlacement, config: RenderConfig, size: number, text: TextMetrics) {
  const l = p.labels;
  const rc = config.rings;
  const line = labelLine(size);
  const r1 = p.r + 2;
  const rb = p.r + rc.labelGap;
  const rl = p.r + rc.labelGap + rc.leaderBend;
  const lo = p.box.y + size;
  const hi = p.box.y + p.box.height - 2;
  const hold = Math.max(0, Math.floor((hi - lo) / line) + 1);
  l.ensurePlaced(l.n);
  for (let s = 0; s < 2; s++) {
    const right = s === 0;
    const idx: number[] = [];
    for (let i = 0; i < l.n; i++) if (l.angle[i] < Math.PI === right) idx.push(i);
    // past what the side holds: kept by order (the block labels first, then the layers' order), then the longer item
    if (idx.length > hold) {
      idx.sort((a, b) => l.order[a] - l.order[b] || l.length[b] - l.length[a] || a - b);
      idx.length = hold;
    }
    const want = (i: number) => p.cy - rb * Math.cos(l.angle[i]);
    idx.sort((a, b) => want(a) - want(b) || a - b);
    const ys = idx.map(want);
    for (let k = 0; k < ys.length; k++) ys[k] = Math.max(ys[k], k ? ys[k - 1] + line : lo);
    for (let k = ys.length - 1; k >= 0; k--) ys[k] = Math.min(ys[k], k + 1 < ys.length ? ys[k + 1] - line : hi);
    for (let k = 0; k < idx.length; k++) {
      const i = idx[k];
      const m = l.m++;
      const t = l.angle[i];
      const y = ys[k];
      const dy = y - p.cy;
      const dx = Math.sqrt(Math.max(0, rl * rl - dy * dy));
      const x = right ? p.cx + Math.max(dx, rc.leaderBend) : p.cx - Math.max(dx, rc.leaderBend);
      l.placed[m] = i;
      l.x[m] = x + (right ? 3 : -3);
      l.y[m] = y + 0.35 * size;
      l.right[m] = right ? 1 : 0;
      l.width[m] = Math.min(RingMetrics.labelRoom, text.measure(config.itemFont, l.text[i]));
      const o = 6 * m;
      l.leader[o] = p.cx + r1 * Math.sin(t);
      l.leader[o + 1] = p.cy - r1 * Math.cos(t);
      l.leader[o + 2] = p.cx + rb * Math.sin(t);
      l.leader[o + 3] = p.cy - rb * Math.cos(t);
      l.leader[o + 4] = x;
      l.leader[o + 5] = y;
    }
  }
}

/** The ring's labels and leaders (on the page), and its title in the hole. */
function drawRingLabels(p: RingPlacement, out: PassWriter, tokens: SeqTokens, config: RenderConfig) {
  const l = p.labels;
  const g = out.args;
  const muted = tokens.colors['msui-text-muted'];
  const ink = tokens.colors['msui-text'];
  const xy = p.path;
  for (let k = 0; k < l.m; k++) {
    const o = 6 * k;
    for (let q = 0; q < 6; q++) xy[q] = l.leader[o + q];
    g[4] = NoColor;
    g[5] = muted;
    g[6] = 1;
    g[7] = 0;
    out.pathArgs(xy, 3, false);
    g[0] = l.x[k];
    g[1] = l.y[k];
    g[2] = l.width[k];
    g[3] = RingMetrics.labelRoom;
    g[4] = ink;
    out.textArgs(l.text[l.placed[k]], l.right[k] ? AnchorStart : AnchorEnd, config.itemFont);
  }
  if (!p.title) return;
  const size = p.titleSize;
  const hole = Math.max(0, p.r - p.depth);
  g[0] = p.cx;
  g[1] = p.cy - 2;
  g[2] = 0;
  g[3] = 2 * hole * 0.9;
  g[4] = ink;
  out.textArgs(p.title, AnchorMiddle, config.labelFont);
  g[0] = p.cx;
  g[1] = p.cy + size + 2;
  g[2] = 0;
  g[3] = 2 * hole * 0.9;
  g[4] = muted;
  out.textArgs(p.subtitle, AnchorMiddle, config.labelFont);
}

/**
 * The window (`render/window.ts`) through the warp: a faint sector over the rings, and on the rim its body and its
 * grips, lit by the hovered part; idle (nothing follows) a quiet outline.
 */
function drawRingWindow(
  p: RingPlacement,
  s: SectionLayout,
  out: PassWriter,
  emphasis: Emphasis | undefined,
  tokens: SeqTokens,
) {
  if (!p.windowShown || !p.strips.n) return;
  p.windowFrame();
  const x0 = p.wx0;
  const x1 = p.wx1;
  const h = emphasis?.hover;
  const lit = h?.target === 'window' && h.section === s.id ? h.part : undefined;
  const w = p.warp.writer(out);
  if (p.windowActive) {
    const g = w.args;
    const ink = tokens.colors['seq-window'];
    // a narrow window grows about its middle (`windowExtent`'s rule, the rim's)
    const wide = x1 - x0 >= WindowMetrics.floor;
    g[0] = wide ? x0 : (x0 + x1) / 2 - WindowMetrics.floor / 2;
    g[1] = 0;
    g[2] = wide ? x1 - x0 : WindowMetrics.floor;
    g[3] = p.depth;
    g[4] = ((ink & 0xffffff00) | Math.round(WindowMetrics.fill * 0.5 * 255)) >>> 0;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    w.rectArgs();
  }
  const g = w.args;
  g[0] = x0;
  g[1] = x1;
  g[2] = -RingMetrics.rim;
  g[3] = RingMetrics.rim;
  drawWindowArgs(w, p.windowActive, lit, tokens);
  p.warp.close(out);
}

export type { NavContext, Viewport };
