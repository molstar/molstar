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
import type { Key } from '@molstar/seqx-core/query';
import { Anchor, grown, NoColor, type PassWriter } from '../draw-lists.js';
import type { Box, Emphasis, Hit, LabelSink } from '../layer-type.js';
import { sectionId } from '../layer-type.js';
import type { TextMetrics } from '../text.js';
import { type SeqTokens, withAlpha } from '../tokens.js';
import type { SectionLayout } from '../tracks.js';
import { drawWindowArgs, floorExtent, WindowFill, windowPart } from '../window.js';
import {
  blockLabel,
  blockLength,
  blockNamePrefix,
  chipSlots,
  isChip,
  type ResolvedBlock,
  runLabel,
  shortBlockLabel,
} from './blocks.js';
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
export const RingMetrics = { labelRoom: 96, rim: 12, margin: 4, groupBand: 16 } as const;

/** A rings section's params as the spec writes them (_spec/README.md "Resolved view"), with their defaults. */
export interface RingParams {
  /** Where the origin is, radians clockwise from the top. */
  readonly rotation: number;
  readonly labels: 'outside' | 'none';
  /** What the hole holds (`section.center`): absent, the render config's `rings.center`. */
  readonly center: 'title' | 'none' | undefined;
}

export function ringParams(params: ResolvedSection['params']): RingParams {
  const r = params.rotation;
  return {
    rotation: typeof r === 'number' && Number.isFinite(r) ? (r * Math.PI) / 180 : 0,
    labels: params.labels === 'none' ? 'none' : 'outside',
    center: params.center === 'title' || params.center === 'none' ? params.center : undefined,
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
  /**
   * Dense groups (`findGroups`): each given label's group (−1: none), the groups' sizes, middle angles, best orders and
   * summed lengths, their members flat from `groupStart`; scratch (the sort, a cursor) kept between layouts.
   */
  ng = 0;
  group = new Int32Array(16);
  sorted = new Int32Array(16);
  cursor = new Int32Array(16);
  members = new Int32Array(16);
  groupStart = new Int32Array(16);
  groupSize = new Int32Array(16);
  groupAngle = new Float64Array(16);
  groupOrder = new Float64Array(16);
  groupLength = new Float64Array(16);
  /** The full set's placed slot of each given label (−1: dropped), for the fan. */
  slot = new Int32Array(16);
  /** The grouped set's: each given label's group (−1: its own item), and each group's hit. */
  of = new Int32Array(16);
  groupHits: Hit[] = [];

  reset() {
    this.n = 0;
    this.m = 0;
    this.ng = 0;
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
  /** Each block's chip slot (px of the outer ring; 0 for a block that is no chip). */
  slots: Float64Array = new Float64Array(4);
  /** The labels as drawn at rest: a dense fan is one grouped label (`layoutLabels`). */
  labels = new RingLabels();
  /** Every label placed as if no group existed, when a fan is dense (m is 0 otherwise): what a hover on a group opens. */
  fan = new RingLabels();
  /** The most members a group's fan opens to (`rings.cluster.fanMax`): a larger group stays one label, its card the list. */
  fanMax = Infinity;
  /** The fan's labels' own angles while their pivoted ones are placed (scratch, kept between layouts). */
  fanFrom = new Float64Array(16);
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
  /** The title's and subtitle's widths, CSS px: the plate under them (chords cross the hole). */
  titleW = 0;
  labelFont = '';
  /** The window: drawn, active (something follows), and its ends in the section's display units (end past start). */
  windowShown = false;
  windowActive = false;
  windowStart = 0;
  windowEnd = 0;
  /** A leader's points, reused. */
  readonly path = new Float64Array(6);
  /** The gap between two blocks as laid, radians (the share rule: at most `blockGap`). */
  gapAngle = 0;
  /** Each block's frame x and width as laid (a block in a run too), by its index. */
  bx = new Float64Array(4);
  bw = new Float64Array(4);
  /** Each block's run (its first block's index) or −1, at the outer radius's zoom. */
  run = new Int32Array(4);
  /** The runs: their first and last blocks and frame extent. */
  runN = 0;
  runFirst = new Int32Array(4);
  runLast = new Int32Array(4);
  runX = new Float64Array(4);
  runW = new Float64Array(4);
  /** The groups' arcs (a section of several assemblies): frame extents, drawn outside the outer ring. */
  groupN = 0;
  groupX = new Float64Array(4);
  groupW = new Float64Array(4);
  /** The groups' names, drawn in their band outside the outer ring (never among the block labels). */
  groupText: string[] = [];
  /**
   * Each group's name's centre on the page (`groupLayout`): on its arc's middle angle, past the arcs and the block
   * labels there (never over an arc or a label).
   */
  groupTX = new Float64Array(4);
  groupTY = new Float64Array(4);
  /** The band outside the outer ring the groups' arcs and names take (0 without groups): the labels start past it. */
  groupBand = 0;

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
    // `xOfUnit` inline (a ring has no groups); a block's units its positions times its scale
    this.fx = a.x0 + (Math.min((v - s.first[i]) / s.scale[i], a.total) - a.u0) * a.scale + turns * this.circumference;
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
    return this.inLabelOf(this.labels, k, x, y);
  }

  /** `inLabel` of a label set's placed label `k` (the grouped set's or the fan's). */
  inLabelOf(l: RingLabels, k: number, x: number, y: number): boolean {
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

  /**
   * The open group's member whose fan label is under a page point (its index among the fan's given labels), or −1: a
   * fan is hit only while its group is open (`openGroupOf` the hover: the draw and the hit read the same input).
   */
  fanAt(x: number, y: number, hover?: Hit): number {
    const g = openGroupOf(this, hover);
    if (g < 0) return -1;
    const f = this.fan;
    for (let q = f.groupStart[g]; q < f.groupStart[g + 1]; q++) {
      const i = f.members[q];
      const k = f.slot[i];
      if (k >= 0 && this.inLabelOf(f, k, x, y)) return i;
    }
    return -1;
  }

  at(p: ArrayLike<number>, out: PlacePoint): boolean {
    const x = p[0];
    const y = p[1];
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height || x < b.x || x >= b.x + b.width) return false;
    if (this.fanAt(x, y) >= 0 || this.labelAt(x, y) >= 0) {
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
    out.offset = s.first[i] + s.axes[i].unitAt(fx) * s.scale[i];
    return true;
  }

  offsetAt(x: number, y: number): number {
    const s = this.strips;
    if (!s.n || !this.warp.unmap(x, y)) return -1;
    const fx = this.warp.px;
    // in a gap between blocks: the nearer block's edge
    const i = s.atX(fx);
    const k = i < 0 ? s.nearestX(fx) : i;
    // the exact centre has no angle: no unit
    return k < 0 ? -1 : s.offsetIn(k, fx);
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

/** A length in bases as a genome's title says it: `3.1 Gb`, `48 Mb`, `16.6 kb`, `950 bp`. */
export function basesText(n: number): string {
  const units: [number, string][] = [
    [1e9, 'Gb'],
    [1e6, 'Mb'],
    [1e3, 'kb'],
  ];
  for (const [k, u] of units) if (n >= k) return `${+(n / k).toPrecision(n / k >= 100 ? 3 : 2)} ${u}`;
  return `${Math.round(n)} bp`;
}

/**
 * A ring's title: its system's name (a section of blocks: its id) and its length; a section of one assembly its name
 * and its blocks' total length in bases (`GRCh38`, `3.1 Gb`), of several their names; the closing run's count after.
 */
function ringTitle(section: ResolvedSection, ctx: PlaceContext, total: number): [string, string] {
  const system = ctx.model?.systems.get(section.axis.system);
  const groups = section.groups;
  if (groups?.length && section.blocks) {
    let bases = 0;
    for (const b of section.blocks) bases += piecesTotal(b.axis.ranges ?? []) * b.members.length;
    const more = section.more ? ` (+${section.more.toLocaleString('en-US')} not shown)` : '';
    return [groups.map((g) => g.label).join(' · '), `${basesText(bases)}${more}`];
  }
  const name = section.blocks ? sectionId(section) : (system?.name ?? sectionId(section));
  const unit = section.blocks ? '' : system?.alphabet === 'protein' ? ' aa' : ' bp';
  return [name, `${Math.round(total).toLocaleString('en-US')}${unit}`];
}

/**
 * A ring's runs at the outer radius's zoom (px a unit): two or more consecutive blocks of one group, no chips, each
 * under `blocks.runPx` (_spec/README.md "Render", Rings layout); each block's run start or −1, undefined with none.
 */
function ringRuns(
  p: RingPlacement,
  blocks: readonly ResolvedBlock[] | undefined,
  config: RenderConfig,
  zoom: number,
): Int32Array | undefined {
  const small = config.blocks.runPx;
  if (!blocks || !(small > 0)) return undefined;
  const n = blocks.length;
  p.run = grown(p.run, n);
  const tiny = (b: number) => !isChip(blocks[b]) && blockLength(blocks[b]) * zoom < small;
  let any = false;
  let b = 0;
  while (b < n) {
    if (!tiny(b)) {
      p.run[b++] = -1;
      continue;
    }
    let e = b + 1;
    while (e < n && tiny(e) && (blocks[e].group ?? -1) === (blocks[b].group ?? -1)) e++;
    const v = e - b > 1 ? b : -1;
    if (v >= 0) any = true;
    for (let k = b; k < e; k++) p.run[k] = v;
    b = e;
  }
  return any ? p.run : undefined;
}

/** A block of a run laid at frame x: its run's extent grown (a run's first block opens it). */
function addRun(p: RingPlacement, start: number, b: number, x: number, w: number): void {
  if (start === b) {
    const r = p.runN++;
    if (r >= p.runFirst.length) {
      p.runFirst = grown(p.runFirst, r + 1);
      p.runLast = grown(p.runLast, r + 1);
      p.runX = grown(p.runX, r + 1);
      p.runW = grown(p.runW, r + 1);
    }
    p.runFirst[r] = b;
    p.runX[r] = x;
    p.runW[r] = 0;
  }
  const r = p.runN - 1;
  p.runLast[r] = b;
  p.runW[r] = x + w - p.runX[r];
}

/** The groups' arcs of a section of several assemblies (EC-210): each its blocks' frame extent. */
function layGroupArcs(
  p: RingPlacement,
  section: ResolvedSection,
  blocks: readonly ResolvedBlock[] | undefined,
  _end: number,
): void {
  p.groupN = 0;
  const groups = section.groups;
  if (!groups || groups.length < 2 || !blocks) return;
  const s = p.strips;
  for (const g of groups) {
    let x0 = Infinity;
    let x1 = -Infinity;
    for (let i = 0; i < s.n; i++) {
      const b = s.block[i];
      if (b < g.first || b >= g.first + g.count) continue;
      x0 = Math.min(x0, s.axes[i].x0);
      x1 = Math.max(x1, s.axes[i].x0 + s.axes[i].width);
    }
    for (let r = 0; r < p.runN; r++)
      if (p.runFirst[r] >= g.first && p.runFirst[r] < g.first + g.count) {
        x0 = Math.min(x0, p.runX[r]);
        x1 = Math.max(x1, p.runX[r] + p.runW[r]);
      }
    if (!(x1 > x0)) continue;
    const k = p.groupN++;
    if (k >= p.groupX.length) {
      p.groupX = grown(p.groupX, k + 1);
      p.groupW = grown(p.groupW, k + 1);
    }
    p.groupX[k] = x0;
    p.groupW[k] = x1 - x0;
    // its name outside, at its middle angle: placed with the blocks' labels, before them
    p.groupText[k] = g.label;
  }
}

/**
 * Each group's name (EC-210): at its arc's middle angle, out along it from the arcs to the first place where its box
 * meets no placed block label, no other group's name and no arc (its nearest point past the arcs), inside the section's
 * box; with none there, the farthest place inside the box along the angle. Laid out with the labels, never on a hover.
 */
function placeGroupNames(p: RingPlacement, text: TextMetrics, font: string): void {
  if (p.groupTX.length < p.groupN) {
    p.groupTX = grown(p.groupTX, p.groupN);
    p.groupTY = grown(p.groupTY, p.groupN);
  }
  const l = p.labels;
  const h = p.titleSize + 2;
  const lh = labelLine(p.labelSize) / 2;
  const box = p.box;
  const clear = p.r + 6;
  for (let k = 0; k < p.groupN; k++) {
    const w = text.measure(font, p.groupText[k]) + 4;
    const t = angleOf(p, p.groupX[k] + p.groupW[k] / 2);
    const sx = Math.sin(t);
    const sy = -Math.cos(t);
    let best = -1;
    let last = -1;
    for (let d = 0; d < 4 * (box.width + box.height); d += 2) {
      const cx = p.cx + d * sx;
      const cy = p.cy + d * sy;
      const x0 = cx - w / 2;
      const y0 = cy - h / 2;
      // inside the section's box
      if (x0 < box.x || x0 + w > box.x + box.width || y0 < box.y || y0 + h > box.y + box.height) {
        if (last >= 0) break;
        continue;
      }
      // its nearest point past the arcs
      const nx = Math.max(x0, Math.min(p.cx, x0 + w)) - p.cx;
      const ny = Math.max(y0, Math.min(p.cy, y0 + h)) - p.cy;
      if (nx * nx + ny * ny < clear * clear) continue;
      last = d;
      let hit = false;
      for (let q = 0; q < l.m && !hit; q++) {
        const lx = l.right[q] ? l.x[q] : l.x[q] - l.width[q];
        const ly = l.y[q] - 0.35 * p.labelSize;
        hit = lx < x0 + w && lx + l.width[q] > x0 && ly - lh < y0 + h && ly + lh > y0;
      }
      for (let q = 0; q < k && !hit; q++) hit = Math.abs(p.groupTX[q] - cx) < w && Math.abs(p.groupTY[q] - cy) < h;
      if (!hit) {
        best = d;
        break;
      }
    }
    const d = best >= 0 ? best : Math.max(0, last);
    p.groupTX[k] = p.cx + d * sx;
    p.groupTY[k] = p.cy + d * sy;
  }
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
    p.groupBand = section.groups && section.groups.length > 1 && outside ? RingMetrics.groupBand : 0;
    // the groups' names past the block labels: a band of their own outside them, as much again
    const names = p.groupBand ? RingMetrics.groupBand : 0;
    const size = fontSize(config.itemFont);
    p.labelSize = size;
    p.titleSize = fontSize(config.labelFont);
    // the diameter fits the box's width (its labels' room aside), capped
    const side = outside
      ? rc.labelGap + p.groupBand + rc.leaderBend + RingMetrics.labelRoom + names
      : RingMetrics.margin + RingMetrics.rim;
    const written = section.params.height;
    const avail = Math.max(48, ctx.width - 2 * config.inset - 2 * side);
    const d = Math.min(rc.maxDiameter, typeof written === 'number' && written > 0 ? written : avail, avail);
    const r = d / 2;
    const margin = outside
      ? rc.labelGap + p.groupBand + rc.leaderBend + size + names
      : RingMetrics.margin + RingMetrics.rim;
    const top = ctx.box.y;
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: d + 2 * margin };
    p.contentHeight = p.box.height;
    p.cx = config.inset + p.box.width / 2;
    p.cy = top + margin + r;
    p.r = r;
    // the rings' depth: their heights and gaps, scaled to keep the hole
    // (a ring's padding is radial: its well is its height plus its padding, scaled with it)
    const pads = ctx.pads;
    let depth = 0;
    for (let j = 0; j < nb; j++)
      depth += (j ? config.trackGap : 0) + ctx.heights[j] + (pads ? pads[2 * j] + pads[2 * j + 1] : 0);
    const room = r * (1 - Math.min(Math.max(rc.hole, 0), 0.95));
    const k = depth > room ? room / depth : 1;
    p.depth = depth * k;
    // the blocks: arcs of the circle less their gaps (a whole cyclic block alone has none)
    const blocks = section.blocks;
    const nBlocks = blocks?.length ?? 1;
    const whole =
      !blocks && axis.cyclic && axis.pieces.length === 2 && axis.pieces[0] === 0 && axis.total === axis.length;
    const circ = TwoPi * r;
    // each block's units over its share of the arc
    let total = 0;
    if (blocks) for (const b of blocks) total += blockLength(b);
    else total = axis.total;
    p.total = total;
    // the gaps by one rule with rows (_spec/README.md "Render", Rings layout): a gap after each block (none inside a run
    // of small blocks at the outer radius, none for a whole cyclic block), a group's boundary `groupGap` degrees taken
    // first, the rest sharing at most `gapShare` of the circle, each at most `blockGap`
    if (blocks) {
      p.bx = grown(p.bx, nBlocks);
      p.bw = grown(p.bw, nBlocks);
    }
    const runs = ringRuns(p, blocks, config, total > 0 ? circ / total : 0);
    let inner = 0;
    let groupBounds = 0;
    if (!whole)
      for (let b = 0; b < nBlocks; b++) {
        const next = (b + 1) % nBlocks;
        if (runs && runs[b] >= 0 && next > 0 && runs[next] === runs[b]) continue;
        if (blocks && nBlocks > 1 && (blocks[b].group ?? -1) !== (blocks[next].group ?? -1)) groupBounds++;
        else inner++;
      }
    const groupGap = (rc.groupGap * Math.PI) / 180;
    const fixed = Math.min(groupBounds * groupGap, TwoPi * rc.gapShare);
    const gapOf = groupBounds ? fixed / groupBounds : 0;
    const gap = Math.min((rc.blockGap * Math.PI) / 180, Math.max(0, TwoPi * rc.gapShare - fixed) / Math.max(1, inner));
    const arc = Math.max(TwoPi - inner * gap - fixed, 1e-6);
    p.gapAngle = gap;
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
    // a chip (a one-position block) a fixed angle, its label's width on the outer ring (at most half the circle
    // together): the polymers share the rest
    let chipPx = 0;
    // the units the chips stand for (a `whole` chip of several, each of them)
    let chips = 0;
    const slots = blocks ? (p.slots = chipSlots(ctx.text, config, blocks, arc * r, p.slots)) : p.slots;
    if (blocks)
      for (let b = 0; b < blocks.length; b++)
        if (slots[b] > 0) {
          chipPx += slots[b];
          chips += blockLength(blocks[b]);
        }
    const shared = Math.max(arc * r - chipPx, 1e-6);
    // from half the first gap past the origin (the gap before the first block, the last one's)
    const before = whole
      ? 0
      : nBlocks > 1 && blocks && (blocks[nBlocks - 1].group ?? -1) !== (blocks[0].group ?? -1)
        ? gapOf
        : gap;
    let x = (before / 2) * r;
    let first = 0;
    p.runN = 0;
    for (let bi = 0; bi < nBlocks; bi++) {
      const b = blocks?.[bi];
      const pieces = b ? b.axis.ranges : axis.pieces;
      const len = !pieces ? 0 : b ? blockLength(b) : piecesTotal(pieces);
      if (!pieces || len <= 0) continue;
      const chip = b !== undefined && slots[bi] > 0;
      const width = chip
        ? slots[bi]
        : chips
          ? (len / Math.max(total - chips, 1e-9)) * shared
          : (len / Math.max(total, 1e-9)) * arc * r;
      const next = (bi + 1) % nBlocks;
      const inRun = runs !== undefined && runs[bi] >= 0;
      const after = whole
        ? 0
        : inRun && next > 0 && runs[next] === runs[bi]
          ? 0
          : blocks && nBlocks > 1 && (blocks[bi].group ?? -1) !== (blocks[next].group ?? -1)
            ? gapOf
            : gap;
      if (blocks) {
        p.bx[bi] = x;
        p.bw[bi] = width;
      }
      if (inRun) addRun(p, runs[bi], bi, x, width);
      // inside a run a block under `readPx` device px is not read nor drawn: the run's band alone (the joint rule)
      if (inRun && width * ctx.dpr < config.blocks.readPx) {
        x += width + after * r;
        first += len;
        continue;
      }
      const a = AxisTransform.on(b ? b.axis : axis, pieces, x, width, ctx.dpr);
      const i = strips.add(a, b ? bi : 0, 0, 0, first);
      if (b && !chip) strips.scale[i] = b.scale ?? 1;
      if (chip) strips.chip[i] = 1;
      strips.warp[i] = WarpKind.polar;
      strips.top[i] = p.box.y;
      strips.height[i] = p.box.height;
      let y = 0;
      for (let j = 0; j < nb; j++) {
        if (j) y += config.trackGap * k;
        const pt = (pads ? pads[2 * j] : 0) * k;
        const pb = (pads ? pads[2 * j + 1] : 0) * k;
        strips.padT[i * nb + j] = pt;
        strips.padB[i * nb + j] = pb;
        strips.bandY[i * nb + j] = y + pt;
        strips.bandH[i * nb + j] = ctx.heights[j] * k;
        y += pt + ctx.heights[j] * k + pb;
      }
      x += width + after * r;
      first += len;
    }
    // the strips' layers, their labels into the sink (placed below, once the ring is laid out)
    const prefix = blocks ? blockNamePrefix(blocks) : '';
    for (let i = 0; i < strips.n; i++) {
      Sink.at(p, i);
      const pill = ctx.layStrip(p, i, outside ? Sink : undefined);
      // a block's own label leads its labels, unless its chip's pill draws it
      const bi = strips.block[i];
      const b = blocks?.[bi];
      if (b && outside && !pill && !(runs && runs[bi] >= 0)) {
        const a = strips.axes[i];
        // whole where its arc holds it, else short (its name less the names' shared prefix: `11`, `X`), as a row's
        const whole = blockLabel(b, config);
        const short =
          whole && ctx.text.measure(config.labelFont, whole) > a.width ? shortBlockLabel(b, config, prefix) : '';
        p.labels.add(short || whole, angleOf(p, a.x0 + a.width / 2), Infinity, -1, undefined);
      }
    }
    // a run's label where it has room, then the groups' names outside the labels (EC-210)
    if (outside)
      for (let k = 0; k < p.runN; k++) {
        const n = p.runLast[k] - p.runFirst[k] + 1;
        p.labels.add(runLabel(n), angleOf(p, p.runX[k] + p.runW[k] / 2), Infinity, -1, undefined);
      }
    if (outside) layoutLabels(p, config, size, ctx.text);
    layGroupArcs(p, section, blocks, x);
    if (outside) placeGroupNames(p, ctx.text, config.labelFont);
    const [title, sub] = (params.center ?? rc.center) === 'title' ? ringTitle(section, ctx, total) : ['', ''];
    p.title = title;
    p.subtitle = sub;
    p.titleW = title ? Math.max(ctx.text.measure(config.labelFont, title), ctx.text.measure(config.labelFont, sub)) : 0;
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
      g[1] = p.strips.bandY[k] - p.strips.padT[k];
      g[2] = a.width;
      g[3] = p.strips.bandH[k] + p.strips.padT[k] + p.strips.padB[k];
      g[4] = tokens.colors['seq-track'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
      return;
    }
    if (pass === 'data') {
      p.warp.refused = 0;
      p.warp.reduced = 0;
      p.warp.runsUsed = 0;
      p.warp.imageRuns = config.rings.imageRuns;
      if (p.runN || p.groupN) drawRingRuns(p, out, tokens);
    } else if (pass === 'window') drawRingWindow(p, s, out, emphasis, tokens);
    else if (pass === 'labels') drawRingLabels(p, out, tokens, config);
    else if (pass === 'tip') drawRingFan(p, out, emphasis, tokens, config);
  },
  hitChrome(s, x, y, hover) {
    const p = s.placement;
    if (!(p instanceof RingPlacement)) return undefined;
    const i = p.fanAt(x, y, hover);
    if (i >= 0) return p.fan.hits[i];
    const k = p.labelAt(x, y);
    if (k >= 0) return p.labels.hits[p.labels.placed[k]];
    if (!p.strips.n || !p.windowActive || !p.windowShown || !p.warp.unmap(x, y)) return undefined;
    const fx = p.warp.px;
    const fy = p.warp.py;
    p.windowFrame();
    const c = p.circumference;
    // the pointer's frame x a turn either way too: a grip reaches across the origin from either side
    const wc = s.config.navigator.window;
    const part =
      windowPart(fx, fy, p.wx0, p.wx1, -RingMetrics.rim, RingMetrics.rim, false, wc) ??
      windowPart(fx + c, fy, p.wx0, p.wx1, -RingMetrics.rim, RingMetrics.rim, false, wc) ??
      windowPart(fx - c, fy, p.wx0, p.wx1, -RingMetrics.rim, RingMetrics.rim, false, wc);
    if (!part) return undefined;
    return { layer: s.section.address, section: s.id, system: s.section.axis.system, target: 'window', part };
  },
};

/**
 * A ring's runs (a band in `seq-wash` over each run's arc through the rings, a hairline at each end) and its groups'
 * arcs (a hairline just outside the outer ring), through the warp in the data pass (_spec/README.md "Render", Rings
 * layout); through `args`, nothing boxed.
 */
function drawRingRuns(p: RingPlacement, page: PassWriter, tokens: SeqTokens): void {
  const out = p.warp.writer(page);
  const g = out.args;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  for (let r = 0; r < p.runN; r++) {
    g[0] = p.runX[r];
    g[1] = 0;
    g[2] = p.runW[r];
    g[3] = p.depth;
    g[4] = tokens.colors['seq-wash'];
    out.rectArgs();
    for (let e = 0; e < 2; e++) {
      g[0] = e ? p.runX[r] + p.runW[r] - 0.75 : p.runX[r];
      g[1] = 0;
      g[2] = 0.75;
      g[3] = p.depth;
      g[4] = tokens.colors['msui-text-muted'];
      out.rectArgs();
    }
  }
  for (let k = 0; k < p.groupN; k++) {
    g[0] = p.groupX[k];
    g[1] = -3;
    g[2] = p.groupW[k];
    g[3] = 1.5;
    g[4] = tokens.colors['msui-text-muted'];
    out.rectArgs();
  }
  p.warp.close(page);
}

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
  label(text: string, start: number, end: number, _element: Key, hit?: Hit, order = 0): void {
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

/** The "n items" texts, one string a count for the process (a layout makes none). */
const ItemsTexts: string[] = [];
const itemsText = (n: number) => (ItemsTexts[n] ??= `${n} items`);

/** The sort's keys and the comparator reading them: a module-level function, no closure made a layout. */
let SortKeys = new Float64Array(64);
const byKey = (a: number, b: number) => SortKeys[a] - SortKeys[b] || a - b;

/**
 * A dense fan as one label (the defaults audit; the author, 2026-11-20: "relayout, but precomputed so the UI doesn't
 * move"): runs of one layer's labels (the same `order`, each with a hit) where more than `over` fall within `span`
 * radians. A run is found sweeping the layer's labels by angle from just after its largest circular gap, so a run
 * across the origin is one (its middle angle modulo 2π). Writes `l.group` (a label's group number, or −1) and the
 * groups' records (`ng` of them: size, middle angle, best order, summed length, and `members`, their labels in the
 * order given, from `groupStart`); returns `ng`. Typed scratch only: nothing is made when no run is dense.
 */
function findGroups(l: RingLabels, over: number, span: number): number {
  l.ng = 0;
  const n = l.n;
  if (!(span > 0) || n <= over) return 0;
  l.sorted = grown(l.sorted, n);
  if (SortKeys.length < n) SortKeys = new Float64Array(Math.max(n, 2 * SortKeys.length));
  let m = 0;
  // the layer first, the angle within: the key's integer part is the layer's, its fraction the angle's share of the circle
  for (let i = 0; i < n; i++) {
    if (l.hits[i] === undefined) continue;
    SortKeys[i] = l.order[i] * 4 + (l.angle[i] / TwoPi) * 0.999999;
    l.sorted[m++] = i;
  }
  if (m <= over) return 0;
  l.group = grown(l.group, n);
  l.group.fill(-1, 0, n);
  const sorted = l.sorted.subarray(0, m);
  sorted.sort(byKey);
  for (let a = 0; a < m; ) {
    let e = a + 1;
    while (e < m && l.order[sorted[e]] === l.order[sorted[a]]) e++;
    const cnt = e - a;
    if (cnt > over) {
      // from just after the layer's largest circular gap, a sweep of runs
      let at = 0;
      let gap = -1;
      for (let k = 0; k < cnt; k++) {
        const x = l.angle[sorted[a + k]];
        const d = k + 1 < cnt ? l.angle[sorted[a + k + 1]] - x : l.angle[sorted[a]] + TwoPi - x;
        if (d > gap) {
          gap = d;
          at = (k + 1) % cnt;
        }
      }
      for (let k = 0; k < cnt; ) {
        const t0 = l.angle[sorted[a + ((at + k) % cnt)]];
        let r = 1;
        let last = t0;
        while (k + r < cnt) {
          const d = (l.angle[sorted[a + ((at + k + r) % cnt)]] - t0 + TwoPi) % TwoPi;
          if (d > span) break;
          last = t0 + d;
          r++;
        }
        if (r > over) {
          const g = l.ng++;
          l.groupSize = grown(l.groupSize, g + 1);
          l.groupAngle = grown(l.groupAngle, g + 1);
          l.groupOrder = grown(l.groupOrder, g + 1);
          l.groupLength = grown(l.groupLength, g + 1);
          let best = Infinity;
          let len = 0;
          for (let q = 0; q < r; q++) {
            const j = sorted[a + ((at + k + q) % cnt)];
            l.group[j] = g;
            best = Math.min(best, l.order[j]);
            len += Number.isFinite(l.length[j]) ? l.length[j] : 0;
          }
          l.groupSize[g] = r;
          l.groupAngle[g] = ((((t0 + last) / 2) % TwoPi) + TwoPi) % TwoPi;
          l.groupOrder[g] = best;
          l.groupLength[g] = len;
        }
        k += r;
      }
    }
    a = e;
  }
  if (l.ng === 0) return 0;
  // the members, flat by group (offsets in `groupStart`), each in the order the labels were given
  l.groupStart = grown(l.groupStart, l.ng + 1);
  l.groupStart[0] = 0;
  for (let g = 0; g < l.ng; g++) l.groupStart[g + 1] = l.groupStart[g] + l.groupSize[g];
  l.members = grown(l.members, l.groupStart[l.ng]);
  l.cursor = grown(l.cursor, l.ng);
  l.cursor.fill(0, 0, l.ng);
  for (let j = 0; j < n; j++) {
    const g = l.group[j];
    if (g >= 0) l.members[l.groupStart[g] + l.cursor[g]++] = j;
  }
  return l.ng;
}

/**
 * Lays a ring's labels out twice when a fan is dense: `p.fan` holds every label placed as if no group existed (what a
 * hover on a group opens: precomputed here, never laid out on a hover), `p.labels` the grouped labels, each group one
 * label at its middle angle with a hit of its own (`target: 'group'`, its members' hits). With no dense run, `p.labels`
 * alone, as before.
 */
function layoutLabels(p: RingPlacement, config: RenderConfig, size: number, text: TextMetrics) {
  const all = p.labels;
  const rc = config.rings;
  p.fan.reset();
  p.fanMax = rc.cluster.fanMax;
  all.ng = 0;
  if (findGroups(all, rc.cluster.over, (rc.cluster.degrees * Math.PI) / 180) === 0)
    return placeLabels(p, all, config, size, text);
  // the full set placed into `p.fan` (swapped in: no copy), then the grouped set built into `p.labels`
  const full = all;
  p.labels = p.fan;
  p.fan = full;
  const grouped = p.labels;
  grouped.reset();
  grouped.groupHits.length = 0;
  // a fan that opens pivots its members wider than their own arc about the group's middle, a line apart at the bend,
  // their leaders from their own angles: placed with the rest, then their own angles back
  const from = (p.fanFrom = grown(p.fanFrom, full.n));
  from.set(full.angle.subarray(0, full.n));
  pivotFans(p, full, labelLine(size) / (p.r + rc.labelGap + p.groupBand));
  placeLabels(p, full, config, size, text, from);
  full.angle.set(from.subarray(0, full.n));
  full.slot = grown(full.slot, full.n);
  full.slot.fill(-1, 0, full.n);
  for (let k = 0; k < full.m; k++) full.slot[full.placed[k]] = k;
  for (let j = 0; j < full.n; j++) {
    const g = full.group[j];
    if (g < 0) {
      grouped.add(full.text[j], full.angle[j], full.length[j], full.order[j], full.hits[j]);
      grouped.of = grown(grouped.of, grouped.n);
      grouped.of[grouped.n - 1] = -1;
    } else {
      // the group's label once, at its first member
      if (full.members[full.groupStart[g]] !== j) continue;
      const first = full.hits[j]!;
      const members: Hit[] = [];
      for (let q = full.groupStart[g]; q < full.groupStart[g + 1]; q++) members.push(full.hits[full.members[q]]!);
      const hit: Hit = {
        layer: first.layer,
        section: first.section,
        system: first.system,
        target: 'group',
        state: `group:${g}`,
        count: members.length,
        members,
      };
      grouped.groupHits[g] = hit;
      grouped.add(itemsText(members.length), full.groupAngle[g], full.groupLength[g], full.groupOrder[g], hit);
      grouped.of = grown(grouped.of, grouped.n);
      grouped.of[grouped.n - 1] = g;
    }
  }
  placeLabels(p, grouped, config, size, text);
}

/**
 * Pivots each group's members that a fan opens to (at most `p.fanMax`) about the group's middle angle, `step` radians
 * apart (a line at the bend) or their own spacing if wider, in their own order along the ring: a dozen leaders fan out
 * instead of running together over the group's own few degrees.
 */
function pivotFans(p: RingPlacement, l: RingLabels, step: number) {
  for (let g = 0; g < l.ng; g++) {
    const r = l.groupSize[g];
    if (r > p.fanMax || r < 2) continue;
    const c = l.groupAngle[g];
    const a = l.groupStart[g];
    // the members by their offset from the middle (small: a dozen; the scratch is the group's slice of `sorted`)
    const off = (j: number) => ((((l.angle[j] - c + Math.PI) % TwoPi) + TwoPi) % TwoPi) - Math.PI;
    const m = l.sorted;
    for (let q = 0; q < r; q++) m[q] = l.members[a + q];
    for (let q = 1; q < r; q++) {
      const j = m[q];
      let k = q - 1;
      while (k >= 0 && off(m[k]) > off(j)) {
        m[k + 1] = m[k];
        k--;
      }
      m[k + 1] = j;
    }
    const d = Math.max(step, (off(m[r - 1]) - off(m[0])) / (r - 1));
    for (let q = 0; q < r; q++) {
      const t = c + (q - (r - 1) / 2) * d;
      l.angle[m[q]] = ((t % TwoPi) + TwoPi) % TwoPi;
    }
  }
}

/**
 * Places a ring's labels: each side (right: angles in [0, π)) sorted by height, spread a line apart within the box
 * (a forward pass then a backward one), past what a side holds dropped by order then length; each label hugging the
 * circle at its leader's outer radius, its leader out along its angle to the bend, then to the label.
 */
function placeLabels(
  p: RingPlacement,
  l: RingLabels,
  config: RenderConfig,
  size: number,
  text: TextMetrics,
  from?: Float64Array,
) {
  const rc = config.rings;
  const line = labelLine(size);
  // past the groups' band: a leader never crosses a group's name
  const r1 = p.r + 2 + p.groupBand;
  const rb = p.r + rc.labelGap + p.groupBand;
  const rl = p.r + rc.labelGap + p.groupBand + rc.leaderBend;
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
      // from the mark (its own angle where a fan pivots it)
      const t0 = from ? from[i] : t;
      l.leader[o] = p.cx + r1 * Math.sin(t0);
      l.leader[o + 1] = p.cy - r1 * Math.cos(t0);
      l.leader[o + 2] = p.cx + rb * Math.sin(t);
      l.leader[o + 3] = p.cy - rb * Math.cos(t);
      l.leader[o + 4] = x;
      l.leader[o + 5] = y;
    }
  }
}

/**
 * The group whose fan is open for a hover (a hit), or −1: the hovered group's own hit, or the hit of one of its members.
 * A pure function of the placement and the hover: the draw (the emphasis the frame is drawn with) and the hit (the
 * viewer's hover model) both ask it, so a hit gives the same answer with no draw between.
 */
export function openGroupOf(p: RingPlacement, h: Hit | undefined): number {
  const f = p.fan;
  if (!h || f.m === 0 || f.ng === 0) return -1;
  const l = p.labels;
  if (h.target === 'group') {
    for (let q = 0; q < l.groupHits.length; q++) {
      const gh = l.groupHits[q];
      if (gh && gh.state === h.state && gh.layer === h.layer) return f.groupSize[q] > p.fanMax ? -1 : q;
    }
    return -1;
  }
  if (!h.element) return -1;
  const key = h.element.key;
  const owner = h.element.owner;
  for (let q = 0; q < f.groupStart[f.ng]; q++) {
    const e = f.hits[f.members[q]]?.element;
    if (e && e.key === key && e.owner === owner) {
      for (let r = 0; r < f.ng; r++)
        if (q >= f.groupStart[r] && q < f.groupStart[r + 1]) return f.groupSize[r] > p.fanMax ? -1 : r;
    }
  }
  return -1;
}

/**
 * The open group's fan, over the labels (the tip pass: a hover redraws it alone, nothing is laid out): the group is the
 * hover's (its hit, or one of its members'); its grouped label is covered by the page's ground and its members' labels
 * drawn where the layout placed them (`p.fan`), each on a ground of its own.
 */
function drawRingFan(
  p: RingPlacement,
  out: PassWriter,
  emphasis: Emphasis | undefined,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const f = p.fan;
  const g = openGroupOf(p, emphasis?.hover);
  if (g < 0) return;
  const l = p.labels;
  const ground = tokens.colors['msui-bg'];
  const ink = tokens.colors['msui-text'];
  const muted = tokens.colors['msui-text-muted'];
  const a = out.args;
  const size = p.labelSize;
  const half = labelLine(size) / 2;
  const cover = (set: RingLabels, k: number) => {
    const x0 = set.right[k] ? set.x[k] : set.x[k] - set.width[k];
    a[0] = x0 - 2;
    a[1] = set.y[k] - 0.35 * size - half;
    a[2] = set.width[k] + 4;
    a[3] = 2 * half;
    a[4] = ground;
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  };
  for (let k = 0; k < l.m; k++) if (l.of[l.placed[k]] === g) cover(l, k);
  const xy = p.path;
  for (let q = f.groupStart[g]; q < f.groupStart[g + 1]; q++) {
    const k = f.slot[f.members[q]];
    if (k < 0) continue;
    cover(f, k);
    const o = 6 * k;
    for (let c = 0; c < 6; c++) xy[c] = f.leader[o + c];
    a[4] = NoColor;
    a[5] = muted;
    a[6] = 1;
    a[7] = 0;
    out.pathArgs(xy, 3, false);
    a[0] = f.x[k];
    a[1] = f.y[k];
    a[2] = f.width[k];
    a[3] = RingMetrics.labelRoom;
    a[4] = ink;
    out.textArgs(f.text[f.members[q]], f.right[k] ? AnchorStart : AnchorEnd, config.itemFont);
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
  // each group's name past the arcs and the block labels at its arc's middle angle (`placeGroupNames`), upright, in
  // the label font and ink (the block labels are the item font, muted)
  for (let k = 0; k < p.groupN; k++) {
    g[0] = p.groupTX[k];
    g[1] = p.groupTY[k] + p.titleSize * 0.35;
    g[2] = 0;
    g[3] = Math.max(0, p.groupW[k]);
    g[4] = ink;
    out.textArgs(p.groupText[k], AnchorMiddle, config.labelFont);
  }
  if (!p.title) return;
  const size = p.titleSize;
  const hole = Math.max(0, p.r - p.depth);
  // a plate in the page's colour under the title and its subtitle: a pairs layer's chords cross the hole
  const w = Math.min(p.titleW, 2 * hole * 0.9) + 8;
  g[0] = p.cx - w / 2;
  g[1] = p.cy - size - 4;
  g[2] = w;
  g[3] = 2 * size + 10;
  g[4] = withAlpha(tokens.colors['msui-bg'], 0.85);
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 4;
  g[8] = 0;
  out.rectArgs();
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

const RingExtent = new Float64Array(4);

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
    // a narrow window grows about its middle (the window's floor, the rim's: no span on a ring's frame)
    const e = RingExtent;
    e[0] = x0;
    e[1] = x1;
    e[2] = Number.NEGATIVE_INFINITY;
    e[3] = Number.POSITIVE_INFINITY;
    floorExtent(e, 0, s.config.navigator.window);
    g[0] = e[0];
    g[1] = 0;
    g[2] = e[1] - e[0];
    g[3] = p.depth;
    g[4] = ((ink & 0xffffff00) | Math.round(WindowFill * 0.5 * 255)) >>> 0;
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
  drawWindowArgs(w, p.windowActive, lit, tokens, s.config.navigator.window);
  p.warp.close(out);
}

export type { NavContext, Viewport };
