/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Blocks in the tracks layout (_spec/README.md "Render", Blocks): a strip a block in view, each on its own system (its
// transform's), translated by nothing. `row` (the default): the blocks side by side along one axis, one viewport over
// them laid end to end, a gap between two (shrinking when zoomed out so positions keep three quarters of the width,
// the gaps in view sharing the rest, EC-102); `stack`: a block a row of the section's tracks, aligned at their first
// position, one zoom (EC-102, EC-115). Each block's label in a row above its tracks. Its lays are kept by the block
// and its first unit, so a pan keeps them. A one-position block is a chip: a fixed slot its label's width, never a
// position's share (zoom does not scale it), where it stands in a row, in a stack on one row the chips share.

import { blockUnits, piecesTotal, type ResolvedSection } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from '../axis.js';
import { fontSize, type RenderConfig } from '../config.js';
import { Anchor, grown, NoColor, type PassName, type PassWriter } from '../draw-lists.js';
import type { LettersLay } from '../letters.js';
import { cutLabel, longestCut } from '../range-labels.js';
import { elide, type TextMetrics } from '../text.js';
import { wellCorners } from '../wells.js';
import type { SectionLayout } from '../tracks.js';
import { PlacePoint, StripPlacement } from './strips.js';
import type { NavContext, PlaceContext } from './types.js';

/** A section's block as the resolved view gives it. */
export type ResolvedBlock = NonNullable<ResolvedSection['blocks']>[number];

// the anchor a number read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;

/**
 * Fixed sizes (CSS px): the gap between two blocks in a row at most (a boundary that reads), the share of the width
 * the gaps in view may take (a quarter: positions keep three quarters, however many blocks), a chip's padding either
 * side of its label (the text chip's rule), the share of the width (or the circle) the chips may take at most (half:
 * the polymers keep the rest) and the least a chip shrinks to past it.
 */
export const BlockMetrics = { gap: 12, gapShare: 0.25, chipPad: 5, chipShare: 0.5, chipMin: 16 } as const;

/**
 * Whether a block is a chip: its form `whole` (a molecule's default: a ligand, an ion, a glycan, every position of it
 * in one slot), or one position.
 */
export function isChip(block: {
  readonly form?: string;
  readonly axis: { readonly ranges?: readonly number[] };
}): boolean {
  return !!block.axis.ranges && (block.form === 'whole' || piecesTotal(block.axis.ranges) === 1);
}

/** A block's display units (a chip one, `blockUnits`; 0 while pending). */
export function blockLength(block: ResolvedBlock): number {
  return blockUnits(block) ?? 0;
}

/** A chip's width, CSS px: its label (its `×N` too) in the label font, `chipPad` either side; zoom never scales it. */
export function chipWidth(text: TextMetrics, config: RenderConfig, label: string): number {
  return text.measure(config.labelFont, label) + 2 * BlockMetrics.chipPad;
}

/**
 * Each block's chip slot (px; 0 for a block that is no chip) within `room` px: its `chipWidth`, all of them shrunk
 * together when they would take more than `chipShare` of the room, never below `chipMin` (their labels then elided).
 */
export function chipSlots(
  text: TextMetrics,
  config: RenderConfig,
  blocks: readonly ResolvedBlock[],
  room: number,
  out: Float64Array,
): Float64Array {
  const slots = grown(out, blocks.length);
  let sum = 0;
  for (let b = 0; b < blocks.length; b++) {
    slots[b] = isChip(blocks[b]) ? chipWidth(text, config, blocks[b].label) : 0;
    sum += slots[b];
  }
  const cap = BlockMetrics.chipShare * room;
  if (sum <= cap) return slots;
  // past the cap: each shrunk by one factor, never below the least (itself less when even the least of them all
  // would pass the cap): the factor whose slots sum to the cap, by bisection
  let n = 0;
  for (let b = 0; b < blocks.length; b++) if (slots[b] > 0) n++;
  const least = Math.min(BlockMetrics.chipMin, cap / n);
  let lo = 0;
  let hi = 1;
  for (let it = 0; it < 40; it++) {
    const k = (lo + hi) / 2;
    let s = 0;
    for (let b = 0; b < blocks.length; b++) if (slots[b] > 0) s += Math.max(least, slots[b] * k);
    if (s > cap) hi = k;
    else lo = k;
  }
  for (let b = 0; b < blocks.length; b++) if (slots[b] > 0) slots[b] = Math.max(least, slots[b] * lo);
  return slots;
}

/**
 * A chip's label as its pill shows it in a slot of `slot` px (less its padding), by the one rule of a label inside a mark
 * (`labelsInside`'s, layouts.md): whole; else cut with an ellipsis keeping at least `labels.minChars` letters (a label
 * of so many or fewer only whole); else none, the pill without text (its label in the label row and its card).
 */
export function chipText(text: TextMetrics, config: RenderConfig, label: string, slot: number): string {
  return fitLabel(text, config, label, slot - 2 * BlockMetrics.chipPad);
}

/** A label in `w` px by the rule of a label inside a mark: whole, cut keeping `labels.minChars` letters, or ''. */
function fitLabel(text: TextMetrics, config: RenderConfig, label: string, w: number): string {
  if (text.measure(config.labelFont, label) <= w) return label;
  const minChars = config.labels.minChars;
  if (label.length <= minChars) return '';
  const cut = longestCut(text, config.labelFont, label, w, Math.max(1, minChars));
  return cut ? cutLabel(label, cut) : '';
}

/**
 * A row's arithmetic (_spec/README.md "Render", Blocks): its blocks' offsets, lengths and chip slots over a data box,
 * placed for any window (`lay`): the frame's strips, a pan's and a zoom's inverse. Continuous in the window: a chip
 * partly in view takes px by the part of it in view, a gap by its boundary's distance to the nearer edge (none at the
 * edge, whole a gap in), so a pan moves nothing by a jump.
 */
export class BlockRow {
  n = 0;
  offset = new Float64Array(4);
  len = new Float64Array(4);
  /** A chip's px a unit (its slot over its positions: one but for a `whole` chip of several), 0 for no chip. */
  chip = new Float64Array(4);
  x0 = 0;
  width = 0;
  /** From `lay`: each block's units in view `[vs, ve)` (none when `ve <= vs`), its x and width, its gap's share. */
  vs = new Float64Array(4);
  ve = new Float64Array(4);
  x = new Float64Array(4);
  w = new Float64Array(4);
  f = new Float64Array(4);
  /** From `lay`: px a position of the blocks that are no chips, the window. */
  scale = 0;
  u0 = 0;
  u1 = 0;

  reset(n: number): void {
    this.n = n;
    this.offset = grown(this.offset, n);
    this.len = grown(this.len, n);
    this.vs = grown(this.vs, n);
    this.ve = grown(this.ve, n);
    this.x = grown(this.x, n);
    this.w = grown(this.w, n);
    this.f = grown(this.f, n);
  }

  /** Places the blocks for window `[u0, u1)`: into `vs`, `ve`, `x`, `w`, `scale`. */
  lay(u0: number, u1: number): void {
    this.u0 = u0;
    this.u1 = u1;
    const width = this.width;
    const span = Math.max(u1 - u0, 1e-9);
    let chipPx = 0;
    let chipU = 0;
    for (let b = 0; b < this.n; b++) {
      const o = this.offset[b];
      const vs = Math.max(u0, o);
      const ve = Math.min(u1, o + this.len[b]);
      this.vs[b] = vs;
      this.ve[b] = ve;
      if (this.chip[b] > 0 && ve > vs) {
        chipPx += this.chip[b] * (ve - vs);
        chipU += ve - vs;
      }
    }
    // a gap by its boundary's distance to the nearer edge at the chips' provisional scale: whole 12 px in
    const prov = Math.max(width - chipPx, 1e-9) / Math.max(span - chipU, 1e-9);
    let inner = 0;
    let cu = 0;
    let cp = 0;
    for (let b = 0; b < this.n; b++) {
      const o = this.offset[b];
      let f = 0;
      if (b > 0 && o > u0 && o < u1) {
        const left = (o - u0 - cu) * prov + cp;
        f = Math.min(1, Math.max(0, Math.min(left, width - left) / BlockMetrics.gap));
      }
      this.f[b] = f;
      inner += f;
      if (this.chip[b] > 0 && this.ve[b] > this.vs[b]) {
        cu += this.ve[b] - this.vs[b];
        cp += this.chip[b] * (this.ve[b] - this.vs[b]);
      }
    }
    // the gaps in view share a quarter of the width, at most 12 px each
    const gap = Math.min(BlockMetrics.gap, (BlockMetrics.gapShare * width) / Math.max(1, inner));
    const scale = Math.max(width - inner * gap - chipPx, 1e-9) / Math.max(span - chipU, 1e-9);
    this.scale = scale;
    let gaps = 0;
    let chipsU = 0;
    let chipsPx = 0;
    for (let b = 0; b < this.n; b++) {
      gaps += this.f[b] * gap;
      const vs = this.vs[b];
      const ve = this.ve[b];
      this.x[b] = this.x0 + (vs - u0 - chipsU) * scale + gaps + chipsPx;
      const chip = this.chip[b] > 0;
      this.w[b] = ve > vs ? (chip ? this.chip[b] * (ve - vs) : (ve - vs) * scale) : 0;
      if (chip && ve > vs) {
        chipsU += ve - vs;
        chipsPx += this.w[b];
      }
    }
  }

  /**
   * Each block in view's edges on the device grid (`dpr`), as laid (a chip's kept): a block's cells, its clip and the columns of another
   * row over the same window (a grid's under a tracks section's) then meet on whole device pixels. A block may be left
   * none wide (its strip is then not placed).
   */
  snap(dpr: number): void {
    const d = dpr > 0 ? dpr : 1;
    for (let b = 0; b < this.n; b++) {
      // a chip keeps its slot (its label's width) and its place
      if (this.ve[b] <= this.vs[b] || this.chip[b] > 0) continue;
      const a = Math.round(this.x[b] * d) / d;
      const e = Math.round((this.x[b] + this.w[b]) * d) / d;
      this.x[b] = a;
      this.w[b] = Math.max(0, e - a);
    }
  }

  /** The x of section unit `u` as last laid (past the window: on at the polymers' scale). */
  xOf(u: number): number {
    if (u < this.u0) return this.x0 - (this.u0 - u) * this.scale;
    if (u > this.u1) return this.x0 + this.width + (u - this.u1) * this.scale;
    for (let b = 0; b < this.n; b++) {
      const vs = this.vs[b];
      const ve = this.ve[b];
      if (ve <= vs || u < vs || (u >= ve && u < this.u1)) continue;
      return this.x[b] + (u - vs) * (this.chip[b] > 0 ? this.chip[b] : this.scale);
    }
    return this.x0 + this.width;
  }

  /** The section unit at x as last laid (in a gap: the next block's first). */
  unitAt(x: number): number {
    for (let b = 0; b < this.n; b++) {
      if (this.ve[b] <= this.vs[b]) continue;
      if (x < this.x[b]) return this.vs[b];
      if (x < this.x[b] + this.w[b])
        return this.vs[b] + (x - this.x[b]) / (this.chip[b] > 0 ? this.chip[b] : this.scale);
    }
    return this.u1;
  }

  /** The px from x to the next block's first when x is in a gap as last laid (0 on a block). */
  gapAfter(x: number): number {
    for (let b = 0; b < this.n; b++) {
      if (this.ve[b] <= this.vs[b]) continue;
      if (x < this.x[b]) return this.x[b] - x;
      if (x < this.x[b] + this.w[b]) return 0;
    }
    return 0;
  }

  /** The window's start, `w` units wide, that puts unit `u` at `x`: a bisection (x of `u` falls as the start grows). */
  startFor(w: number, u: number, x: number): number {
    let lo = u - w;
    let hi = u;
    for (let k = 0; k < 52 && hi - lo > 1e-9; k++) {
      const mid = (lo + hi) / 2;
      this.lay(mid, mid + w);
      if (this.xOf(u) > x) lo = mid;
      else hi = mid;
    }
    return (lo + hi) / 2;
  }
}

/** A section of blocks placed: its strips by x (a row) or by y (a stack), each strip's units its block's. */
export class BlockPlacement extends StripPlacement {
  stacked = false;
  /** The label row above each block's tracks: its height and the labels' baseline in it. */
  labelRow = 0;
  /** Each strip's block's label (and the block's first unit: `strips.first`). */
  labels: string[] = [];
  /** Each strip's label as drawn: elided to its strip's width ('' when nothing fits: dropped); a chip's in its slot. */
  shown: string[] = [];
  /** A row's arithmetic (pan and zoom through its inverse). */
  readonly row = new BlockRow();
  private slots: Float64Array = new Float64Array(4);

  /** Whether strip `i` starts its row (a stack's chips share one): its row's fill and track labels drawn once. */
  rowStart(i: number): boolean {
    const s = this.strips;
    return i === 0 || s.top[i] !== s.top[i - 1];
  }

  /** The strip at a page point: by x in a row, by its rows' extent in a stack (the chips' row by x too); −1 when none. */
  private stripOf(x: number, y: number): number {
    const s = this.strips;
    if (!this.stacked) return s.atX(x);
    let row = -1;
    for (let i = 0; i < s.n; i++) {
      if (y < s.top[i] || y >= s.top[i] + s.height[i]) continue;
      if (!s.chip[i]) return i;
      if (x >= s.axes[i].x0 && x < s.axes[i].x0 + s.axes[i].width) return i;
      if (row < 0) row = i;
    }
    return row;
  }

  override at(p: ArrayLike<number>, out: PlacePoint): boolean {
    const x = p[0];
    const y = p[1];
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height) return false;
    const s = this.strips;
    const i = this.stripOf(x, y);
    const first = s.n ? s.axes[0].x0 : 0;
    if (i < 0) {
      out.kind = x < first ? 'gutter' : 'gap';
      return true;
    }
    const a = s.axes[i];
    out.strip = i;
    out.x = x;
    out.y = y;
    for (let j = 0; j < s.bands; j++) {
      const k = i * s.bands + j;
      if (s.bandH[k] > 0 && y >= s.bandY[k] && y < s.bandY[k] + s.bandH[k]) {
        out.band = j;
        break;
      }
    }
    if (x < a.x0) out.kind = 'gutter';
    else if (x >= a.x0 + a.width) out.kind = this.stacked ? 'gap' : 'gutter';
    else if (out.band < 0) out.kind = 'gap';
    else {
      out.kind = 'data';
      out.offset = s.first[i] + a.unitAt(x);
    }
    return true;
  }

  override offsetAt(x: number, y: number): number {
    const s = this.strips;
    if (!s.n) return -1;
    let i = this.stripOf(x, y);
    if (i < 0 && !this.stacked) i = s.nearestX(x);
    else if (i < 0) {
      // a stack's nearest row by y
      let best = Infinity;
      for (let k = 0; k < s.n; k++) {
        const d = Math.min(Math.abs(y - s.top[k]), Math.abs(y - s.top[k] - s.height[k]));
        if (d < best) [best, i] = [d, k];
      }
    }
    return s.offsetIn(i, x);
  }

  override place(offset: number, out: PlacePoint): boolean {
    return this.strips.placeAcross(offset, out);
  }

  /** The row's chip slots for the section's blocks (kept: a placement a frame). */
  slotsFor(ctx: PlaceContext, blocks: readonly ResolvedBlock[], room: number): Float64Array {
    this.slots = chipSlots(ctx.text, ctx.config, blocks, room, this.slots);
    return this.slots;
  }
}

/**
 * A pan of a row through its inverse (_spec/README.md "Render", Blocks): the unit `dx` px right of a reference point
 * comes to it (the window's middle less half the move, both in view), so what is under the pointer moves with it over
 * chips and gaps; undefined past a box's width (the caller pans linearly). With the pointer's px `x` (from the axis's
 * left) the reference is the pointer, within the box the move keeps in view: the content under it follows it, not the
 * window's middle.
 */
export function rowPan(row: BlockRow, v: Viewport, dx: number, x?: number): Viewport | undefined {
  if (!(Math.abs(dx) < row.width) || !dx) return undefined;
  const w = v.end - v.start;
  row.lay(v.start, v.end);
  let ref = row.x0 + row.width / 2 - dx / 2;
  if (x !== undefined && Number.isFinite(x)) {
    const lo = row.x0 + Math.max(0, -dx);
    const hi = row.x0 + row.width - Math.max(0, dx);
    ref = Math.min(Math.max(row.x0 + x, lo), Math.max(lo, hi));
  }
  // a point in a gap has no unit of its own: both points moved on to the next block's first px
  ref += row.gapAfter(ref + dx);
  const u = row.unitAt(ref + dx);
  const start = row.startFor(w, u, ref);
  return { ...v, start, end: start + w };
}

/** A zoom of a row about unit `at`, its px kept where it is (chips entering or leaving the view move it not). */
export function rowZoom(row: BlockRow, c: NavContext, v: Viewport, factor: number, at: number): Viewport {
  row.lay(v.start, v.end);
  const x = row.xOf(at);
  const total = c.total;
  const w = Math.min(Math.max((v.end - v.start) / factor, Math.min(c.width / c.maxScale, total)), total);
  const start = row.startFor(w, at, x);
  return { ...v, start, end: start + w };
}

/** Places a section's blocks in a row or a stack over the section's window (`ctx.axis`'s). */
export function placeBlocks(section: ResolvedSection, ctx: PlaceContext, prev?: BlockPlacement): BlockPlacement {
  const p = prev ?? new BlockPlacement();
  const blocks = section.blocks!;
  const { config, axis } = ctx;
  const strips = p.strips;
  const nb = ctx.bands.length;
  strips.reset(nb);
  p.labels.length = 0;
  p.shown.length = 0;
  p.stacked = section.params.arrange === 'stack';
  p.labelRow = fontSize(config.labelFont) + 4;
  const top = ctx.box.y;
  const u0 = axis.u0;
  const u1 = axis.u1;
  const x0 = axis.x0;
  const width = axis.width;
  // the chips' slots: at most half the width together
  const slots = p.slotsFor(ctx, blocks, width);
  const row = p.row;
  row.reset(blocks.length);
  row.x0 = x0;
  row.width = width;
  for (let b = 0; b < blocks.length; b++) {
    row.offset[b] = blocks[b].offset ?? 0;
    row.len[b] = blockLength(blocks[b]);
  }
  row.chip = grown(row.chip, blocks.length);
  // a chip's slot shared by its positions (a `whole` chip of several: each a part of its one slot)
  for (let b = 0; b < blocks.length; b++) row.chip[b] = slots[b] > 0 ? slots[b] / Math.max(1, row.len[b]) : 0;
  // a row's gaps and chips take their px from the width (the row's x is not linear over them, EC-102)
  // the blocks' edges on the device grid, as a grid's under them (its columns align with these)
  if (!p.stacked) {
    row.lay(u0, u1);
    row.snap(ctx.dpr);
  }
  // a stack: no gaps, one scale
  const scale = Math.max(width, 1e-9) / (u1 - u0);
  // the rows' height: the label row, then the tracks
  let rows = p.labelRow;
  const pads = ctx.pads;
  for (let j = 0; j < nb; j++)
    rows += (j ? config.trackGap : 0) + ctx.heights[j] + (pads ? pads[2 * j] + pads[2 * j + 1] : 0);
  let y = top;
  // a stack: the other blocks a row each, then the chips' one row, side by side from the first position's x
  let stackRows = 0;
  let chipRow = -1;
  let chipX = x0;
  for (let round = 0; round < (p.stacked ? 2 : 1); round++)
    for (let b = 0; b < blocks.length; b++) {
      const block = blocks[b];
      const chip = slots[b] > 0;
      if (p.stacked && chip !== (round === 1)) continue;
      let rowTop = y;
      if (p.stacked && (!chip || chipRow < 0)) {
        if (stackRows++) y += config.sectionGap;
        rowTop = y;
        y += rows;
        if (chip) chipRow = rowTop;
      } else if (p.stacked) rowTop = chipRow;
      if (!block.axis.ranges) continue;
      const o = p.stacked ? 0 : row.offset[b];
      let a: AxisTransform;
      if (p.stacked && chip) {
        // a stack's chips are on their row whatever the window, each its slot over its one position
        a = AxisTransform.on(block.axis, block.axis.ranges, chipX, slots[b], ctx.dpr, { start: 0, end: row.len[b] });
        chipX += slots[b] + BlockMetrics.gap;
      } else if (p.stacked) {
        const vs = Math.max(u0, 0);
        const ve = Math.min(u1, row.len[b]);
        if (ve <= vs) continue;
        a = AxisTransform.on(block.axis, block.axis.ranges, x0 + (vs - u0) * scale, (ve - vs) * scale, ctx.dpr, {
          start: vs,
          end: ve,
        });
      } else {
        const vs = row.vs[b];
        const ve = row.ve[b];
        // in view but none wide once snapped: no strip (nothing of it is drawn or hit)
        if (ve <= vs || row.w[b] <= 0) continue;
        // a chip partly in view takes px by the part of it in view (a pan moves it continuously)
        a = AxisTransform.on(block.axis, block.axis.ranges, row.x[b], row.w[b], ctx.dpr, {
          start: vs - o,
          end: ve - o,
        });
      }
      const i = strips.add(a, b, 0, 0, o);
      if (chip) strips.chip[i] = 1;
      p.labels[i] = block.label;
      // a label that does not fit its block is elided, or dropped when nothing of it does; a chip's by the pill's rule
      // (whole, or cut keeping `labels.minChars` letters, or none)
      p.shown[i] = chip
        ? fitLabel(ctx.text, config, block.label, Math.min(a.width, slots[b]))
        : elide(ctx.text, config.labelFont, block.label, a.width);
      strips.top[i] = rowTop;
      strips.height[i] = rows;
      let by = rowTop + p.labelRow;
      for (let j = 0; j < nb; j++) {
        if (j) by += config.trackGap;
        // a band inside its well (_spec/README.md "Render", Tracks layout)
        const pt = pads ? pads[2 * j] : 0;
        const pb = pads ? pads[2 * j + 1] : 0;
        strips.padT[i * nb + j] = pt;
        strips.padB[i * nb + j] = pb;
        strips.bandY[i * nb + j] = by + pt;
        strips.bandH[i * nb + j] = ctx.heights[j];
        by += pt + ctx.heights[j] + pb;
      }
    }
  const height = p.stacked ? y - top : rows;
  p.contentHeight = height;
  p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height };
  for (let i = 0; i < strips.n; i++) ctx.layStrip(p, i);
  return p;
}

/** A block's label above its tracks (the labels pass), faded at its strip's width; a stack's rows filled as tracks. */
export function blockChrome(
  s: SectionLayout,
  p: BlockPlacement,
  pass: PassName,
  out: PassWriter,
  strip: number,
  band: number,
): void {
  const strips = p.strips;
  if (pass === 'data' && band >= 0 && strip > 0 && p.stacked) {
    // a row's ground once, before its first strip's layers: never over a strip drawn before it (a chip's pill)
    if (!p.rowStart(strip)) return;
    const k = strip * strips.bands + band;
    const well = s.tracks[band]?.well;
    // its well, as the first row's (_spec/README.md "Render", Shapes): through `args`, nothing boxed
    const g = out.args;
    g[0] = well?.x ?? 0;
    g[1] = strips.bandY[k] - strips.padT[k];
    g[2] = well?.width ?? 0;
    g[3] = strips.bandH[k] + strips.padT[k] + strips.padB[k];
    wellCorners(g, 7, s.config);
    g[4] = s.tokens.colors['seq-track'];
    g[5] = NoColor;
    g[6] = 0;
    out.rectArgs();
    return;
  }
  if (pass !== 'labels' || band >= 0) return;
  const g = out.args;
  // the label font's size as the layout read it (its row less 4 px): a draw parses no font (a regex's match allocates)
  const size = p.labelRow - 4;
  for (let i = 0; i < strips.n; i++) {
    // a chip's label is in its letters' pill where the pill is drawn
    if (strips.chip[i] && hasPill(s, i)) continue;
    const label = p.shown[i];
    if (!label) continue;
    const a = strips.axes[i];
    g[0] = a.x0;
    g[1] = strips.top[i] + size;
    g[2] = 0;
    g[3] = Math.max(0, a.width);
    g[4] = s.tokens.colors['msui-text-muted'];
    out.textArgs(label, AnchorStart, s.config.labelFont);
  }
}

/** Whether a chip's strip draws its label in a pill: a letters layer that drew one on a band with height. */
export function hasPill(s: SectionLayout, strip: number): boolean {
  const nb = s.placement.strips.bands;
  for (let j = 0; j < nb; j++) {
    const k = strip * nb + j;
    if (!(s.placement.strips.bandH[k] > 0)) continue;
    const lays = s.lays[k];
    for (let i = 0; lays && i < lays.length; i++) {
      if (lays[i].type.name !== 'letters') continue;
      // drawn: a chip's label and a letter laid out (none while its tile is pending, at bars, or with no letter)
      const l = lays[i].lay as LettersLay;
      if (l.chip !== '' && l.n > 0) return true;
    }
  }
  return false;
}
