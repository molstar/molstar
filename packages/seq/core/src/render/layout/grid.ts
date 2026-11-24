/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The grid layout (_spec/README.md "Render", The grid layout): a section of two axes, one strip, the plot, whose
// transform is x and which carries a second transform `strip.y` down its frame (read by a layer as `ctx.axisY`). Its
// tracks' bands overlay the plot (later tracks over earlier, the last one's layers hit first). Its x is the data box
// (the view's x); its height is the data box's width over `aspect` (1: square), capped by the view's room below it
// less a tenth (that room never below `heights.grid`), at least two rows, and fixed: no zoom or pan of x or y changes it
// (cells are rectangles when the two show different spans); `aspect: none` leaves it to the tracks' heights. Its viewport is the tracks' window plus `y` when y has a window of its own. A grid over a
// section of blocks (its x and y the blocks' display axis, the chains side by side) is a strip a pair of blocks in view:
// each its x block's transform and its y block's (the blocks laid along each axis as a row of blocks is, gaps between
// them; y's a second section's blocks over its window when `follow_y` names one, EC-209), its bands that pair's cell, its cells clipped at the block's edges; a chain boundary is the blocks' gap with no
// line in it (a hairline only where two blocks touch with no gap), the blocks' labels on both axes.

import { AxisTransform, type Viewport } from '../axis.js';
import type { RenderConfig } from '../config.js';
import { Anchor, grown, type PassName, type PassWriter } from '../draw-lists.js';
import type { SectionLayout } from '../tracks.js';
import { BlockRow, type ResolvedBlock, rowBlocks } from './blocks.js';
import type { MemberLabels } from './grid-members.js';
import { MemberHooks } from '../members-hooks.js';
import { type PlacePoint, StripPlacement } from './strips.js';
import { clampWindow, sameWindow, TracksNav } from './tracks.js';
import type { ResolvedSection } from '@molstar/seqx-core/view';
import type { NavContext, PlaceContext, SectionLayoutType, SectionNav } from './types.js';

// the anchors as numbers read once: a draw reads no property of an imported object
const AnchorStart: number = Anchor.start;
const AnchorEnd: number = Anchor.end;

/** The least plot height: two rows (CSS px). */
const MinHeight = 40;

/** The plot's share of the view's height at most: its room less a tenth. */
const RoomShare = 0.9;

/** The grid's placement: its strips', its bands all the plot (`at` gives the last band: the topmost track). */
export class GridPlacement extends StripPlacement {
  /** The y axis's first and last positions in view, 1-based, as drawn at the plot's top and bottom (made at place). */
  yFirst = '';
  yLast = '';
  // overlaid bands: `at` gives the last with height (the topmost track)
  override readonly overlay = true;
  /** A grid over blocks: the blocks laid along x and along y (null otherwise). */
  xRow: BlockRow | null = null;
  yRow: BlockRow | null = null;
  /** Each y block's transform down the plot (made once a place, reused by every x block's strips). */
  readonly yAxes: (AxisTransform | undefined)[] = [];
  /** The plot's top and height on the page: its well's less the padding (_spec/README.md "Render", Tracks layout). */
  plotTop = 0;
  plotHeight = 0;
  /** The padding above and below the plot inside its well. */
  padTop = 0;
  padBottom = 0;
  /** Categorical axes' labels (EC-35): x's room above the plot, their advances (made by the chunk when it lands). */
  members: MemberLabels | undefined;

  override at(p: ArrayLike<number>, out: PlacePoint): boolean {
    if (!this.xRow) return super.at(p, out);
    const x = p[0];
    const y = p[1];
    const b = this.box;
    out.strip = -1;
    out.band = -1;
    out.kind = 'empty';
    out.offset = -1;
    if (y < b.y || y >= b.y + b.height) return false;
    const s = this.strips;
    for (let i = 0; i < s.n; i++) {
      const a = s.axes[i];
      if (x < a.x0 || x >= a.x0 + a.width || y < s.top[i] || y >= s.top[i] + s.height[i]) continue;
      out.strip = i;
      out.x = x;
      out.y = y;
      // overlaid: the topmost band with height
      for (let j = s.bands - 1; j >= 0; j--)
        if (s.bandH[i * s.bands + j] > 0) {
          out.band = j;
          break;
        }
      out.kind = out.band < 0 ? 'gap' : 'data';
      if (out.band >= 0) out.offset = s.first[i] + a.unitAt(x) * s.scale[i];
      return true;
    }
    out.kind = s.n && x < s.axes[0].x0 ? 'gutter' : 'gap';
    return true;
  }

  override offsetAt(x: number, y: number): number {
    if (!this.xRow) return super.offsetAt(x, y);
    const r = this.xRow;
    return r.unitAt(Math.min(Math.max(x, r.x0), r.x0 + r.width - 1e-6));
  }
}

/** A row of a section's blocks over `[x0, x0 + width)` (a grid's axis across or down) for window `[u0, u1)`. */
function blockRow(
  row: BlockRow,
  blocks: readonly ResolvedBlock[],
  x0: number,
  width: number,
  u0: number,
  u1: number,
  dpr: number,
  config: RenderConfig,
) {
  // the tracks' row's rules (groups, gaps, runs), so the plot's columns align with the blocks above
  rowBlocks(row, blocks, config);
  // no chips on a plot's axes: every block a position's share
  row.chip = grown(row.chip, blocks.length);
  row.chip.fill(0);
  row.x0 = x0;
  row.width = width;
  row.lay(u0, u1);
  // each block's edges on the device grid, as the tracks' row above: its cells' clip and its tiles' edges meet exactly
  row.snap(dpr);
}

/**
 * A grid over a section of blocks (_spec/README.md "Render", The grid layout): a strip a pair of blocks in view, x's
 * block across and y's down, each laid as a row of blocks over its window (y follows x: the same window down the
 * plot's height), the pair's cell its bands' box.
 */
function placeBlocksGrid(
  section: ResolvedSection,
  ctx: PlaceContext,
  p: GridPlacement,
  top: number,
  h: number,
  pt: number,
  pb: number,
) {
  const blocks = section.blocks!;
  // y over a second section of blocks (`follow_y`, EC-209), its window that section's; else x's blocks and window
  const yBlocks = gridYBlocks(section);
  const strips = p.strips;
  const x = ctx.axis;
  const xr = (p.xRow ??= new BlockRow());
  const yr = (p.yRow ??= new BlockRow());
  blockRow(xr, blocks, x.x0, x.width, x.u0, x.u1, ctx.dpr, ctx.config);
  const yv = ctx.y?.viewport ?? { start: x.u0, end: x.u1 };
  blockRow(yr, yBlocks, top, h, yv.start, yv.end, ctx.dpr, ctx.config);
  const ys = p.yAxes;
  ys.length = yBlocks.length;
  for (let c = 0; c < yBlocks.length; c++) {
    const oc = yr.offset[c];
    const yb = yBlocks[c];
    ys[c] =
      yr.ve[c] > yr.vs[c] && yr.w[c] > 0 && yb.axis.ranges
        ? AxisTransform.on(yb.axis, yb.axis.ranges!, yr.x[c], yr.w[c], ctx.dpr, {
            start: (yr.vs[c] - oc) / (yb.scale ?? 1),
            end: (yr.ve[c] - oc) / (yb.scale ?? 1),
          })
        : undefined;
  }
  for (let b = 0; b < blocks.length; b++) {
    // a block in view whose snapped width is none has no strip: nothing of it is drawn or hit
    if (xr.ve[b] <= xr.vs[b] || xr.w[b] <= 0 || !blocks[b].axis.ranges) continue;
    const ob = xr.offset[b];
    const sb = blocks[b].scale ?? 1;
    const ax = AxisTransform.on(blocks[b].axis, blocks[b].axis.ranges!, xr.x[b], xr.w[b], ctx.dpr, {
      start: (xr.vs[b] - ob) / sb,
      end: (xr.ve[b] - ob) / sb,
    });
    for (let c = 0; c < yBlocks.length; c++) {
      const ay = ys[c];
      if (!ay) continue;
      const oc = yr.offset[c];
      const i = strips.add(ax, b, 0, 0, ob);
      strips.scale[i] = sb;
      strips.ys[i] = ay;
      strips.yBlock[i] = c;
      strips.yFirst[i] = oc;
      strips.top[i] = yr.x[c];
      strips.height[i] = yr.w[c];
      // the plot's well round its cells: the padding above the first row and below the last
      const padT = yr.x[c] <= top + 0.5 ? pt : 0;
      const padB = yr.x[c] + yr.w[c] >= top + h - 0.5 ? pb : 0;
      for (let j = 0; j < ctx.bands.length; j++) {
        strips.bandY[i * strips.bands + j] = yr.x[c];
        strips.bandH[i * strips.bands + j] = yr.w[c];
        strips.padT[i * strips.bands + j] = padT;
        strips.padB[i * strips.bands + j] = padB;
      }
    }
  }
}

/** A grid over blocks' y blocks: a second section's (`follow_y`), else its x blocks. */
export function gridYBlocks(section: ResolvedSection): readonly ResolvedBlock[] {
  return section.y?.blocks ?? section.blocks!;
}

/** A grid's y window: its own, clamped within its units and at least a unit. */
function clampY(total: number, y: { start: number; end: number }): { start: number; end: number } {
  const w = Math.min(Math.max(y.end - y.start, 1), total);
  const s = Math.min(Math.max(y.start, 0), total - w);
  return s === y.start && s + w === y.end ? y : { start: s, end: s + w };
}

/**
 * The grid's arithmetic: x as the tracks' window; y, when it has a window of its own (`v.y`, `NavContext.totalY`),
 * panned by dy (px of the plot's height), zoomed with x about its centre, fit whole.
 */
export const GridNav: SectionNav = {
  pan(c, v, dx, dy) {
    const x = clampWindow(c, {
      ...v,
      start: v.start + (dx * (v.end - v.start)) / c.width,
      end: v.end + (dx * (v.end - v.start)) / c.width,
    });
    if (!v.y || c.totalY === undefined || !dy) return x;
    const d = (dy * (v.y.end - v.y.start)) / Math.max(1, c.height);
    return { ...x, y: clampY(c.totalY, { start: v.y.start + d, end: v.y.end + d }) };
  },
  zoom(c, v, factor, at) {
    const x = TracksNav.zoom(c, v, factor, at);
    if (!v.y || c.totalY === undefined) return x;
    const m = (v.y.start + v.y.end) / 2;
    return { ...x, y: clampY(c.totalY, { start: m - (m - v.y.start) / factor, end: m + (v.y.end - m) / factor }) };
  },
  fit(c, v) {
    const x = { ...v, start: 0, end: c.total };
    return c.totalY !== undefined && v.y ? { ...x, y: { start: 0, end: c.totalY } } : x;
  },
  reveal: TracksNav.reveal,
  step: TracksNav.step,
  clamp(c: NavContext, v: Viewport) {
    const x = clampWindow(c, v);
    if (!v.y || c.totalY === undefined) return x;
    const y = clampY(c.totalY, v.y);
    return y === v.y ? x : { ...x, y };
  },
};

/** Two grid viewports place alike: their windows, and their own y windows. */
function sameGrid(a: Viewport | undefined, b: Viewport | undefined): boolean {
  if (!sameWindow(a, b)) return false;
  const ya = a?.y;
  const yb = b?.y;
  return ya === yb || (!!ya && !!yb && ya.start === yb.start && ya.end === yb.end);
}

export const GridLayout: SectionLayoutType<GridPlacement> = {
  name: 'grid',
  fillBands: true,
  axesY: true,
  place(section, ctx: PlaceContext, _viewport, prev) {
    const p = prev instanceof GridPlacement ? prev : new GridPlacement();
    const strips = p.strips;
    const { config } = ctx;
    strips.reset(ctx.bands.length);
    // the plot inside its well: the bands' largest padding above and below (they share the plot)
    let pt = 0;
    let pb = 0;
    for (let j = 0; ctx.pads && j < ctx.bands.length; j++) {
      pt = Math.max(pt, ctx.pads[2 * j]);
      pb = Math.max(pb, ctx.pads[2 * j + 1]);
    }
    // a categorical x: its labels' room above the plot, the box grown by it (fixed: no zoom moves it)
    const c = MemberHooks.chrome;
    if (c && !p.members && (ctx.axis.members || ctx.y?.members)) p.members = new c.MemberLabels();
    const room = section.blocks || !c || !p.members ? 0 : c.placeMemberLabels(p.members, ctx.axis, ctx.text, config);
    const well = ctx.box.y;
    const top = well + room + pt;
    const width = ctx.axis.width;
    const aspect = section.params.aspect;
    let h: number;
    const x = ctx.axis;
    const y = ctx.y;
    if (aspect === 'none') {
      h = 0;
      for (const t of ctx.heights) h = Math.max(h, t);
    } else {
      const a = typeof aspect === 'number' && aspect > 0 ? aspect : 1;
      // the room below the section's top less a tenth, never below a grid track's default height: a view sized to its
      // content (whose height follows the plot's) keeps a stable plot
      const room = Math.max(config.heights.grid, (ctx.height - top) * RoomShare);
      // x is the view's x (the data box, aligned with the sections it shares x with); the box is the data box's width
      // over `a`, capped by the room, and never follows the windows (a plot's box is fixed): a zoom or a pan of either
      // axis leaves it, and every track's place, where it is; cells are rectangles when x and y show different spans
      h = Math.max(MinHeight, Math.min(width / a, room));
    }
    p.yFirst = p.yLast = '';
    p.plotTop = top;
    p.plotHeight = h;
    p.padTop = pt;
    p.padBottom = pb;
    if (section.blocks) {
      placeBlocksGrid(section, ctx, p, top, h, pt, pb);
      p.contentHeight = h + pt + pb;
      p.box = { x: config.inset, y: well, width: ctx.width - 2 * config.inset, height: h + pt + pb };
      for (let k = 0; k < strips.n; k++) ctx.layStrip(p, k);
      return p;
    }
    p.xRow = p.yRow = null;
    const i = strips.add(x, 0, 0, 0, x.u0);
    for (let j = 0; j < ctx.bands.length; j++) {
      strips.bandY[i * strips.bands + j] = top;
      strips.bandH[i * strips.bands + j] = h;
      strips.padT[i * strips.bands + j] = pt;
      strips.padB[i * strips.bands + j] = pb;
    }
    if (y?.pieces.length) {
      const ya = AxisTransform.on(y.axis, y.pieces, top, h, ctx.dpr, y.viewport);
      if (y.members) ya.members = y.members;
      strips.ys[i] = ya;
      // a categorical y names its members down the label column instead
      const first = ya.members ? -1 : ya.elementAtUnit(ya.u0);
      const last = ya.members ? -1 : ya.elementAtUnit(Math.max(ya.u0, ya.u1 - 1e-6));
      if (first >= 0) p.yFirst = String(first + 1);
      if (last >= 0) p.yLast = String(last + 1);
    }
    strips.top[i] = well + room;
    strips.height[i] = h + pt + pb;
    p.contentHeight = room + h + pt + pb;
    p.box = { x: config.inset, y: well, width: ctx.width - 2 * config.inset, height: room + h + pt + pb };
    ctx.layStrip(p, i);
    return p;
  },
  sameKey: sameGrid,
  nav: GridNav,
  chrome(s, pass, out, _emphasis, strip, band) {
    if (pass === 'background' && strip === -1) {
      tintGaps(s, out);
      return;
    }
    if (pass !== 'labels' || strip !== -1) return;
    drawChrome(s, pass, out);
  },
  hitChrome(s, x, y) {
    const p = s.placement;
    const st = p.strips;
    if (!(p instanceof GridPlacement) || p.xRow || !st.n || !p.members) return undefined;
    return MemberHooks.chrome?.hitMemberLabel(s, p.members, st.axes[0], st.ys[0], p.plotTop, p.plotHeight, x, y);
  },
};

/**
 * The plot's labels: the first track's label at its top in the label column, the y axis's first and last positions
 * (1-based) at the label column's right edge, top and bottom, muted.
 */
function drawChrome(s: SectionLayout, _pass: PassName, out: PassWriter) {
  const strips = s.placement.strips;
  if (!strips.n) return;
  const { config, tokens } = s;
  const g = out.args;
  const t = s.tracks[0];
  const p = s.placement;
  if (p instanceof GridPlacement && p.xRow && p.yRow) {
    drawBlocksChrome(s, p, out);
    return;
  }
  const top = p instanceof GridPlacement ? p.plotTop : strips.top[0];
  // categorical axes (EC-35): the members' labels; the track's label above the plot then, where the x labels leave room
  const ay = strips.ys[0];
  const members = p instanceof GridPlacement && (strips.axes[0].members || ay?.members);
  if (members && p.members) MemberHooks.chrome?.drawMemberLabels(s, p.members, strips.axes[0], ay, p.plotTop, out);
  const labelAt = members && ay?.members ? (p.box.y + 12 <= top - 4 ? p.box.y + 12 : -1) : top + 12;
  if (t?.labelText && labelAt >= 0) {
    g[0] = config.inset;
    g[1] = labelAt;
    g[2] = t.labelWidth;
    g[3] = config.labelWidth;
    g[4] = tokens.colors['msui-text-muted'];
    out.textArgs(t.labelText, AnchorStart, config.labelFont);
  }
  if (!(p instanceof GridPlacement)) return;
  const right = s.axis ? s.axis.x0 - config.labelGap / 2 : config.inset + config.labelWidth;
  for (let k = 0; k < 2; k++) {
    const text = k === 0 ? p.yFirst : p.yLast;
    if (!text) continue;
    g[0] = right;
    g[1] = k === 0 ? top + 28 : top + p.plotHeight - 4;
    g[2] = 0;
    g[3] = config.labelWidth;
    g[4] = tokens.colors['msui-text-muted'];
    out.textArgs(text, AnchorEnd, config.itemFont);
  }
}

/**
 * A grid over two sections of blocks (`follow_y`): the gaps between its blocks, down the plot between x blocks and
 * across it between y blocks, tinted a step off the cells' ground (`msui-line`: darker in light, lighter in dark) under
 * the data, so a block-pair cell reads as a cell with no line in a gap (the single-section grid, the PAE plot, has none).
 * Each gap is exactly the space between two blocks in view.
 */
function tintGaps(s: SectionLayout, out: PassWriter) {
  const p = s.placement;
  if (!(p instanceof GridPlacement) || !p.xRow || !p.yRow || !s.section.y?.blocks) return;
  const xr = p.xRow;
  const yr = p.yRow;
  const ink = s.tokens.colors['msui-line'];
  const top = p.plotTop;
  const h = p.plotHeight;
  for (let b = 1; b < s.section.blocks!.length; b++) {
    if (xr.ve[b] <= xr.vs[b] || xr.ve[b - 1] <= xr.vs[b - 1]) continue;
    const x0 = xr.x[b - 1] + xr.w[b - 1];
    if (!(xr.x[b] > x0)) continue;
    gap(out, x0, top, xr.x[b] - x0, h, ink);
  }
  const yBlocks = s.section.y.blocks;
  for (let b = 1; b < yBlocks.length; b++) {
    if (yr.ve[b] <= yr.vs[b] || yr.ve[b - 1] <= yr.vs[b - 1]) continue;
    const y0 = yr.x[b - 1] + yr.w[b - 1];
    if (!(yr.x[b] > y0)) continue;
    gap(out, xr.x0, y0, xr.width, yr.x[b] - y0, ink);
  }
}

/** A gap's tint through the writer's `args` (nothing boxed): a rect without stroke. */
function gap(out: PassWriter, x: number, y: number, w: number, h: number, ink: number) {
  const g = out.args;
  g[0] = x;
  g[1] = y;
  g[2] = w;
  g[3] = h;
  g[4] = ink;
  g[5] = 0;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * A grid over blocks' chrome: a chain boundary is the blocks' gap with no line in it, a `seq-break` hairline only where
 * two blocks touch with no gap (down the plot at an x boundary, across it at a y one); the first track's label
 * at the plot's top in the label column; and each y block's label at the label column's right edge at the block's top
 * (below the track's label when it starts there), drawn when the block is at least 14 px tall. The x blocks' labels are
 * their section's, above.
 */
function drawBlocksChrome(s: SectionLayout, p: GridPlacement, out: PassWriter) {
  const { config, tokens } = s;
  const g = out.args;
  const blocks = s.section.blocks!;
  const yBlocks = gridYBlocks(s.section);
  const xr = p.xRow!;
  const yr = p.yRow!;
  const top = p.plotTop;
  const h = p.plotHeight;
  const ink = tokens.colors['seq-break'];
  const muted = tokens.colors['msui-text-muted'];
  const right = s.axis ? s.axis.x0 - config.labelGap / 2 : config.inset + config.labelWidth;
  // a boundary is the blocks' gap, with no line in it; a hairline a device pixel wide only where two blocks in view
  // touch with less than a device pixel between them (a boundary at the window's edge)
  const px = 1 / (s.axis?.dpr || 1);
  for (let b = 1; b < blocks.length; b++) {
    const at = touching(xr, b, px);
    if (at >= 0) out.fill(at - px / 2, top, px, h, ink);
  }
  for (let b = 1; b < yBlocks.length; b++) {
    const atY = touching(yr, b, px);
    if (atY >= 0) out.fill(xr.x0, atY - px / 2, xr.width, px, ink);
  }
  const t = s.tracks[0];
  if (t?.labelText) {
    g[0] = config.inset;
    g[1] = top + 12;
    g[2] = t.labelWidth;
    g[3] = config.labelWidth;
    g[4] = muted;
    out.textArgs(t.labelText, AnchorStart, config.labelFont);
  }
  for (let b = 0; b < yBlocks.length; b++) {
    const label = yBlocks[b].label;
    if (yr.ve[b] > yr.vs[b] && yr.w[b] >= 14) {
      g[0] = right;
      g[1] = yr.x[b] + (t && yr.x[b] - top < 20 ? 28 : 12);
      g[2] = 0;
      g[3] = config.labelWidth;
      g[4] = muted;
      out.textArgs(label, AnchorEnd, config.itemFont);
    }
  }
}

/** Where block `b` meets the block before it with no gap (less than `px` between them), −1 otherwise. */
function touching(r: BlockRow, b: number, px: number): number {
  if (r.ve[b] <= r.vs[b] || r.ve[b - 1] <= r.vs[b - 1]) return -1;
  const end = r.x[b - 1] + r.w[b - 1];
  return r.x[b] - end < px ? (end + r.x[b]) / 2 : -1;
}
