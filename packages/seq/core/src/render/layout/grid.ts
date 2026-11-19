/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The grid layout (_spec/README.md "Render", The grid layout): a section of two axes, one strip, the plot, whose
// transform is x and which carries a second transform `strip.y` down its frame (read by a layer as `ctx.axisY`). Its
// tracks' bands overlay the plot (later tracks over earlier, the last one's layers hit first). Its height is the data
// box's width over `aspect` (1: square cells), capped by the view's room below it less a tenth (never below
// `heights.grid`), at least two rows; `aspect: none` leaves it to the tracks' heights. Its viewport is the tracks' window plus `y` when y has a window of its own.

import { AxisTransform, type Viewport } from '../axis.js';
import { Anchor, type PassName, type PassWriter } from '../draw-lists.js';
import type { SectionLayout } from '../tracks.js';
import { StripPlacement } from './strips.js';
import { clampWindow, sameWindow, TracksNav } from './tracks.js';
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
    const top = ctx.box.y;
    const width = ctx.axis.width;
    const aspect = section.params.aspect;
    let h: number;
    if (aspect === 'none') {
      h = 0;
      for (const x of ctx.heights) h = Math.max(h, x);
    } else {
      const a = typeof aspect === 'number' && aspect > 0 ? aspect : 1;
      // the room below the section's top less a tenth, never below a grid track's default height: a view sized to its
      // content (whose height follows the plot's) keeps a stable plot
      const room = Math.max(config.heights.grid, (ctx.height - top) * RoomShare);
      h = Math.max(MinHeight, Math.min(width / a, room));
    }
    const i = strips.add(ctx.axis, 0, 0, 0, ctx.axis.u0);
    for (let j = 0; j < ctx.bands.length; j++) {
      strips.bandY[i * strips.bands + j] = top;
      strips.bandH[i * strips.bands + j] = h;
    }
    const y = ctx.y;
    p.yFirst = p.yLast = '';
    if (y?.pieces.length) {
      const ya = AxisTransform.on(y.axis, y.pieces, top, h, ctx.dpr, y.viewport);
      strips.ys[i] = ya;
      const first = ya.elementAtUnit(ya.u0);
      const last = ya.elementAtUnit(Math.max(ya.u0, ya.u1 - 1e-6));
      if (first >= 0) p.yFirst = String(first + 1);
      if (last >= 0) p.yLast = String(last + 1);
    }
    strips.top[i] = top;
    strips.height[i] = h;
    p.contentHeight = h;
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: h };
    ctx.layStrip(p, i);
    return p;
  },
  sameKey: sameGrid,
  nav: GridNav,
  chrome(s, pass, out, _emphasis, strip, band) {
    if (pass !== 'labels' || strip !== -1) return;
    drawChrome(s, pass, out);
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
  const top = strips.top[0];
  if (t) {
    g[0] = config.inset;
    g[1] = top + 12;
    g[2] = t.labelWidth;
    g[3] = config.labelWidth;
    g[4] = tokens.colors['msui-text-muted'];
    out.textArgs(t.track.label, AnchorStart, config.labelFont);
  }
  const p = s.placement;
  if (!(p instanceof GridPlacement)) return;
  const right = s.axis ? s.axis.x0 - config.labelGap / 2 : config.inset + config.labelWidth;
  for (let k = 0; k < 2; k++) {
    const text = k === 0 ? p.yFirst : p.yLast;
    if (!text) continue;
    g[0] = right;
    g[1] = k === 0 ? top + 28 : top + strips.height[0] - 4;
    g[2] = 0;
    g[3] = config.labelWidth;
    g[4] = tokens.colors['msui-text-muted'];
    out.textArgs(text, AnchorEnd, config.itemFont);
  }
}
