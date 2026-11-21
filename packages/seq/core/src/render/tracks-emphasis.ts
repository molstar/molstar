/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location, LocationPart } from '@molstar/seqx-core/query';
import { axisBreaks, breakHolding } from '@molstar/seqx-core/view';
import { AxisTransform, Cut } from './axis.js';
import type { RenderConfig } from './config.js';
import { NoColor, type PassWriter, Round } from './draw-lists.js';
import type { Box, Emphasis, EmphasisHighlight, Hit, HitMark, LayerType } from './layer-type.js';
import { walkStripes } from './highlight-stripes.js';
import { partOn, rectOver } from './marks.js';
import type { SeqTokens } from './tokens.js';
import type { SectionLayout, ViewLayout } from './tracks-types.js';

// The view's emphasis marks (_spec/README.md "Render", Emphasis): the highlights' tints, bands and outlines in
// stripes, the selection's tint and outline, the hover carried from a linked section and the keyboard cursor, a box at
// a time on a strip's axis.

/** The tints under the data: the highlight channels', then the selection's. */
export function drawTints(
  out: PassWriter,
  s: SectionLayout,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  if (emphasis.highlights) drawHighlights(out, axis, box, emphasis, config);
  if (EmphasisPart !== Part.outlines)
    drawSelectionTint(out, axis, box, emphasis, tokens, config, projectedSides(s, axis, emphasis));
}

/**
 * Which of the emphasis a box draws: all of it, its tints alone (inside a track's well clip), or its outlines alone
 * (after the clip: an outline around a mark is never clipped).
 */
export const Part = { all: 0, tints: 1, outlines: 2 } as const;
let EmphasisPart: number = Part.all;

/** Sets which of the emphasis a box draws (`Part`): the view's draw, around a track's well clip. */
export function emphasisPart(part: number): void {
  EmphasisPart = part;
}

/** Whether a section draws the selection's projected parts: not a selected rect's sides where it shows both systems. */
function projectedSides(s: SectionLayout, axis: AxisTransform, emphasis: Emphasis): boolean {
  return !rectOver(emphasis.selection, axis.system, s.section.y?.axis.system);
}

/**
 * A section's emphasis marks of the writer's pass over `box` (a strip's band column, or one band's box) on the strip's
 * axis: the highlights' tints (over the data), the selection's tint and outline, the hover carried here from a linked
 * section (on the hover's own section only over its own-system letters tracks, and only for a letters hit lifted
 * through a codon map: a residue's codon on the DNA track beside it), the keyboard cursor.
 */
export function drawSectionEmphasis(
  out: PassWriter,
  s: SectionLayout,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  layout: ViewLayout,
  hovered = true,
) {
  const { tokens, config } = layout;
  const over = config.emphasis.selectionPlace === 'over';
  switch (out.pass.name) {
    case 'highlight':
      // over the data, below the selection's pass; under it, in the background pass
      if (over && emphasis.highlights) drawHighlights(out, axis, box, emphasis, config);
      break;
    case 'selection': {
      const sides = projectedSides(s, axis, emphasis);
      if (over && EmphasisPart !== Part.outlines) drawSelectionTint(out, axis, box, emphasis, tokens, config, sides);
      if (EmphasisPart !== Part.tints) drawSelectionOutline(out, axis, box, emphasis, tokens, config, sides);
      break;
    }
    case 'hover':
      if (emphasis.projected?.hover && (emphasis.projected.hoverSection !== s.section.address || !hovered))
        drawProjectedHover(out, axis, box, emphasis.projected.hover, tokens, config);
      break;
    case 'cursor':
      // a cursor on a cell of two axes is its layer's on its own section (it alone knows its y); the others mark its
      // column, unless it stands on a row's header (the column it left is no place of it)
      if (
        emphasis.cursor?.system === axis.system &&
        emphasis.cursor.part !== 'header' &&
        (emphasis.cursor.ySystem === undefined || emphasis.cursor.section !== s.section.address)
      )
        drawCursor(out, axis, box, emphasis.cursor, tokens, config);
      break;
  }
}

// The cursor and the emphasis marks take the layout, not its numbers: a number passed to a call that is not inlined
// is boxed (nothing allocated once warm). What a shared mark function needs besides its objects (a fill, a stroke, a
// width, the corners, an offset) goes through `Paint`, a typed array: its numbers are never boxed; colors at an
// opacity are made inline (`withAlpha` would box its result).

/**
 * The keyboard cursor: a thin outline of its cell over the box, in a neutral ink (`emphasis.cursorInk`); a squircle by
 * the cell over the box (every emphasis a squircle), inside it by half its width (concentric), square on a side the
 * view cuts and below the minimum.
 */
function drawCursor(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  cursor: { readonly position: number },
  tokens: SeqTokens,
  config: RenderConfig,
) {
  if (!axis.cell(cursor.position)) return;
  const w = config.outline.cursor;
  // `snap` inline, the rect through `args`: an emphasis redraw allocates nothing once warm
  const dpr = axis.dpr;
  let x0 = Math.round(axis.sx0 * dpr) / dpr;
  let x1 = Math.round(axis.sx1 * dpr) / dpr;
  // a side past the view's edge is cut there: square
  const round = (x0 >= axis.x0 ? Round.left : 0) | (x1 <= axis.x0 + axis.width ? Round.right : 0);
  x0 = Math.max(x0, axis.x0);
  x1 = Math.min(x1, axis.x0 + axis.width);
  const cw = x1 - x0;
  const h = box.height;
  // `squircleRadius` inline (a number returned from a call is boxed); the outline's less half its width
  const r =
    cw < config.squircleMin || h < config.squircleMin
      ? 0
      : Math.min(h * config.squircleRatio, config.squircleMax, cw / 2, h / 2);
  const g = out.args;
  g[0] = x0 + w / 2;
  g[1] = box.y + w / 2;
  g[2] = Math.max(cw - w, 0.5);
  g[3] = h - w;
  g[4] = NoColor;
  g[5] = tokens.colors[config.emphasis.cursorInk];
  g[6] = w;
  g[7] = r > 0 && round !== 0 ? Math.max(r - w / 2, 0) : 0;
  g[8] = round;
  out.rectArgs();
}

/**
 * What `paintRanges` draws with: the fill, the stroke and its width, whether the ends are squircles where they really
 * end (1) or square (0), how far the stroke stands outside the range (0: inside it, concentric with the fill), and
 * which of the range's ends are no real end (a highlight's run cut where the stripes change: 1 its start, 2 its end).
 */
const Paint = new Float64Array(6);
/** Whether a range wholly in a hidden stretch lights its break (a highlight's later runs of it do not: once). */
let LitBreaks = true;

/**
 * A part's ranges on a box: per range and piece of the axis in view, a rect over the box's height snapped to the device
 * grid, filled with `Paint[0]` (none: `NoColor`), then outlined with `Paint[1]` at `Paint[2]` px, inside it or
 * `Paint[4]` px outside it, concentric. With `Paint[3]`, an end the range really has is a squircle (a third of the
 * height, capped) and an end cut by the view or a piece is square (the frame rule).
 */
function paintRanges(out: PassWriter, axis: AxisTransform, box: Box, part: LocationPart, config: RenderConfig) {
  const fill = Paint[0];
  const stroke = Paint[1];
  const w = Paint[2];
  const rounded = Paint[3] !== 0;
  const outset = Paint[4];
  const dpr = axis.dpr;
  const g = out.args;
  const h = box.height;
  for (let i = 0; i + 1 < part.ranges.length; i += 2) {
    axis.piecesIn(part.ranges[i], part.ranges[i + 1]);
    const ka = axis.ka;
    const kb = axis.kb;
    let drew = false;
    for (let k = ka; k <= kb; k++) {
      if (!axis.spanIn(k, part.ranges[i], part.ranges[i + 1])) continue;
      drew = true;
      const x0 = Math.round(axis.sx0 * dpr) / dpr;
      const x1 = Math.max(Math.round(axis.sx1 * dpr) / dpr, x0 + 1 / dpr);
      // a run's ends that are no real end are cut (on the screen's sides: a reversed axis's start on the right)
      const open = Paint[5];
      let cut = axis.scut;
      if ((open & 1) !== 0) cut |= axis.dir === -1 ? Cut.right : Cut.left;
      if ((open & 2) !== 0) cut |= axis.dir === -1 ? Cut.left : Cut.right;
      const bw = x1 - x0;
      // `squircleRadius` inline (a number returned from a call is boxed)
      const r =
        !rounded || bw < config.squircleMin || h < config.squircleMin
          ? 0
          : Math.min(h * config.squircleRatio, config.squircleMax, bw / 2, h / 2);
      const round = r > 0 ? ((cut & Cut.left) === 0 ? Round.left : 0) | ((cut & Cut.right) === 0 ? Round.right : 0) : 0;
      if (fill !== NoColor) {
        g[0] = x0;
        g[1] = box.y;
        g[2] = bw;
        g[3] = h;
        g[4] = fill;
        g[5] = NoColor;
        g[6] = 0;
        g[7] = round ? r : 0;
        g[8] = round;
        out.rectArgs();
      }
      if (stroke !== NoColor) {
        // inside: inset by half the width; outside: out by the offset and half the width; the radius follows
        const d = outset > 0 ? -(outset + w / 2) : w / 2;
        g[0] = x0 + d;
        g[1] = box.y + d;
        g[2] = Math.max(bw - 2 * d, 0.5);
        g[3] = Math.max(h - 2 * d, 0.5);
        g[4] = NoColor;
        g[5] = stroke;
        g[6] = w;
        g[7] = round ? Math.max(r - d, 0) : 0;
        g[8] = round;
        out.rectArgs();
      }
    }
    // a range wholly in a hidden stretch lights its break: what a break hides is not dropped
    if (!drew && LitBreaks) litBreak(out, axis, box, part.ranges[i]);
  }
}

/**
 * The break hiding position `p`, when in view: a bar two device pixels wide across the box in the paint's fill, else its
 * stroke (read from `Paint`: a color passed to a call that is not inlined is boxed).
 */
function litBreak(out: PassWriter, axis: AxisTransform, box: Box, p: number) {
  const b = axisBreaks(axis);
  const i = b.length ? breakHolding(axis, p) : -1;
  if (i < 0 || b[i] < axis.u0 || b[i] > axis.u1) return;
  const px = 1 / axis.dpr;
  const g = out.args;
  g[0] = Math.round(axis.xOfUnit(b[i]) * axis.dpr) / axis.dpr - px;
  g[1] = box.y;
  g[2] = 2 * px;
  g[3] = box.height;
  g[4] = Paint[0] !== NoColor ? Paint[0] : Paint[1];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The selection's accent tint over a box's positions: its own part at `emphasis.selectionTint`, else the part
 * projected here at that times `emphasis.projected` (the same look, lighter). Squircle ends where it really ends.
 */
function drawSelectionTint(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
  projected: boolean,
) {
  const own = partOn(emphasis.selection, axis.system);
  const part = own ?? (projected ? partOn(emphasis.projected?.selection, axis.system) : undefined);
  if (!part) return;
  const accent = tokens.colors['msui-accent'];
  const a = config.emphasis.selectionTint * (own ? 1 : config.emphasis.projected);
  Paint[0] = ((accent & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
  Paint[1] = NoColor;
  Paint[2] = 0;
  Paint[3] = 1;
  Paint[4] = 0;
  Paint[5] = 0;
  paintRanges(out, axis, boxed(box), part, config);
}

/**
 * The selection's one thin outline around each range on a box, inside its tint: the accent at
 * `emphasis.selectionOutline`, its own part, else the projected one, lighter.
 */
function drawSelectionOutline(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
  projected: boolean,
) {
  const own = partOn(emphasis.selection, axis.system);
  const part = own ?? (projected ? partOn(emphasis.projected?.selection, axis.system) : undefined);
  if (!part) return;
  const accent = tokens.colors['msui-accent'];
  const a = ((accent & 0xff) / 255) * config.emphasis.selectionOutline * (own ? 1 : config.emphasis.projected);
  Paint[0] = NoColor;
  Paint[1] = ((accent & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
  Paint[2] = config.outline.selection;
  Paint[3] = 1;
  Paint[4] = 0;
  Paint[5] = 0;
  paintRanges(out, axis, boxed(box), part, config);
}

/**
 * The hover carried onto a box from a linked section: the hover's outline around each range, `emphasis.hoverOffset`
 * outside it (so it never touches the selection's), lighter.
 */
function drawProjectedHover(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  projected: Location,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const part = partOn(projected, axis.system);
  if (!part) return;
  const hover = tokens.colors['seq-hover'];
  const a = ((hover & 0xff) / 255) * config.emphasis.projected;
  Paint[0] = NoColor;
  Paint[1] = ((hover & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
  Paint[2] = config.outline.hover;
  // every emphasis a squircle (_spec/README.md "Render", Shapes): round where it really ends
  Paint[3] = 1;
  Paint[4] = config.emphasis.hoverOffset;
  Paint[5] = 0;
  paintRanges(out, axis, boxed(box), part, config);
}

/**
 * The mark an extension's type writes its hovered item's bounds into: an instance of its own class (an object literal
 * of this shape shares its map with every other such literal, whose fields may hold other values: a number written into
 * one is then boxed, a warm hover allocating).
 */
class MarkBox implements HitMark {
  x = 0.5;
  y = 0.5;
  width = 0.5;
  height = 0.5;
  radius = 0.5;
}
const HoverMark: HitMark = new MarkBox();

/**
 * The hover of an item of an extension's layer type (_spec/README.md "Render", Emphasis): the outline around the mark
 * the type reports (`LayerType.mark`), `emphasis.hoverOffset` outside it and concentric (its radius plus the offset),
 * in the hover ink, as a built-in's hover looks. Round on the sides the view does not cut, square where the mark meets
 * the strip's edge. Nothing allocated: the mark is written into one reused box.
 */
export function drawExtensionHover(
  out: PassWriter,
  type: LayerType<any, any>,
  lay: unknown,
  hover: Hit,
  axis: AxisTransform,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const key = hover.element;
  if (!type.mark || !key) return;
  const m = HoverMark;
  m.radius = 0;
  if (!type.mark(lay, key, m)) return;
  return drawMarkOutline(out, m, axis, tokens, config);
}

/** An outline around a mark's bounds (`drawExtensionHover`). */
export function drawMarkOutline(
  out: PassWriter,
  m: HitMark,
  axis: AxisTransform,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const half = config.outline.hover / 2 + config.emphasis.hoverOffset;
  const r = m.radius ?? 0;
  const g = out.args;
  g[0] = m.x - half;
  g[1] = m.y - half;
  g[2] = m.width + 2 * half;
  g[3] = m.height + 2 * half;
  g[4] = NoColor;
  g[5] = tokens.colors['seq-hover'];
  g[6] = config.outline.hover;
  g[7] = r > 0 ? r + half : 0;
  // a side the strip cuts is square
  g[8] = (m.x > axis.x0 ? Round.left : 0) | (m.x + m.width < axis.x0 + axis.width ? Round.right : 0);
  out.rectArgs();
}

/** A mutable box: one stripe of a band (reused: a stripe makes no object). */
class StripeBox implements Box {
  x = 0;
  y = 0;
  width = 0;
  height = 0;
}
const Stripe = new StripeBox();
/**
 * A box copied into one of `paintRanges`' own (a strip's box, a stripe): it reads one class of box, so its loads stay
 * monomorphic and a warm redraw boxes no number (a load merged over boxes of several shapes did, 16 bytes a call).
 */
const PaintBox = new StripeBox();
function boxed(box: Box): Box {
  PaintBox.x = box.x;
  PaintBox.y = box.y;
  PaintBox.width = box.width;
  PaintBox.height = box.height;
  return PaintBox;
}
/** The one `[start, end)` range of a highlight run, as a part (reused). */
const RunRanges = new Int32Array(2);
const RunPart: LocationPart = { system: '', ranges: RunRanges };
/** The channels with a part on this axis that take a stripe (fills), and those that do not (outlines): by index. */
const HParts: (LocationPart | undefined)[] = [];
const HIndex: number[] = [];
const HOutlines: number[] = [];
/** What `emitRun` draws on (set before `walkStripes`: it is given no closure). */
let HOut!: PassWriter;
let HAxis!: AxisTransform;
let HBox!: Box;
let HList!: readonly EmphasisHighlight[];
let HConfig!: RenderConfig;

/**
 * The highlight channels on a box, in order (_spec/README.md "Render", Emphasis): where one channel alone lights a
 * position it takes the whole band; where several do, the band is split into equal stripes, one a channel present, in
 * channel order, at most `emphasis.highlightStripes` (the rest share the last one). A channel's `style` is drawn in its
 * stripe: `tint` (its color at `emphasis.highlightTint`) or `band` (a bar along the stripe's foot); an `outline` takes
 * no stripe (it has no fill): one outline over the channel's own ranges on the whole band, above the fills.
 */
function drawHighlights(out: PassWriter, axis: AxisTransform, box: Box, emphasis: Emphasis, config: RenderConfig) {
  const hs = emphasis.highlights;
  if (!hs) return;
  let n = 0;
  let o = 0;
  for (let i = 0; i < hs.length; i++) {
    const part = partOn(hs[i].location, axis.system);
    if (!part) continue;
    if (hs[i].style === 'outline') HOutlines[o++] = i;
    else {
      HParts[n] = part;
      HIndex[n] = i;
      n++;
    }
  }
  // the fills (tints, bands) and the outlines drawn apart when a well's clip holds the fills alone
  const fills = EmphasisPart === Part.outlines ? 0 : n;
  if (EmphasisPart === Part.tints) o = 0;
  if (fills === 1) paintHighlight(out, axis, box, hs[HIndex[0]], HParts[0]!, 0, 1, config);
  else if (fills > 1) {
    HOut = out;
    HAxis = axis;
    HBox = box;
    HList = hs;
    HConfig = config;
    walkStripes(HParts, n, Math.max(1, config.emphasis.highlightStripes | 0), emitRun);
  }
  for (let k = 0; k < o; k++) {
    const h = hs[HOutlines[k]];
    paintHighlight(out, axis, box, h, partOn(h.location, axis.system)!, 0, 1, config);
  }
  for (let c = 0; c < n; c++) HParts[c] = undefined;
}

/** The ends of the run `emitRun` paints that are no real ends (`Paint[5]`); 0 for a channel's whole ranges. */
let RunOpen = 0;

/** `walkStripes`' run: a channel's `[start, end)` in its stripe. A hidden range lights its break once, at its first run. */
function emitRun(
  channel: number,
  start: number,
  end: number,
  slot: number,
  count: number,
  first: number,
  last: number,
) {
  RunRanges[0] = start;
  RunRanges[1] = end;
  LitBreaks = first === 1;
  // a run's ends where the stripes change are no real ends: square there
  RunOpen = (first === 1 ? 0 : 1) | (last === 1 ? 0 : 2);
  paintHighlight(HOut, HAxis, HBox, HList[HIndex[channel]], RunPart, slot, count, HConfig);
  RunOpen = 0;
  LitBreaks = true;
}

/** `v` snapped to the device pixels of an axis. */
const snap = (v: number, dpr: number) => Math.round(v * dpr) / dpr;

/** One channel's part in stripe `slot` of `stripes` of the box, in its style; stripe edges on device pixels. */
function paintHighlight(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  h: EmphasisHighlight,
  part: LocationPart,
  slot: number,
  stripes: number,
  config: RenderConfig,
) {
  const em = config.emphasis;
  const dpr = axis.dpr;
  Stripe.x = box.x;
  Stripe.width = box.width;
  Stripe.y = box.y;
  Stripe.height = box.height;
  if (stripes > 1) {
    // one edge between two stripes, so they neither overlap nor leave a seam
    const y0 = snap(box.y + (slot * box.height) / stripes, dpr);
    const y1 = snap(box.y + ((slot + 1) * box.height) / stripes, dpr);
    if (y1 <= y0) return;
    Stripe.y = y0;
    Stripe.height = y1 - y0;
  }
  const c = h.color;
  const style = h.style ?? 'tint';
  if (style === 'band') {
    const bottom = Stripe.y + Stripe.height;
    const bh = Math.min(Math.max(snap(em.highlightBand, dpr), 1 / dpr), Stripe.height);
    Stripe.y = bottom - bh;
    Stripe.height = bh;
    Paint[0] = c;
    Paint[1] = NoColor;
    Paint[2] = 0;
  } else if (style === 'outline') {
    Paint[0] = NoColor;
    Paint[1] = c;
    Paint[2] = config.outline.selection;
  } else {
    const a = ((c & 0xff) / 255) * em.highlightTint;
    Paint[0] = ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
    Paint[1] = NoColor;
    Paint[2] = 0;
  }
  // every style a squircle by its stripe's height (a `band` bar below the minimum: square), round where it really ends
  Paint[3] = 1;
  Paint[4] = 0;
  Paint[5] = RunOpen;
  paintRanges(out, axis, Stripe, part, config);
}
