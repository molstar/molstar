/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text layout (_spec/README.md "Render", Text layout): a section's display units in lines of whole groups (a gap
// after every `group` units, bridged by spans), a strip a line in view translated down the page, its bands the
// section's tracks under a row of group labels, a band no layer used on a line collapsed there; numbers per line in a
// gutter left (`start`) or right (`end`); the lines laid out lazily (only those in view), their heights kept per
// section (full until laid out: the text's `end` is an estimate), the reader's place held by the anchor (`viewport.start`,
// the display unit at the top of the box, fractional within its line) across a reflow. A text section in a page caps at
// `section.maxPx` and scrolls inside its box; an export lays out every line (`figure`).

import { issueOf, type Value } from '@molstar/core/repr/tree';
import { DefaultNumbering, findNumbering, type NumberingOf } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedSection } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from '../axis.js';
import { fontSize, type RenderConfig } from '../config.js';
import { Anchor, type PassWriter } from '../draw-lists.js';
import { RenderIssues } from '../issues.js';
import { rulerBaseline } from '../ruler.js';
import type { TextMetrics } from '../text.js';
import type { SectionLayout } from '../tracks.js';
import type { Hit } from '../layer-type.js';
import { blockLength, chipWidth, isChip, type ResolvedBlock } from './blocks.js';
import { PlacePoint, StripPlacement, type Strips } from './strips.js';
import type { NavContext, PlaceContext, SectionLayoutType, SectionNav, StepDir } from './types.js';

// anchors read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;
const AnchorMiddle: number = Anchor.middle;
const AnchorEnd: number = Anchor.end;

/** A text section's params as the spec writes them (_spec/README.md "Resolved view"), with their defaults. */
export interface TextParams {
  readonly group: number;
  readonly numbers: 'start' | 'end' | 'both' | 'none';
  readonly groupLabels: 'above' | 'inline' | 'none';
  readonly align: 'start' | 'center';
  readonly maxWidth?: number;
  readonly height?: number;
}

export function textParams(params: Readonly<Record<string, Value>>): TextParams {
  const num = (k: string) => (typeof params[k] === 'number' ? (params[k] as number) : undefined);
  const str = <T extends string>(k: string, of: readonly T[], d: T): T =>
    of.includes(params[k] as T) ? (params[k] as T) : d;
  const maxWidth = num('max_width');
  const height = num('height');
  return {
    group: Math.max(0, Math.floor(num('group') ?? 10)),
    numbers: str('numbers', ['start', 'end', 'both', 'none'], 'start'),
    groupLabels: str('group_labels', ['above', 'inline', 'none'], 'above'),
    align: str('align', ['start', 'center'], 'start'),
    ...(maxWidth !== undefined ? { maxWidth } : {}),
    ...(height !== undefined ? { height } : {}),
  };
}

/** The least px a position a text shows: its letters' threshold (a text shows letters, EC-116). */
export function textMinScale(config: RenderConfig): number {
  return Math.max(config.density.letters, (config.letterSizes[0] ?? 0) * config.monoAdvance + 2);
}

/** A text section's line breaking: arithmetic on cell widths, gaps and measured labels. */
export interface TextGeometry {
  readonly total: number;
  /** px a position (the letter size) and the gap after a group (px; a group of 0: none). */
  readonly cell: number;
  readonly group: number;
  readonly gap: number;
  readonly unitsPerLine: number;
  /**
   * Each line's first display unit when lines do not all hold `unitsPerLine` (`text.pieceBreak: 'line'`: a piece edge
   * starts a line); undefined: line `l` starts at `l · unitsPerLine`.
   */
  readonly starts?: Float64Array;
  /** With `starts`: the pieces' first display units (what a new breaking starts lines at). */
  readonly pieces?: Float64Array;
  /** The width lines may take (px), and the part of a group's gap an inline group label takes (px). */
  readonly avail: number;
  readonly inline: number;
  /** Lines in the section, and those laid out at most (`text.maxLines`). */
  readonly lines: number;
  readonly laid: number;
  /** A line's x on the page and the width of a full one. */
  readonly lineX: number;
  readonly lineWidth: number;
  /** The row of group labels above the bands, the bands' full tops and heights, a full line's height and gap. */
  readonly labelRow: number;
  readonly bandTops: Float64Array;
  readonly bandHeights: readonly number[];
  /** Each band's padding inside its well, `[top, bottom]` a band (_spec/README.md "Render", Tracks layout). */
  readonly bandPads: Float64Array;
  readonly lineHeight: number;
  readonly lineGap: number;
  readonly numberFont: string;
  readonly params: TextParams;
  /** A section of blocks: each starts a line, its chips flowing like words (what a zoom breaks again by). */
  readonly blocks?: TextBlocks;
}

/**
 * A text's blocks (_spec/README.md "Render", Blocks): each block's first unit, and a chip's slot (px, its label's
 * width, at most a line; 0 for a block that is no chip), and the space between two chips (a share of a cell).
 */
export interface TextBlocks {
  readonly offsets: Float64Array;
  readonly slots: Float64Array;
  readonly gap: number;
  /**
   * A block's label inline on its first line (EC-95), px: its label in the label font and a word's gap, the cells it
   * takes from that line (a chip's is its pill: 0). Its continuation lines have none.
   */
  readonly lead: Float64Array;
}

/** The units block `b`'s inline label takes from its first line at `cell` px a unit (at most all but one). */
export function leadUnits(tb: TextBlocks, b: number, upl: number, cell: number): number {
  const px = tb.lead[b];
  return px > 0 && cell > 0 ? Math.max(0, Math.min(upl - 1, Math.ceil(px / cell))) : 0;
}

/** A text's blocks for its section: their offsets and chips' slots. */
export function textBlocks(
  blocks: readonly ResolvedBlock[],
  text: TextMetrics,
  config: RenderConfig,
  avail: number,
): TextBlocks {
  const offsets = new Float64Array(blocks.length);
  const slots = new Float64Array(blocks.length);
  const lead = new Float64Array(blocks.length);
  for (let b = 0; b < blocks.length; b++) {
    offsets[b] = blocks[b].offset ?? 0;
    const chip = isChip(blocks[b]);
    slots[b] = chip ? Math.min(avail, chipWidth(text, config, blocks[b].label)) : 0;
    // its label inline, with a word's gap after it (half the label font's size)
    if (!chip && blocks[b].label)
      lead[b] = text.measure(config.labelFont, blocks[b].label) + fontSize(config.labelFont) / 2;
  }
  return { offsets, slots, gap: config.text.groupGap, lead };
}

/**
 * Each line's first unit in a text of blocks: every block starts a line and is cut into lines of `upl`; a run of
 * chips flows like words, as many on a line as fit `avail` with `gap` between two.
 */
export function blockStarts(
  total: number,
  upl: number,
  avail: number,
  gap: number,
  tb: TextBlocks,
  cell = 0,
): Float64Array {
  const out: number[] = [];
  let line = -1;
  const n = tb.offsets.length;
  for (let b = 0; b < n; b++) {
    const o = tb.offsets[b];
    const end = b + 1 < n ? tb.offsets[b + 1] : total;
    if (tb.slots[b] > 0) {
      if (line >= 0 && line + gap + tb.slots[b] <= avail) line += gap + tb.slots[b];
      else {
        out.push(o);
        line = tb.slots[b];
      }
      continue;
    }
    line = -1;
    // its first line holds its inline label and fewer units (EC-95); the others are whole lines
    out.push(o);
    for (let u = o + upl - leadUnits(tb, b, upl, cell); u < end; u += upl) out.push(u);
  }
  return Float64Array.from(out);
}

/** The block holding section unit `u` (its index in `offsets`). */
function blockOf(tb: TextBlocks, u: number): number {
  let lo = 0;
  let hi = tb.offsets.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (tb.offsets[mid] <= u) lo = mid + 1;
    else hi = mid;
  }
  return Math.max(0, lo - 1);
}

/** The width of `n` units of a line: their cells and the gaps between their groups. */
export function lineWidthOf(g: Pick<TextGeometry, 'cell' | 'group' | 'gap'>, n: number): number {
  if (n <= 0) return 0;
  return n * g.cell + (g.group > 0 ? (Math.ceil(n / g.group) - 1) * g.gap : 0);
}

/** How many units a line of `avail` px holds: whole groups of `group` cells with `gap` after each (none: cells). */
export function unitsPerLineOf(avail: number, cell: number, group: number, gap: number): number {
  if (group > 0) return Math.max(1, Math.floor((avail + gap) / (group * cell + gap))) * group;
  return Math.max(1, Math.floor(avail / cell));
}

/** Line `l`'s first display unit. */
export function lineStart(g: Pick<TextGeometry, 'starts' | 'unitsPerLine'>, l: number): number {
  return g.starts ? g.starts[l] : l * g.unitsPerLine;
}

/** The display unit past line `l`'s last. */
export function lineEnd(g: Pick<TextGeometry, 'starts' | 'unitsPerLine' | 'total'>, l: number): number {
  if (!g.starts) return Math.min(g.total, (l + 1) * g.unitsPerLine);
  return l + 1 < g.starts.length ? g.starts[l + 1] : g.total;
}

/** The line holding display unit `u` (fractional: within its line). */
export function lineOf(g: Pick<TextGeometry, 'starts' | 'unitsPerLine'>, u: number): number {
  const s = g.starts;
  if (!s) return Math.max(0, Math.floor(u / g.unitsPerLine));
  let lo = 0;
  let hi = s.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (s[mid] <= u) lo = mid + 1;
    else hi = mid;
  }
  return Math.max(0, lo - 1);
}

/** A display unit as a line and the fraction of that line before it: what an anchor reads as (into `out`). */
function lineFrac(g: TextGeometry, u: number, out: { line: number; frac: number }): void {
  const l = Math.min(g.lines - 1, lineOf(g, u));
  const a = lineStart(g, l);
  out.line = l;
  out.frac = Math.max(0, Math.min(1, (u - a) / Math.max(1, lineEnd(g, l) - a)));
}

/**
 * Each line's first unit when a piece edge starts a line (`pieceBreak: 'line'`): each piece (`pieces`, each its first
 * display unit) cut into lines of `upl`; undefined when one piece, or none breaks.
 */
function pieceStarts(total: number, upl: number, pieces: ArrayLike<number> | undefined): Float64Array | undefined {
  if (!pieces || pieces.length < 2) return undefined;
  const out: number[] = [];
  for (let k = 0; k < pieces.length; k++) {
    const a = pieces[k];
    const b = k + 1 < pieces.length ? pieces[k + 1] : total;
    for (let u = a; u < b; u += upl) out.push(u);
  }
  return Float64Array.from(out);
}

/**
 * Breaks a section's `total` display units into lines for a width, a scale and its params; with `text.pieceBreak:
 * 'line'`, each of `pieces` (a piece's first display unit) starts a line.
 */
export function textGeometry(
  total: number,
  params: TextParams,
  width: number,
  scale: number | undefined,
  heights: readonly number[],
  config: RenderConfig,
  text: TextMetrics,
  pieces?: ArrayLike<number>,
  blocks?: readonly ResolvedBlock[],
  pads?: readonly number[],
): TextGeometry {
  const t = config.text;
  const cell = Math.max(textMinScale(config), scale ?? t.scale);
  const numberFont = config.labelFont;
  const numberW = text.measure(numberFont, String(Math.max(1, total))) + t.numberGap;
  const start = params.numbers === 'start' || params.numbers === 'both' ? numberW : 0;
  const end = params.numbers === 'end' || params.numbers === 'both' ? numberW : 0;
  const dataX = config.inset + config.labelWidth + config.labelGap + start;
  let avail = Math.max(cell, width - dataX - config.inset - end);
  if (params.maxWidth !== undefined) avail = Math.min(avail, params.maxWidth);
  const group = params.group;
  // an inline group label takes part in breaking: the gap holds it
  const inline = params.groupLabels === 'inline' && group > 0 ? numberW : 0;
  const gap = group > 0 ? t.groupGap * cell + inline : 0;
  const unitsPerLine = unitsPerLineOf(avail, cell, group, gap);
  const full = lineWidthOf({ cell, group, gap }, Math.min(unitsPerLine, Math.max(1, total)));
  const lineX = params.align === 'center' ? dataX + Math.max(0, (avail - full) / 2) : dataX;
  // a section of blocks: each starts a line, its chips flowing like words; else a piece edge may start one
  const tb = blocks ? textBlocks(blocks, text, config, avail) : undefined;
  const starts = tb
    ? blockStarts(total, unitsPerLine, avail, tb.gap * cell, tb, cell)
    : t.pieceBreak === 'line'
      ? pieceStarts(total, unitsPerLine, pieces)
      : undefined;
  const lines = starts ? Math.max(1, starts.length) : Math.max(1, Math.ceil(total / unitsPerLine));
  const labelRow = params.groupLabels === 'above' && group > 0 ? fontSize(numberFont) + 2 : 0;
  const bandTops = new Float64Array(heights.length);
  const bandPads = new Float64Array(2 * heights.length);
  if (pads) bandPads.set(pads.slice(0, bandPads.length));
  let y = labelRow;
  for (let j = 0; j < heights.length; j++) {
    if (j) y += t.bandGap;
    // a band inside its well: its padding above and below
    bandTops[j] = y + bandPads[2 * j];
    y += bandPads[2 * j] + heights[j] + bandPads[2 * j + 1];
  }
  return {
    total,
    cell,
    group,
    gap,
    unitsPerLine,
    ...(starts ? { starts, pieces: tb ? tb.offsets : Float64Array.from(pieces!) } : {}),
    ...(tb ? { blocks: tb } : {}),
    avail,
    inline,
    lines,
    laid: Math.min(lines, t.maxLines),
    lineX,
    lineWidth: full,
    labelRow,
    bandTops,
    bandHeights: heights,
    bandPads,
    lineHeight: y,
    lineGap: t.lineGap,
    numberFont,
    params,
  };
}

/** A line's bands at their full tops and heights, each inside its well (its padding). */
function setBands(strips: Strips, i: number, nb: number, g: TextGeometry): void {
  for (let j = 0; j < nb; j++) {
    const k = i * nb + j;
    strips.bandY[k] = g.bandTops[j];
    strips.bandH[k] = g.bandHeights[j];
    strips.padT[k] = g.bandPads[2 * j];
    strips.padB[k] = g.bandPads[2 * j + 1];
  }
}

/** A number, a group label or a track's label: what the text layout's chrome draws, in a strip's frame. */
interface TextMark {
  readonly strip: number;
  readonly text: string;
  readonly x: number;
  readonly y: number;
  readonly anchor: number;
  readonly width: number;
  readonly maxWidth: number;
  readonly font: string;
}

/** A text section placed: its lines' geometry, their heights (laid out, else full), the anchor, the chrome's marks. */
export class TextPlacement extends StripPlacement {
  geometry: TextGeometry | undefined;
  /** Each line's height with its gap after: measured once laid out, full before. */
  heights = new Float64Array(0);
  /** The geometry the heights were kept for: a reflow starts them again. */
  private heightsKey = '';
  /** The anchor (`viewport.start`, clamped), the px scrolled (`top`), and the window followers read. */
  anchor = 0;
  window: Viewport = { start: 0, end: 0 };
  marks: TextMark[] = [];
  /** Whether the section scrolls inside its box: its strips are clipped to it. */
  override clip = true;

  /** Keeps the line heights for a geometry; a new one fills them full. */
  keep(g: TextGeometry): void {
    const key = `${g.total}|${g.unitsPerLine}|${g.lines}|${g.lineHeight}|${g.lineGap}`;
    if (key === this.heightsKey && this.heights.length === g.laid) return;
    this.heightsKey = key;
    if (this.heights.length !== g.laid) this.heights = new Float64Array(g.laid);
    this.heights.fill(g.lineHeight + g.lineGap);
  }

  /** The px scrolled for an anchor through the placed lines (`Placement.scrollTopAt`). */
  scrollTopAt(anchor: number): number {
    return this.geometry ? topOf(this, anchor) : this.top;
  }

  /** The px from the first line's top to line `l`'s. */
  lineTop(l: number): number {
    let y = 0;
    const h = this.heights;
    for (let i = 0; i < l && i < h.length; i++) y += h[i];
    return y;
  }

  /** The line at px `y` from the first line's top, and the fraction of it above `y`: into `lineAt` and `frac`. */
  lineAt = 0;
  frac = 0;
  locate(y: number): void {
    const h = this.heights;
    let top = 0;
    let l = 0;
    while (l < h.length - 1 && top + h[l] <= y) top += h[l++];
    this.lineAt = l;
    this.frac = h.length ? Math.max(0, Math.min(1, (y - top) / h[l])) : 0;
  }
}

/** A line's transform: the section's pieces windowed to the line, its groups, at the line's x. */
function lineAxis(base: AxisTransform, g: TextGeometry, u0: number, u1: number): AxisTransform {
  return new AxisTransform(
    base.system,
    base.pieces,
    g.lineX,
    lineWidthOf(g, u1 - u0),
    base.dpr,
    { start: u0, end: u1 },
    {
      dir: base.dir,
      ...(base.length !== undefined ? { length: base.length } : {}),
      cyclic: base.cyclic,
      // groups counted from the line's first unit (a piece edge may start a line off a group's boundary)
      ...(g.group > 0 ? { groups: { size: g.group, gap: g.gap, origin: u0 } } : {}),
    },
  );
}

/** A block's line: its pieces windowed to its units `[b0, b1)`, its groups from the line's first, at the line's x. */
function blockLineAxis(
  block: ResolvedBlock,
  dpr: number,
  g: TextGeometry,
  b0: number,
  b1: number,
  dx = 0,
): AxisTransform {
  const ax = block.axis;
  return new AxisTransform(
    ax.system,
    ax.ranges ?? [],
    g.lineX + dx,
    lineWidthOf(g, b1 - b0),
    dpr,
    { start: b0, end: b1 },
    {
      dir: ax.dir,
      ...(ax.length !== undefined ? { length: ax.length } : {}),
      cyclic: ax.cyclic,
      ...(g.group > 0 ? { groups: { size: g.group, gap: g.gap, origin: b0 } } : {}),
    },
  );
}

/**
 * A line of chips from block `b` (before section unit `u1`): each chip a strip its slot wide over its positions,
 * the chips side by side with a word's gap between, the line's bands full (no band collapses on it).
 */
function chipLine(
  p: TextPlacement,
  ctx: PlaceContext,
  g: TextGeometry,
  tb: TextBlocks,
  blocks: readonly ResolvedBlock[],
  b: number,
  u1: number,
  y: number,
): void {
  const strips = p.strips;
  const nb = ctx.bands.length;
  let x = g.lineX;
  for (let k = b; k < blocks.length && tb.offsets[k] < u1; k++) {
    const block = blocks[k];
    if (!block.axis.ranges || !isChip(block)) continue;
    const a = AxisTransform.on(block.axis, block.axis.ranges, x, tb.slots[k], ctx.dpr, {
      start: 0,
      end: blockLength(block),
    });
    const i = strips.add(a, k, 0, y, tb.offsets[k]);
    strips.chip[i] = 1;
    setBands(strips, i, nb, g);
    strips.height[i] = g.lineHeight;
    ctx.layStrip(p, i);
    x += tb.slots[k] + tb.gap * g.cell;
  }
}

/**
 * The numbering a text's numbers print (_spec/README.md "Render", Text layout): the first of its layers that names one
 * (`numbering`, as a ruler does), on the section's system; `ordinal` when none does or the system has none of that id.
 */
export function textNumbering(section: ResolvedSection, model: StaticModel | undefined): NumberingOf | undefined {
  if (!model) return undefined;
  for (const t of section.tracks)
    for (const l of t.layers) {
      const id = l.params.numbering;
      if (typeof id !== 'string' || !id || id === DefaultNumbering) continue;
      const system = model.system(section.axis.system);
      return system ? findNumbering(system, id) : undefined;
    }
  return undefined;
}

/** A display unit's label in a numbering (its 1-based ordinal by default), as a number's text. */
function numberAt(a: AxisTransform, u: number, n: NumberingOf | undefined): string {
  const p = a.elementAtUnit(u);
  if (p < 0) return '';
  return n ? (n.numbering.labelAt(p - n.shift) ?? '') : String(p + 1);
}

export const TextLayout: SectionLayoutType<TextPlacement> = {
  name: 'text',
  fillBands: false,
  heights: (config) => config.text.heights,
  place(section, ctx, viewport, prev) {
    const p = prev ?? new TextPlacement();
    const { config } = ctx;
    const params = textParams(section.params);
    const g = textGeometry(
      ctx.axis.total,
      params,
      ctx.width,
      viewport?.scale,
      ctx.heights,
      config,
      ctx.text,
      ctx.axis.offsets,
      section.blocks,
      ctx.pads,
    );
    p.geometry = g;
    p.lineBlocks = g.blocks !== undefined;
    p.keep(g);
    if (g.lines > g.laid)
      ctx.report(
        issueOf(RenderIssues, 'seq.render.text.lines', `The text holds ${g.lines} lines`, {
          address: section.address,
          hint: RenderIssues.codes['seq.render.text.lines'].hint,
        }),
      );
    // the numbers' labels: a layer's numbering, else ordinals
    const numbering = textNumbering(section, ctx.model);
    const strips = p.strips;
    const nb = ctx.bands.length;
    strips.reset(nb);
    const marks: TextMark[] = [];
    p.marks = marks;
    const anchor = Math.max(0, Math.min(viewport?.start ?? 0, Math.max(0, g.total - 1e-9)));
    p.anchor = anchor;
    lineFrac(g, anchor, LineFrac);
    const first = Math.min(g.laid - 1, LineFrac.line);
    const frac = LineFrac.frac;
    const content = () => p.lineTop(g.laid);
    const top = ctx.box.y;
    // the box: the section's `height`, else its content's up to `section.maxPx`; an export shows every line
    const figure = ctx.figure === true;
    // its content's up to `section.maxPx` (a host sizing its canvas to the content reads it), or its `height`
    const boxOf = () => (figure ? content() : (params.height ?? Math.min(content(), config.section.maxPx)));
    let boxH = boxOf();
    const from = figure ? 0 : first;
    let y = figure ? top : top - frac * p.heights[first];
    const labelSize = fontSize(config.labelFont);
    for (let l = from; l < g.laid && (figure || y < top + boxH); l++) {
      const u0 = lineStart(g, l);
      const u1 = lineEnd(g, l);
      // a line of chips: each a strip its slot wide, flowing like words (_spec/README.md "Render", Blocks)
      const tb = g.blocks;
      const b = tb ? blockOf(tb, u0) : 0;
      if (tb && tb.slots[b] > 0) {
        chipLine(p, ctx, g, tb, section.blocks!, b, u1, y);
        p.heights[l] = g.lineHeight + g.lineGap;
        y += p.heights[l];
        continue;
      }
      // a block's first line: its label inline before its letters (EC-95)
      const lead = tb && u0 === tb.offsets[b] ? leadUnits(tb, b, g.unitsPerLine, g.cell) * g.cell : 0;
      const a = tb
        ? blockLineAxis(section.blocks![b], ctx.dpr, g, u0 - tb.offsets[b], u1 - tb.offsets[b], lead)
        : lineAxis(ctx.axis, g, u0, u1);
      // the line's units on its axis (a block's: from its first)
      const du = a.u0 - u0;
      const i = strips.add(a, b, 0, y, u0);
      const k0 = i * nb;
      setBands(strips, i, nb, g);
      ctx.layStrip(p, i);
      // a band no layer used on this line collapses (EC-96): the bands below move up and are laid out there
      let h = g.lineHeight;
      if (config.text.collapseEmpty) {
        let moved = false;
        let by = g.labelRow;
        let any = false;
        for (let j = 0; j < nb; j++) {
          const k = k0 + j;
          const keep = strips.used[k] > 0;
          if (keep && any) by += config.text.bandGap;
          // a kept band keeps its well (its padding); a collapsed one has none
          const pt = keep ? g.bandPads[2 * j] : 0;
          const pb = keep ? g.bandPads[2 * j + 1] : 0;
          if (strips.bandY[k] !== by + pt || (!keep && strips.bandH[k] > 0)) moved = true;
          strips.bandY[k] = by + pt;
          strips.bandH[k] = keep ? g.bandHeights[j] : 0;
          strips.padT[k] = pt;
          strips.padB[k] = pb;
          if (keep) {
            by += pt + g.bandHeights[j] + pb;
            any = true;
          }
        }
        h = by;
        if (moved) ctx.layStrip(p, i);
      }
      // the first line's bands carry their tracks' labels: each at its band's middle, pushed down to a label's own line
      // height below the one before when a thin band (a few px) would put them on top of each other
      if (l === 0) {
        let prev = Number.NEGATIVE_INFINITY;
        for (let j = 0; j < nb; j++) {
          if (strips.bandH[k0 + j] <= 0) continue;
          // centred on the band's well (its padding with it)
          const k = k0 + j;
          const mid = strips.bandY[k] + (strips.bandH[k] + strips.padB[k] - strips.padT[k]) / 2;
          const by = Math.max(mid + labelSize * 0.35, prev + labelSize);
          prev = by;
          const m = mark(i, ctx.bands[j].label, config.inset, by, AnchorStart, config.labelFont, ctx.text);
          marks.push({ ...m, maxWidth: config.labelWidth });
        }
      }
      strips.height[i] = h;
      p.heights[l] = h + g.lineGap;
      // the chrome's marks: the numbers, the group labels, the tracks' labels on the first line
      const k = strips.bandH[k0] > 0 ? k0 : k0 + firstShown(strips.bandH, k0, nb);
      // a ruler's band: the numbers on the ruler's label baseline (one row of numbers across the line), else at the middle
      const layers = ctx.bands[k - k0].layers;
      const ruled = layers.length > 0 && layers.every((l) => l.type === 'ruler');
      const baseline = ruled
        ? rulerBaseline(strips.bandY[k], strips.bandH[k])
        : strips.bandY[k] + strips.bandH[k] / 2 + labelSize * 0.35;
      if (lead > 0) {
        const m = mark(i, section.blocks![b].label, g.lineX, baseline, AnchorStart, config.labelFont, ctx.text);
        marks.push({ ...m, maxWidth: lead });
      }
      if (params.numbers === 'start' || params.numbers === 'both') {
        // in the gutter before the line, as an unlabeled line's: a first line's inline label is after it (EC-95)
        const s = numberAt(a, u0 + du, numbering);
        marks.push(mark(i, s, g.lineX - config.text.numberGap, baseline, AnchorEnd, g.numberFont, ctx.text));
      }
      if (params.numbers === 'end' || params.numbers === 'both') {
        const s = numberAt(a, u1 - 1 + du, numbering);
        marks.push(mark(i, s, a.x0 + a.width + config.text.numberGap, baseline, AnchorStart, g.numberFont, ctx.text));
      }
      // a number a group: above its last position, or inline in the gap after it (groups start at a line's start)
      if (g.group > 0 && params.groupLabels !== 'none')
        for (let e = u0 + g.group; e <= u1; e += g.group) {
          if (!a.cell(a.elementAtUnit(e - 1 + du))) continue;
          const s = numberAt(a, e - 1 + du, numbering);
          if (params.groupLabels === 'above')
            marks.push(mark(i, s, (a.sx0 + a.sx1) / 2, g.labelRow - 2, AnchorMiddle, config.itemFont, ctx.text));
          else if (e < u1) marks.push(mark(i, s, a.sx1 + g.gap / 2, baseline, AnchorMiddle, config.itemFont, ctx.text));
        }
      y += p.heights[l];
    }
    boxH = boxOf();
    p.contentHeight = content();
    p.top = p.lineTop(first) + frac * p.heights[first];
    p.scrollY = top;
    // the units in view (lines of blocks: from their first units)
    const n = strips.n;
    const last = n
      ? p.lineBlocks
        ? strips.first[n - 1] + strips.axes[n - 1].u1 - strips.axes[n - 1].u0
        : strips.axes[n - 1].u1
      : 0;
    p.window = { start: n ? (p.lineBlocks ? strips.first[0] : strips.axes[0].u0) : 0, end: last };
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: boxH };
    return p;
  },
  followWindow(section, axis, viewport, ctx) {
    const params = textParams(section.params);
    const g = textGeometry(
      axis.total,
      params,
      ctx.width,
      viewport?.scale,
      ctx.heights,
      ctx.config,
      ctx.text,
      axis.offsets,
      section.blocks,
      ctx.pads,
    );
    const pitch = g.lineHeight + g.lineGap;
    const boxH = params.height ?? Math.min(g.laid * pitch, ctx.config.section.maxPx);
    lineFrac(g, Math.max(0, Math.min(viewport?.start ?? 0, g.total)), LineFrac);
    const first = LineFrac.line;
    const lines = Math.max(1, Math.ceil(boxH / pitch + LineFrac.frac));
    return { start: lineStart(g, first), end: lineEnd(g, Math.min(g.lines - 1, first + lines - 1)) };
  },
  sameKey(a, b) {
    if (a === b) return true;
    if (!a || !b) return false;
    return a.start === b.start && (a.scale ?? 0) === (b.scale ?? 0) && (a.top ?? 0) === (b.top ?? 0);
  },
  get nav() {
    return TextNav;
  },
  hitChrome(s, x, y) {
    // a line's gutter (its band's label, its numbers) hits as the band's row (_spec/README.md "Render", Text layout)
    GutterAt[0] = x;
    GutterAt[1] = y;
    if (!s.placement.at(GutterAt, GutterPoint) || GutterPoint.kind !== 'gutter' || GutterPoint.band < 0)
      return undefined;
    const track = s.tracks[GutterPoint.band]?.track;
    // one record reused: a hover over the gutter makes none (read it, never keep it)
    const h = GutterHit as { -readonly [K in keyof Hit]: Hit[K] };
    h.layer = track?.address ?? s.section.address;
    h.section = s.id;
    h.system = s.axis?.system ?? s.section.axis.system;
    return h;
  },
  chrome(s, pass, out, _emphasis, strip, band) {
    // a pending section (no axis: an alignment's source not open) has the pending placement, no lines: nothing to draw
    if (pass !== 'labels' || band !== 0 || s.axis === undefined) return;
    if (!(s.placement instanceof TextPlacement)) throw new Error(`The text layout's chrome got another placement`);
    const marks = s.placement.marks;
    const g = out.args;
    const color = s.tokens.colors['msui-text-muted'];
    for (let i = 0; i < marks.length; i++) {
      const m = marks[i];
      if (m.strip !== strip) continue;
      g[0] = m.x;
      g[1] = m.y;
      g[2] = m.width;
      g[3] = m.maxWidth;
      g[4] = color;
      out.textArgs(m.text, m.anchor, m.font);
    }
  },
};

/** The point a gutter's hit reads, and a unit's line and fraction, reused. */
const GutterPoint = new PlacePoint();
/** The page point `GutterPoint` is placed from (`Placement.at` reads it from an array). */
const GutterAt = new Float64Array(2);
const GutterHit: Hit = { layer: '', section: '', system: '', target: 'row' };
const LineFrac = { line: 0, frac: 0 };

function firstShown(h: Float64Array, k0: number, nb: number): number {
  for (let j = 0; j < nb; j++) if (h[k0 + j] > 0) return j;
  return 0;
}

function mark(
  strip: number,
  text: string,
  x: number,
  y: number,
  anchor: number,
  font: string,
  metrics: TextMetrics,
): TextMark {
  return { strip, text, x, y, anchor, width: metrics.measure(font, text), maxWidth: Infinity, font };
}

/** The text placement a nav context carries, else undefined. */
function placed(c: NavContext): TextPlacement | undefined {
  return c.placement instanceof TextPlacement && c.placement.geometry ? c.placement : undefined;
}

/** The px from the first line's top to an anchor (a display unit, fractional within its line). */
function topOf(p: TextPlacement, anchor: number): number {
  if (!p.heights.length) return 0;
  lineFrac(p.geometry!, anchor, LineFrac);
  const line = Math.min(p.heights.length - 1, LineFrac.line);
  return p.lineTop(line) + LineFrac.frac * p.heights[line];
}

/** The anchor at px `y` from the first line's top. */
function anchorAt(p: TextPlacement, y: number): number {
  p.locate(y);
  const g = p.geometry!;
  const a = lineStart(g, p.lineAt);
  return a + p.frac * (lineEnd(g, p.lineAt) - a);
}

/**
 * The text layout's arithmetic (_spec/README.md "Render", Text layout): a pan scrolls lines (vertical), a zoom steps the
 * letter size keeping the line under the pointer at its y, a reveal brings a line not in view to the top, up and down
 * step a line (none past the first or the last: the cursor then crosses to the section above or below); the anchor is
 * clamped to the content and the scale to the letters' threshold and `maxScale`. Each reads the viewport it is given
 * (its anchor `start`), the placement only for the lines' breaking and heights.
 */
export const TextNav: SectionNav = {
  pan(c, v, _dx, dy) {
    const p = placed(c);
    if (!p || !dy) return v;
    const start = anchorAt(p, topOf(p, v.start) + dy);
    return TextNav.clamp(c, { ...v, start, end: start + (v.end - v.start) });
  },
  zoom(c, v, factor, at) {
    const p = placed(c);
    const scale = Math.max(c.minScale ?? 0, Math.min(c.maxScale, (v.scale ?? p?.geometry?.cell ?? 12) * factor));
    if (!p) return { ...v, scale };
    const g = p.geometry!;
    // the line holding `at` keeps its y: its distance in lines from the anchor, in the new breaking (exact: the
    // breaking's own arithmetic at the new cell)
    const cell = Math.max(c.minScale ?? 0, scale);
    const gap = g.group > 0 ? ((g.gap - g.inline) / g.cell) * cell + g.inline : 0;
    const upl1 = unitsPerLineOf(Math.max(cell, g.avail), cell, g.group, gap);
    lineFrac(g, v.start, LineFrac);
    const lines = lineOf(g, at) - LineFrac.line - LineFrac.frac;
    // the new breaking's lines (pieces breaking lines: their starts again at the new width)
    const starts = g.blocks
      ? blockStarts(g.total, upl1, g.avail, g.blocks.gap * cell, g.blocks, cell)
      : g.starts
        ? pieceStarts(g.total, upl1, g.pieces)
        : undefined;
    const next = { starts, unitsPerLine: upl1 };
    const l1 = lineOf(next, at) - lines;
    const whole = Math.floor(l1);
    const a = lineStart(next, Math.max(0, whole));
    const b = next.starts ? (whole + 1 < next.starts.length ? next.starts[whole + 1] : g.total) : a + upl1;
    const start = whole < 0 ? 0 : Math.max(0, a + (l1 - whole) * (b - a));
    return TextNav.clamp(c, { ...v, start, end: start + (v.end - v.start), scale });
  },
  fit(_c, v) {
    const { scale: _s, ...rest } = v;
    return { ...rest, start: 0, end: v.end - v.start };
  },
  reveal(c, v, start, end) {
    const p = placed(c);
    if (!p) return v;
    const g = p.geometry!;
    // in view: from the anchor's top to the box's bottom; the range's lines from its first line's top to its last's
    // bottom (no more than a box: a long range shows its first lines)
    const top = topOf(p, v.start);
    const first = Math.min(p.heights.length - 1, lineOf(g, start));
    const last = Math.min(p.heights.length - 1, lineOf(g, Math.max(start, end - 1)));
    // a line the box's bottom cuts, or past it, was laid in part or not at all (its unseen bands collapsed): read at its
    // whole height, so the least move never leaves it short of the box
    const bottom = top + p.box.height;
    let y0 = p.lineTop(first);
    let y1 = p.lineTop(last) + p.heights[last];
    if (y1 > bottom) {
      let y = 0;
      for (let l = 0; l <= last; l++) {
        if (l === first) y0 = y;
        const h = p.heights[l];
        y += y + h > bottom ? Math.max(h, g.lineHeight + g.lineGap) : h;
      }
      y1 = y;
    }
    y1 = Math.min(y1 - g.lineGap, y0 + p.box.height);
    if (y0 >= top && y1 <= top + p.box.height) return v;
    // above: its first line to the top; below: its last line's bottom to the box's bottom
    const to = y0 < top ? y0 : y1 - p.box.height;
    return TextNav.clamp(c, { ...v, start: anchorAt(p, Math.max(0, to)) });
  },
  focus(c, v, start) {
    const g = placed(c)?.geometry;
    return TextNav.clamp(c, { ...v, start: g ? lineStart(g, lineOf(g, Math.max(0, start))) : Math.max(0, start) });
  },
  step(c, offset, dir: StepDir) {
    const last = c.total - 1;
    if (dir === 'home') return 0;
    if (dir === 'end') return last;
    if (dir === 'left' || dir === 'right') {
      const next = offset + (dir === 'left' ? -1 : 1);
      return next >= 0 && next <= last ? next : -1;
    }
    const g = placed(c)?.geometry;
    if (!g) return -1;
    // the same column on the line above or below, its last unit when that line is shorter
    const l = lineOf(g, offset);
    const to = l + (dir === 'up' ? -1 : 1);
    if (to < 0 || to >= g.lines) return -1;
    const a = lineStart(g, to);
    return Math.min(a + (offset - lineStart(g, l)), lineEnd(g, to) - 1, last);
  },
  clamp(c, v) {
    const p = placed(c);
    let start = Math.max(0, v.start);
    if (p) {
      // the anchor no lower than where the last line's bottom meets the box's bottom
      const content = p.lineTop(p.geometry!.laid);
      start = Math.min(start, anchorAt(p, Math.max(0, content - p.box.height)));
    } else start = Math.min(start, Math.max(0, c.total - 1));
    const scale = v.scale === undefined ? undefined : Math.max(c.minScale ?? 0, Math.min(c.maxScale, v.scale));
    if (start === v.start && scale === v.scale) return v;
    return { ...v, start, ...(scale !== undefined ? { scale } : {}) };
  },
};

export type { PassWriter, ResolvedSection, SectionLayout, PlaceContext };
