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
// `text.maxHeight` and scrolls inside its box; an export lays out every line (`figure`).

import { issueOf, type Value } from '@molstar/core/repr/tree';
import type { ResolvedSection } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from '../axis.js';
import { fontSize, type RenderConfig } from '../config.js';
import { Anchor, type PassWriter } from '../draw-lists.js';
import { RenderIssues } from '../issues.js';
import type { TextMetrics } from '../text.js';
import type { SectionLayout } from '../tracks.js';
import { StripPlacement } from './strips.js';
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
  readonly lineHeight: number;
  readonly lineGap: number;
  readonly numberFont: string;
  readonly params: TextParams;
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

/** Breaks a section's `total` display units into lines for a width, a scale and its params. */
export function textGeometry(
  total: number,
  params: TextParams,
  width: number,
  scale: number | undefined,
  heights: readonly number[],
  config: RenderConfig,
  text: TextMetrics,
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
  const lines = Math.max(1, Math.ceil(total / unitsPerLine));
  const labelRow = params.groupLabels === 'above' && group > 0 ? fontSize(numberFont) + 2 : 0;
  const bandTops = new Float64Array(heights.length);
  let y = labelRow;
  for (let j = 0; j < heights.length; j++) {
    if (j) y += t.bandGap;
    bandTops[j] = y;
    y += heights[j];
  }
  return {
    total,
    cell,
    group,
    gap,
    unitsPerLine,
    avail,
    inline,
    lines,
    laid: Math.min(lines, t.maxLines),
    lineX,
    lineWidth: full,
    labelRow,
    bandTops,
    bandHeights: heights,
    lineHeight: y,
    lineGap: t.lineGap,
    numberFont,
    params,
  };
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
  top = 0;
  window: Viewport = { start: 0, end: 0 };
  marks: TextMark[] = [];
  /** Whether the section scrolls inside its box: its strips are clipped to it. */
  readonly clip = true;

  /** Keeps the line heights for a geometry; a new one fills them full. */
  keep(g: TextGeometry): void {
    const key = `${g.total}|${g.unitsPerLine}|${g.lineHeight}|${g.lineGap}`;
    if (key === this.heightsKey && this.heights.length === g.laid) return;
    this.heightsKey = key;
    if (this.heights.length !== g.laid) this.heights = new Float64Array(g.laid);
    this.heights.fill(g.lineHeight + g.lineGap);
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
      ...(g.group > 0 ? { groups: { size: g.group, gap: g.gap } } : {}),
    },
  );
}

/** A display unit's 1-based ordinal on the section's system, as a number's text. */
function numberAt(a: AxisTransform, u: number): string {
  const p = a.elementAtUnit(u);
  return p < 0 ? '' : String(p + 1);
}

export const TextLayout: SectionLayoutType<TextPlacement> = {
  name: 'text',
  fillBands: false,
  heights: (config) => config.text.heights,
  place(section, ctx, viewport, prev) {
    const p = prev ?? new TextPlacement();
    const { config } = ctx;
    const params = textParams(section.params);
    const g = textGeometry(ctx.axis.total, params, ctx.width, viewport?.scale, ctx.heights, config, ctx.text);
    p.geometry = g;
    p.keep(g);
    if (g.lines > g.laid)
      ctx.report(
        issueOf(RenderIssues, 'seq.render.text.lines', `The text holds ${g.lines} lines`, {
          address: section.address,
          hint: RenderIssues.codes['seq.render.text.lines'].hint,
        }),
      );
    const strips = p.strips;
    const nb = ctx.bands.length;
    strips.reset(nb);
    const marks: TextMark[] = [];
    p.marks = marks;
    const upl = g.unitsPerLine;
    const anchor = Math.max(0, Math.min(viewport?.start ?? 0, Math.max(0, g.total - 1e-9)));
    p.anchor = anchor;
    const first = Math.min(g.laid - 1, Math.floor(anchor / upl));
    const frac = (anchor - first * upl) / upl;
    const content = () => p.lineTop(g.laid);
    const top = ctx.box.y;
    // the box: the section's `height`, else its content's up to `text.maxHeight`; an export shows every line
    const figure = ctx.figure === true;
    // its content's up to `text.maxHeight` (a host sizing its canvas to the content reads it), or its `height`
    const boxOf = () => (figure ? content() : (params.height ?? Math.min(content(), config.text.maxHeight)));
    let boxH = boxOf();
    const from = figure ? 0 : first;
    let y = figure ? top : top - frac * p.heights[first];
    const labelSize = fontSize(config.labelFont);
    for (let l = from; l < g.laid && (figure || y < top + boxH); l++) {
      const u0 = l * upl;
      const u1 = Math.min(g.total, u0 + upl);
      const a = lineAxis(ctx.axis, g, u0, u1);
      const i = strips.add(a, 0, 0, y, u0);
      const k0 = i * nb;
      for (let j = 0; j < nb; j++) {
        strips.bandY[k0 + j] = g.bandTops[j];
        strips.bandH[k0 + j] = g.bandHeights[j];
      }
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
          if (strips.bandY[k] !== by || (!keep && strips.bandH[k] > 0)) moved = true;
          strips.bandY[k] = by;
          strips.bandH[k] = keep ? g.bandHeights[j] : 0;
          if (keep) {
            by += g.bandHeights[j];
            any = true;
          }
        }
        h = by;
        if (moved) ctx.layStrip(p, i);
      }
      strips.height[i] = h;
      p.heights[l] = h + g.lineGap;
      // the chrome's marks: the numbers, the group labels, the tracks' labels on the first line
      const k = strips.bandH[k0] > 0 ? k0 : k0 + firstShown(strips.bandH, k0, nb);
      const baseline = strips.bandY[k] + strips.bandH[k] / 2 + labelSize * 0.35;
      if (params.numbers === 'start' || params.numbers === 'both') {
        const s = numberAt(a, u0);
        marks.push(mark(i, s, a.x0 - config.text.numberGap, baseline, AnchorEnd, g.numberFont, ctx.text));
      }
      if (params.numbers === 'end' || params.numbers === 'both') {
        const s = numberAt(a, u1 - 1);
        marks.push(mark(i, s, a.x0 + a.width + config.text.numberGap, baseline, AnchorStart, g.numberFont, ctx.text));
      }
      // a number a group: above its last position, or inline in the gap after it (groups start at a line's start)
      if (g.group > 0 && params.groupLabels !== 'none')
        for (let e = u0 + g.group; e <= u1; e += g.group) {
          if (!a.cell(a.elementAtUnit(e - 1))) continue;
          const s = numberAt(a, e - 1);
          if (params.groupLabels === 'above')
            marks.push(mark(i, s, (a.sx0 + a.sx1) / 2, g.labelRow - 2, AnchorMiddle, config.itemFont, ctx.text));
          else if (e < u1) marks.push(mark(i, s, a.sx1 + g.gap / 2, baseline, AnchorMiddle, config.itemFont, ctx.text));
        }
      if (l === 0)
        for (let j = 0; j < nb; j++) {
          if (strips.bandH[k0 + j] <= 0) continue;
          const by = strips.bandY[k0 + j] + strips.bandH[k0 + j] / 2 + labelSize * 0.35;
          const m = mark(i, ctx.bands[j].label, config.inset, by, AnchorStart, config.labelFont, ctx.text);
          marks.push({ ...m, maxWidth: config.labelWidth });
        }
      y += p.heights[l];
    }
    boxH = boxOf();
    p.contentHeight = content();
    p.top = p.lineTop(first) + frac * p.heights[first];
    const last = strips.n ? strips.axes[strips.n - 1].u1 : 0;
    p.window = { start: strips.n ? strips.axes[0].u0 : 0, end: last };
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: boxH };
    return p;
  },
  followWindow(section, axis, viewport, ctx) {
    const params = textParams(section.params);
    const g = textGeometry(axis.total, params, ctx.width, viewport?.scale, ctx.heights, ctx.config, ctx.text);
    const pitch = g.lineHeight + g.lineGap;
    const boxH = params.height ?? Math.min(g.laid * pitch, ctx.config.text.maxHeight);
    const anchor = Math.max(0, Math.min(viewport?.start ?? 0, g.total));
    const first = Math.floor(anchor / g.unitsPerLine);
    const frac = anchor / g.unitsPerLine - first;
    const lines = Math.ceil(boxH / pitch + frac);
    const start = first * g.unitsPerLine;
    return { start, end: Math.min(g.total, start + Math.max(1, lines) * g.unitsPerLine) };
  },
  sameKey(a, b) {
    if (a === b) return true;
    if (!a || !b) return false;
    return a.start === b.start && (a.scale ?? 0) === (b.scale ?? 0) && (a.top ?? 0) === (b.top ?? 0);
  },
  get nav() {
    return TextNav;
  },
  chrome(s, pass, out, _emphasis, strip, band) {
    if (pass !== 'labels' || band !== 0) return;
    const p = s.placement as TextPlacement;
    const marks = p.marks;
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
  const g = p.geometry!;
  const line = Math.min(p.heights.length - 1, Math.floor(anchor / g.unitsPerLine));
  if (line < 0) return 0;
  return p.lineTop(line) + (anchor / g.unitsPerLine - line) * p.heights[line];
}

/** The anchor at px `y` from the first line's top. */
function anchorAt(p: TextPlacement, y: number): number {
  p.locate(y);
  return (p.lineAt + p.frac) * p.geometry!.unitsPerLine;
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
    const upl0 = g.unitsPerLine;
    const upl1 = unitsPerLineOf(Math.max(cell, g.avail), cell, g.group, gap);
    const lines = Math.floor(at / upl0) - v.start / upl0;
    const start = Math.max(0, Math.floor(at / upl1) * upl1 - lines * upl1);
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
    const first = Math.min(p.heights.length - 1, Math.floor(start / g.unitsPerLine));
    const last = Math.min(p.heights.length - 1, Math.floor(Math.max(start, end - 1) / g.unitsPerLine));
    const y0 = p.lineTop(first);
    const y1 = Math.min(p.lineTop(last) + p.heights[last] - g.lineGap, y0 + p.box.height);
    if (y0 >= top && y1 <= top + p.box.height) return v;
    // above: its first line to the top; below: its last line's bottom to the box's bottom
    const to = y0 < top ? y0 : y1 - p.box.height;
    return TextNav.clamp(c, { ...v, start: anchorAt(p, Math.max(0, to)) });
  },
  focus(c, v, start) {
    const p = placed(c);
    const upl = p?.geometry?.unitsPerLine ?? 1;
    return TextNav.clamp(c, { ...v, start: Math.floor(Math.max(0, start) / upl) * upl });
  },
  step(c, offset, dir: StepDir) {
    const last = c.total - 1;
    if (dir === 'home') return 0;
    if (dir === 'end') return last;
    const p = placed(c);
    const upl = p?.geometry?.unitsPerLine ?? 1;
    const d = dir === 'left' ? -1 : dir === 'right' ? 1 : dir === 'up' ? -upl : upl;
    const next = offset + d;
    if (next >= 0 && next <= last) return next;
    // down from the line above the last, shorter: the last unit
    if (dir === 'down' && Math.floor(offset / upl) < Math.floor(last / upl)) return last;
    return -1;
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
