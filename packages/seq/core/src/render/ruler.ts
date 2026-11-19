/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { findNumbering, type NumberingOf, parseLabel, type QuerySystem } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { fontSize } from './config.js';
import { Anchor, drawPass, grown, NoColor, type PassName, type PassWriter, passDraws } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { type Box, type Emphasis, type Hit, type LayerType, type LayoutContext, sectionId } from './layer-type.js';
import type { SeqTokens } from './tokens.js';

// the anchor a number read once (see ranges.ts)
const AnchorMiddle: number = Anchor.middle;

// The `ruler` layer type (_spec/README.md "Render", Layer types): the section's axis labelled by a numbering. Ticks at
// the labels that are multiples of a step of 1, 2 or 5 × 10ⁿ (the least whose px apart reach the label gap); insertion
// labels (100A, 111.1) are minor ticks, labeled where there is room; a label that does not fit is skipped, never shrunk.
// It labels the positions where they are, so a reversed axis counts down and a collapsed one labels its pieces.

/** The ruler's own position loop covers this many positions a piece at most; past it labels are looked up by number. */
const LoopLimit = 4096;
/** Insertion labels are minor ticks from this many px a position (a tick at each is a comb below it). */
const MinorScale = 4;
/** The labels' font: tabular figures (`mono`). */
const RulerFont = '10px mono';
const MajorHeight = 6;
const MinorHeight = 3;

export interface RulerPrep {
  /** The numbering id asked for (`ordinal` when absent). */
  readonly id: string;
  /** The system's numbering of that id; undefined: unknown (reported) and `ordinal` is used. */
  readonly numbering?: NumberingOf;
  readonly system?: QuerySystem;
}

export interface RulerLay {
  key: string;
  system: string;
  section: string;
  box: Box;
  axis: AxisTransform;
  tokens: SeqTokens;
  font: string;
  fontSize: number;
  /** The ticks in view: their center x, whether major, and their label (empty: none drawn) with its measured width. */
  n: number;
  x: Float32Array;
  major: Uint8Array;
  labels: string[];
  widths: Float32Array;
  /** Whether the section prints a label at each group's end (`group_labels`): the ruler leaves those positions to it. */
  grouped: '' | 'above' | 'inline';
  /** Labels drawn: indices into the ticks, with the room taken. */
  shown: Uint8Array;
}

export const RulerType: LayerType<RulerPrep, RulerLay> = {
  name: 'ruler',

  prepare(layer: ResolvedLayer, model: StaticModel): RulerPrep {
    const p = layer.params.numbering;
    const id = typeof p === 'string' && p ? p : 'ordinal';
    const system = layer.data ? model.system(layer.data.system) : undefined;
    const numbering = system ? findNumbering(system, id) : undefined;
    return { id, ...(numbering ? { numbering } : {}), ...(system ? { system } : {}) };
  },

  layout(prep, ctx: LayoutContext, _data, _overlay, prev?: RulerLay): RulerLay {
    const { axis, box, config, tokens } = ctx;
    const lay: RulerLay = prev ?? {
      key: ctx.key,
      system: axis.system,
      section: sectionId(ctx.section),
      box,
      axis,
      tokens,
      font: RulerFont,
      fontSize: fontSize(RulerFont),
      n: 0,
      x: new Float32Array(64),
      major: new Uint8Array(64),
      labels: [],
      widths: new Float32Array(64),
      shown: new Uint8Array(64),
      grouped: '',
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.font = RulerFont;
    lay.fontSize = fontSize(RulerFont);
    lay.n = 0;
    const gl = axis.gsize > 0 ? ctx.section.params.group_labels : undefined;
    lay.grouped = axis.gsize > 0 && gl !== 'none' ? (gl === 'inline' ? 'inline' : 'above') : '';
    if (!prep.numbering && prep.id !== 'ordinal' && prep.id !== 'index')
      ctx.report(
        issueOf(RenderIssues, 'seq.render.numbering', `The system has no numbering ${prep.id}`, {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.numbering'].hint,
        }),
      );
    const numbering = prep.numbering ?? (prep.system ? findNumbering(prep.system, 'ordinal') : undefined);
    const step = labelStep(axis.scale, config.ruler.labelGap);
    const ranges = axis.visible([]);
    for (let r = 0; r + 1 < ranges.length; r += 2) tickRange(lay, numbering, ranges[r], ranges[r + 1], step);
    fit(lay, ctx);
    return lay;
  },

  draw(lay, _pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(RulerDraws, lay, out, emphasis);
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return { layer: lay.key, section: lay.section, position, system: lay.system };
  },
};

/** The least step of 1, 2, 5 × 10ⁿ whose ticks are at least `gap` px apart at `scale` px a position. */
export function labelStep(scale: number, gap: number): number {
  const want = gap / Math.max(scale, 1e-9);
  const e = 10 ** Math.floor(Math.log10(want));
  for (const m of [1, 2, 5, 10]) if (m * e >= want) return m * e;
  return 10 * e;
}

/** The next tick's index, the arrays grown to hold it. */
function nextTick(lay: RulerLay): number {
  const j = lay.n++;
  if (j >= lay.x.length) {
    lay.x = grown(lay.x, j + 1);
    lay.major = grown(lay.major, j + 1);
    lay.widths = grown(lay.widths, j + 1);
    lay.shown = grown(lay.shown, j + 1);
  }
  return j;
}

/** A tick at position `p` with a label, when its cell is in view. */
function tickAt(lay: RulerLay, p: number, label: string, major: boolean) {
  const axis = lay.axis;
  if (!axis.cell(p)) return;
  const x = (axis.sx0 + axis.sx1) / 2;
  if (x < lay.box.x || x > lay.box.x + lay.box.width) return;
  const j = nextTick(lay);
  lay.x[j] = x;
  lay.major[j] = major ? 1 : 0;
  // a group's end number stands over its last cell: the position carries that label, the tick none (text layout, Numbers)
  lay.labels[j] = lay.grouped && axis.endsGroup(axis.unitAt(x), lay.grouped === 'inline') ? '' : label;
}

/** The ticks of the positions `[a, b)` in view (one piece): the labels that are multiples of the step, the insertions. */
function tickRange(lay: RulerLay, of: NumberingOf | undefined, a: number, b: number, step: number) {
  if (!of || b <= a) return;
  const { numbering, shift } = of;
  const builtIn = numbering.id === 'ordinal' || numbering.id === 'index';
  if (b - a > LoopLimit || builtIn) {
    // labels by number: every multiple of the step between the piece's first and last label
    const l0 = parseLabel(numbering.labelAt(a - shift) ?? '').int;
    const l1 = parseLabel(numbering.labelAt(b - 1 - shift) ?? '').int;
    if (l0 === undefined || l1 === undefined) return;
    const lo = Math.min(l0, l1);
    const hi = Math.max(l0, l1);
    if ((hi - lo) / step > 4 * LoopLimit) return;
    for (let n = Math.ceil(lo / step) * step; n <= hi; n += step) {
      const label = String(n);
      const i = builtIn ? n - (numbering.id === 'index' ? 0 : 1) : numbering.indexOf(label);
      const p = i + shift;
      if (i >= 0 && p >= a && p < b) tickAt(lay, p, label, true);
    }
    return;
  }
  const minors = lay.axis.scale >= MinorScale;
  for (let p = a; p < b; p++) {
    const label = numbering.labelAt(p - shift);
    if (label === undefined) continue;
    const parsed = parseLabel(label);
    if (parsed.int === undefined) tickAt(lay, p, label, true);
    else if (parsed.kind === 0) {
      if (parsed.int % step === 0) tickAt(lay, p, label, true);
    } else if (minors) tickAt(lay, p, label, false);
  }
}

/** Which labels are drawn: majors first, then minors, each only where it clears every one drawn and the data box. */
function fit(lay: RulerLay, ctx: LayoutContext) {
  const { box } = lay;
  const gap = Math.max(4, ctx.config.ruler.labelGap / 8);
  const taken: number[] = [];
  for (let pass = 1; pass >= 0; pass--)
    for (let j = 0; j < lay.n; j++) {
      lay.shown[j] = pass === 1 ? 0 : lay.shown[j];
      if (lay.major[j] !== pass || !lay.labels[j]) continue;
      const w = (lay.widths[j] = ctx.text.measure(lay.font, lay.labels[j]));
      const a = lay.x[j] - w / 2;
      const b = lay.x[j] + w / 2;
      if (a < box.x || b > box.x + box.width) continue;
      let free = true;
      for (let k = 0; k + 1 < taken.length && free; k += 2) free = b + gap <= taken[k] || a >= taken[k + 1] + gap;
      if (!free) continue;
      taken.push(a, b);
      lay.shown[j] = 1;
    }
}

/** Ticks on the box's bottom edge (toward the data) and labels above them. */
function drawTicks(lay: RulerLay, out: PassWriter) {
  const { box, tokens } = lay;
  const px = 1 / lay.axis.dpr;
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    const h = lay.major[j] ? MajorHeight : MinorHeight;
    g[0] = Math.round(lay.x[j] * lay.axis.dpr) / lay.axis.dpr - px / 2;
    g[1] = box.y + box.height - h;
    g[2] = px;
    g[3] = h;
    g[4] = tokens.colors['seq-ruler'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

function drawLabels(lay: RulerLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box } = lay;
  const g = out.args;
  for (let j = 0; j < lay.n; j++) {
    if (!lay.shown[j]) continue;
    g[0] = lay.x[j];
    g[1] = box.y + box.height - MajorHeight - 3;
    g[2] = lay.widths[j];
    g[3] = Number.POSITIVE_INFINITY;
    g[4] = lay.tokens.colors['seq-ruler'];
    out.textArgs(lay.labels[j], AnchorMiddle, lay.font);
  }
}

const RulerDraws = passDraws<RulerLay, Emphasis | undefined>({ data: drawTicks, labels: drawLabels });
