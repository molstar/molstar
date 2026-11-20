/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TextMetrics } from './text.js';
import { inkOn, type SeqTokens } from './tokens.js';

// Labels inside their marks (_spec/README.md "Render", Labels): the ranges layer's, placed at layout, drawn by the
// labels pass. The order: whole inside the first mark it fits; else beside in the track's label lane (the layer's);
// else inside cut with an ellipsis while at least `labels.minChars` letters are kept; else dropped.

/** A label's inset inside its mark, both sides (3 px a side), CSS px. */
export const LabelInset = 6;
/** None inside a mark narrower than this, where no text fits, CSS px. */
export const LabelMinPx = 8;
/**
 * The narrowest a glyph is, in ems of the item font (a lower bound: `i`, `l`, `.` are about 0.22 em in sans fonts): a
 * label of n glyphs is at least n of these wide, so one too long for a narrow mark is skipped unmeasured.
 */
const MinGlyphEm = 0.15;

/** The narrowest a glyph of `font` is (CSS px); 0 when its size is not in px (every label is then measured). */
export function minGlyphPx(font: string): number {
  const m = /(\d+(?:\.\d+)?)px/.exec(font);
  return m ? Number(m[1]) * MinGlyphEm : 0;
}

/** What a label inside reads of a ranges layer's layout: its marks, its items and its box. */
export interface LabelMarks {
  m: number;
  mx0: Float32Array;
  mx1: Float32Array;
  mitem: Int32Array;
  mtip: Uint8Array;
  /** Each mark's label: its width (0: none) and the letters kept of a cut one (0: whole). */
  mlabel: Float32Array;
  mcut: Int32Array;
  mink: Uint32Array;
  keys: (string | number)[];
  labels: (string | undefined)[];
  fills: number[];
  style: number;
  font: string;
  laneHeight: number;
  gap: number;
  box: { readonly x: number; readonly width: number };
  /** Scratch: the keys of the items labeled (kept, cleared each layout). */
  labeled: Set<string | number>;
}

/** The px a mark's arrowheads take on the sides `sides` names (1: its left, 2: its right, 3: both). */
export function tipRoom(lay: LabelMarks, j: number, sides: number): number {
  const tips = lay.mtip[j] & sides;
  if (!tips) return 0;
  const head = Math.min((lay.laneHeight - lay.gap) * 0.9, lay.mx1[j] - lay.mx0[j]);
  return tips === 3 ? 2 * head : head;
}

/** Where a mark's label starts: past its left arrowhead, within the box (a mark cut at the edge keeps its label in). */
export function labelX(lay: LabelMarks, j: number): number {
  return Math.max(lay.mx0[j] + tipRoom(lay, j, 1), lay.box.x) + LabelInset / 2;
}

/** The room a mark has for its label: its extent within the box, less its arrowheads (an arrowhead takes no text). */
function labelRoom(lay: LabelMarks, j: number): number {
  const x1 = Math.min(lay.mx1[j] - tipRoom(lay, j, 2), lay.box.x + lay.box.width);
  return x1 - labelX(lay, j) + LabelInset / 2;
}

/** Cut labels by label, the letters kept an index (`Spark…` once, whatever the layout): cleared when large. */
const Cuts = new Map<string, string[]>();

/** A label's first `n` letters and an ellipsis, made once. */
function cutLabel(label: string, n: number): string {
  let cuts = Cuts.get(label);
  if (!cuts) {
    if (Cuts.size > 4096) Cuts.clear();
    cuts = [];
    Cuts.set(label, cuts);
  }
  return (cuts[n] ??= `${label.slice(0, n)}…`);
}

/** The text a mark's label draws: its item's label whole, or cut. */
export function labelText(lay: LabelMarks, j: number): string {
  const label = lay.labels[lay.mitem[j]] as string;
  return lay.mcut[j] > 0 ? cutLabel(label, lay.mcut[j]) : label;
}

/**
 * Labels inside the marks of items not labeled yet, each item once (in the first of its marks, in reading order on the
 * page, its label fits): whole (`minChars` 0), or cut with an ellipsis keeping at least `minChars` letters (a label of
 * `minChars` letters or fewer only whole; the longest cut that fits). Allocation-free once warm: a cut is made once.
 */
export function labelsInside(lay: LabelMarks, text: TextMetrics, tokens: SeqTokens, minChars: number): void {
  const glyph = minGlyphPx(lay.font);
  const labeled = lay.labeled;
  for (let j = 0; j < lay.m; j++) {
    const i = lay.mitem[j];
    if (i < 0) continue;
    const label = lay.labels[i];
    if (!label || labeled.has(lay.keys[i])) continue;
    const w = labelRoom(lay, j) - LabelInset;
    if (w + LabelInset < LabelMinPx) continue;
    let m = 0;
    let cut = 0;
    if (!minChars) {
      // only whole: one surely too long for its mark is not measured (dense items)
      if (label.length * glyph > w) continue;
      m = text.measure(lay.font, label);
      if (m > w) continue;
    } else {
      if (label.length <= minChars || (minChars + 1) * glyph > w) continue;
      // the most letters kept that fit with the ellipsis, at least `minChars`
      let lo = minChars - 1;
      let hi = label.length - 1;
      while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if (text.measure(lay.font, cutLabel(label, mid)) <= w) lo = mid;
        else hi = mid - 1;
      }
      if (lo < minChars) continue;
      cut = lo;
      m = text.measure(lay.font, cutLabel(label, cut));
    }
    lay.mlabel[j] = m;
    lay.mcut[j] = cut;
    lay.mink[j] = lay.style ? tokens.colors['msui-text'] : inkOn(lay.fills[i], tokens);
    labeled.add(lay.keys[i]);
  }
}
