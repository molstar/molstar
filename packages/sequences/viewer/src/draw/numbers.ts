/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Numbers as text (seqspec:sources#levels-and-the-zoom-gate), part of the renderer API: one direct-mapped cache of
// their strings, bounded and shared by every context (a collision replaces the slot, so it never grows), and the
// wrapped layout's group labels above each line (`aboveEnd`), which the composer draws for a `labels` entry;
// `inlineStart` labels are the `letters` layer's. Allocation-free once warm.

import { GROUP_LABEL_PAD, groupLabelPx, groupOf } from '../runs.js';
import { SkinToken, TextAlign, TextBaseline, TextBatch } from './backend.js';
import type { DrawView } from './view.js';

const SIZE = 4096;
const keys = new Float64Array(SIZE).fill(Number.NaN);
const prefixes: string[] = new Array<string>(SIZE).fill('');
const texts: string[] = new Array<string>(SIZE).fill('');
const copies = new Int32Array(SIZE);
/**
 * `n` as text, `prefix` and a space before it when given ("GFP 12"), a copy past the first after it ("48 (copy 2)");
 * cached, so a warm frame makes no string.
 */
export function numberText(n: number, prefix = '', copy = 1): string {
  const h = (Math.imul((n | 0) ^ prefix.length ^ (copy << 12), 0x9e3779b1) >>> 20) & (SIZE - 1);
  if (keys[h] === n && prefixes[h] === prefix && copies[h] === copy) return texts[h];
  keys[h] = n;
  prefixes[h] = prefix;
  copies[h] = copy;
  const text = prefix ? `${prefix} ${n}` : String(n);
  texts[h] = copy > 1 ? `${text} (copy ${copy})` : text;
  return texts[h];
}

const labels = new TextBatch(64);
/**
 * Wrapped, `aboveEnd`: in the band above a line (`v` the `labels` entry, `onlyBlock` and the line set), each group's
 * last number (the run's labels, formatted when it was prepared), right-aligned to its last letter's right edge; a
 * label that would come within 2 px of the previous one, or past the data area's left edge, is left out.
 */
export function drawGroupLabels(v: DrawView): void {
  const b = v.onlyBlock;
  const run = b >= 0 ? v.runs[b] : undefined;
  if (!run || run.inline) return;
  const lo = v.lineStart;
  const hi = v.lineStop;
  const starts = run.startShown;
  const dpr = v.dpr;
  labels.reset();
  labels.style(v.text.font('mono', groupLabelPx(v.fontPx)), TextAlign.right, TextBaseline.alphabetic);
  labels.begin(SkinToken.inkQuiet);
  let clear = Math.round(v.left * dpr);
  for (let g = groupOf(run, lo); g < run.segments && starts[g] < hi; g++) {
    const last = (g + 1 < run.segments ? starts[g + 1] : run.count) - 1;
    if (last < lo || last >= hi || run.labelWidths[g] === 0) continue;
    const x = v.xe(b, last);
    if (x - run.labelWidths[g] * dpr < clear) continue;
    const k = labels.add(run.labels[g]);
    labels.x[k] = x;
    labels.y[k] = Math.round((v.top + v.height - 2) * dpr);
    clear = x + Math.round(2 * dpr);
  }
  if (labels.n > 0) v.backend.text(labels);
}

/**
 * Wrapped, `inlineStart`: each group's first number in shown [lo, hi) of block b, as a superscript 1 px before its
 * first letter (`x` places letters past it), in `out` (its style the caller's: right-aligned, middle).
 */
export function inlineGroupLabels(v: DrawView, b: number, lo: number, hi: number, out: TextBatch): void {
  const run = v.runs[b];
  if (!run?.inline) return;
  const starts = run.startShown;
  let g = groupOf(run, lo);
  if (starts[g] < lo) g++;
  const px = groupLabelPx(v.fontPx);
  // Caps are about 0.7 em: the centers differ by half the difference of the cap heights.
  const y = Math.round((v.top + v.height / 2 - 0.35 * (v.fontPx - px)) * v.dpr) + 0.5;
  out.begin(SkinToken.inkQuiet);
  for (; g < run.segments && starts[g] < hi; g++) {
    if (run.labelWidths[g] === 0) continue;
    const k = out.add(run.labels[g]);
    out.x[k] = v.xd(b, starts[g]) - Math.round(GROUP_LABEL_PAD * v.dpr);
    out.y[k] = y;
  }
}
