/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PassWriter } from './draw-lists.js';
import { Anchor } from './draw-lists.js';
import type { TextMetrics } from './text.js';

// The names of a rows system in the label column (_spec/README.md "Render", Row names): a row's name beside it, at the
// items' font, the text ink, cut with an ellipsis past its room, a name a row from the label minimum on. The rows
// layer's labels (an alignment's sequences) and the matrix layer's states (a VCF's samples, a methylation matrix's
// samples) are one part: a layer lays the names of its visible rows with `addName` in its layout and draws them with
// `drawNames` in the `labels` pass.

const AnchorStart: number = Anchor.start;

/** A name (already cut to its room by the caller) added to a layer's lists: its text, baseline `y` and measured width. */
export function addName(
  text: TextMetrics,
  font: string,
  label: string,
  y: number,
  texts: string[],
  ys: number[],
  ws: number[],
): void {
  texts.push(label);
  ys.push(y);
  ws.push(text.measure(font, label));
}

/** The names as text at `x` from each baseline, `max` px of room, in `ink` (the caller's clip bounds the band). */
export function drawNames(
  out: PassWriter,
  texts: readonly string[],
  ys: readonly number[],
  ws: readonly number[],
  font: string,
  ink: number,
  x: number,
  max: number,
): void {
  const g = out.args;
  for (let k = 0; k < texts.length; k++) {
    g[0] = x;
    g[1] = ys[k];
    g[2] = ws[k];
    g[3] = max;
    g[4] = ink;
    out.textArgs(texts[k], AnchorStart, font);
  }
}
