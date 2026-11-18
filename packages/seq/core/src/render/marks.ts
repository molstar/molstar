/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AxisTransform } from './axis.js';
import { NoColor, type PassWriter } from './draw-lists.js';
import type { Box, Emphasis } from './layer-type.js';
import { type SeqTokens, withAlpha } from './tokens.js';

// Marks every layer type draws alike: spans as bands (placeholders, bars), and the hover's outline of a cell.

/** A rect per `x0, x1, cut` triple, over the band `y` to `y + h`, at least a device pixel wide. */
export function drawSpans(out: PassWriter, spans: readonly number[], y: number, h: number, fill: number, minW = 1) {
  for (let i = 0; i + 2 < spans.length; i += 3) out.fill(spans[i], y, Math.max(spans[i + 1] - spans[i], minW), h, fill);
}

/** What the hover's outline of a cell reads from a layer's layout. */
export interface CellLay {
  readonly axis: AxisTransform;
  readonly box: Box;
  readonly tokens: SeqTokens;
  readonly hoverWidth: number;
}

/**
 * The hover pass of a layer of cells: the outline of the hovered position's cell, inside the box, when the hover is on
 * the layer and the cell in view. It reads its numbers itself (the layout's, the hit's position): a number passed to a
 * call that is not inlined is boxed (nothing allocated once warm).
 */
export function hoverCell(lay: CellLay & { readonly key: string }, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || h.position === undefined) return;
  const { axis, box } = lay;
  const width = lay.hoverWidth;
  if (!axis.cell(h.position)) return;
  // `snap` inline, the rect through `args`: a hover allocates nothing once warm (this function is compiled while it
  // returns early, in the full draws, and a call left out of it then would box its numbers)
  const dpr = axis.dpr;
  const a = Math.max(Math.round(axis.sx0 * dpr) / dpr, box.x);
  const b = Math.min(Math.round(axis.sx1 * dpr) / dpr, box.x + box.width);
  if (b <= a) return;
  const half = width / 2;
  const g = out.args;
  g[0] = a + half;
  g[1] = box.y + half;
  g[2] = Math.max(b - a - width, 0.5);
  g[3] = box.height - width;
  g[4] = NoColor;
  g[5] = lay.tokens.colors['seq-hover'];
  g[6] = width;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** The failed placeholders, a conflict tint over each span: with none, nothing is drawn and nothing made. */
export function drawFailed(out: PassWriter, spans: readonly number[], y: number, h: number, tokens: SeqTokens) {
  if (spans.length) drawSpans(out, spans, y, h, withAlpha(tokens.colors['seq-conflict'], 0.25));
}
