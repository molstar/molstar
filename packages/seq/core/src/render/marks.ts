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
  /** How far the hover's outline stands outside the cell (`emphasis.hoverOffset`). */
  readonly hoverOffset: number;
}

/**
 * The hover pass of a layer of cells: the outline of the hovered position's cell, outside it by `hoverOffset`, when the hover is on
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
  let a = Math.round(axis.sx0 * dpr) / dpr;
  let b = Math.round(axis.sx1 * dpr) / dpr;
  // a cell narrower than a device pixel (a far zoom) snaps to nothing: it takes the pixel it is in, so the hover never
  // blinks off as the pointer moves over cells that round to an empty span (_spec/README.md "Render", Emphasis)
  if (b - a < 1 / dpr) {
    a = Math.floor(((axis.sx0 + axis.sx1) / 2) * dpr) / dpr;
    b = a + 1 / dpr;
  }
  a = Math.max(a, box.x);
  b = Math.min(b, box.x + box.width);
  if (b <= a) return;
  // outside the cell by the offset, concentric: it never touches the selection's outline drawn on the cell
  const d = width / 2 + lay.hoverOffset;
  const g = out.args;
  g[0] = a - d;
  g[1] = box.y - d;
  g[2] = b - a + 2 * d;
  g[3] = box.height + 2 * d;
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
