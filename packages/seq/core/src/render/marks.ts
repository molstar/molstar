/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rangesLocation } from '@molstar/seqx-core/query';
import type { ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import type { RenderConfig } from './config.js';
import { NoColor, type PassWriter, Round } from './draw-lists.js';
import type { Box, Emphasis, Hit } from './layer-type.js';
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
  /** The squircle rule the outline takes (`squircleRatio`, `squircleMax`, `squircleMin`). */
  readonly config: RenderConfig;
}

/**
 * The hover pass of a layer of cells: the outline of the hovered position's cell, outside it by `hoverOffset`, when the hover is on
 * the layer and the cell in view; a squircle by the cell's size (_spec/README.md "Render", Shapes: every emphasis a
 * squircle), square below `squircleMin` and on a side the view cuts, its radius the cell's plus the offset
 * (concentric). It reads its numbers itself (the layout's, the hit's position): a number passed to a call that is not
 * inlined is boxed (nothing allocated once warm).
 */
export function hoverCell(lay: CellLay & { readonly key: string }, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || h.position === undefined) return;
  const { axis, box } = lay;
  const width = lay.hoverWidth;
  // a block's strip has its own system: the position is another system's on any other
  if (h.system !== axis.system || !axis.cell(h.position)) return;
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
  // a side the view's edge cuts is square (the frame rule)
  const round = (a >= box.x ? Round.left : 0) | (b <= box.x + box.width ? Round.right : 0);
  a = Math.max(a, box.x);
  b = Math.min(b, box.x + box.width);
  if (b <= a) return;
  // outside the cell by the offset, concentric: it never touches the selection's outline drawn on the cell
  const d = width / 2 + lay.hoverOffset;
  // the cell's squircle (`squircleRadius` inline: a number returned from a call is boxed), the outline's its plus `d`
  const c = lay.config;
  const cw = b - a;
  const bh = box.height;
  const r =
    cw < c.squircleMin || bh < c.squircleMin ? 0 : Math.min(bh * c.squircleRatio, c.squircleMax, cw / 2, bh / 2);
  const g = out.args;
  g[0] = a - d;
  g[1] = box.y - d;
  g[2] = cw + 2 * d;
  g[3] = bh + 2 * d;
  g[4] = NoColor;
  g[5] = lay.tokens.colors['seq-hover'];
  g[6] = width;
  g[7] = r > 0 && round !== 0 ? r + d : 0;
  g[8] = round;
  out.rectArgs();
}

/** The failed placeholders, a conflict tint over each span: with none, nothing is drawn and nothing made. */
export function drawFailed(out: PassWriter, spans: readonly number[], y: number, h: number, tokens: SeqTokens) {
  if (spans.length) drawSpans(out, spans, y, h, withAlpha(tokens.colors['seq-conflict'], 0.25));
}

/** The share of a block's height an arrow's body takes (its head takes all of it). */
export const ArrowBody = 0.6;

const arrowPoly = new Float64Array(14);

/**
 * A block with a pointed end, as one path (the tail square): the box `args[0..3]` holds (x, y, width, height), its tip
 * on the right (`args[8]` > 0) or the left, its head `args[7]` long (at most the width), fill, stroke and stroke width
 * `args[4..6]`. The body is `ArrowBody` of the height, centered; the head is all of it.
 */
export function arrowArgs(out: PassWriter) {
  const a = out.args;
  const x = a[0];
  const y = a[1];
  const w = a[2];
  const h = a[3];
  const head = Math.min(a[7], w);
  const inset = (h * (1 - ArrowBody)) / 2;
  const right = a[8] > 0;
  const p = arrowPoly;
  // the tail, then along the body to the head's base, out to its corner, the tip, back
  const tail = right ? x : x + w;
  const base = right ? x + w - head : x + head;
  const tip = right ? x + w : x;
  p[0] = tail;
  p[1] = y + inset;
  p[2] = base;
  p[3] = y + inset;
  p[4] = base;
  p[5] = y;
  p[6] = tip;
  p[7] = y + h / 2;
  p[8] = base;
  p[9] = y + h;
  p[10] = base;
  p[11] = y + h - inset;
  p[12] = tail;
  p[13] = y + h - inset;
  a[7] = 0;
  out.pathArgs(p, 7, true);
}

/**
 * A break as a hit (`target: 'break'`): the section's position at its edge (the last shown before it in display order),
 * the hidden stretch as its extent on the section's system (two ranges across a cyclic origin).
 */
export function breakHit(
  s: { readonly section: Pick<ResolvedSection, 'address'>; readonly id: string },
  axis: AxisTransform,
  x: number,
  hiddenStart: number,
  hiddenEnd: number,
  count?: number,
): Hit {
  const length = axis.length ?? hiddenStart;
  const ranges = hiddenEnd >= hiddenStart ? [hiddenStart, hiddenEnd] : [0, hiddenEnd, hiddenStart, length];
  const position = axis.positionAt(Math.min(x, axis.x0 + axis.width - 1e-6) - 1e-6);
  return {
    layer: s.section.address,
    section: s.id,
    system: axis.system,
    target: 'break',
    position: position >= 0 ? position : axis.dir === -1 ? hiddenEnd : hiddenStart - 1,
    extent: rangesLocation(axis.system, ranges),
    ...(count !== undefined ? { count } : {}),
  };
}
