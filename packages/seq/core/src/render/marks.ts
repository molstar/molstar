/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AxisTransform } from './axis.js';
import { NoColor, type PassWriter } from './draw-lists.js';
import type { Box } from './layer-type.js';

// Marks every layer type draws alike: spans as bands (placeholders, bars), and the hover's outline of a cell.

/** A rect per `x0, x1, cut` triple, over the band `y` to `y + h`, at least a device pixel wide. */
export function drawSpans(out: PassWriter, spans: readonly number[], y: number, h: number, fill: number, minW = 1) {
  for (let i = 0; i + 2 < spans.length; i += 3) out.rect(spans[i], y, Math.max(spans[i + 1] - spans[i], minW), h, fill);
}

/** The hover's outline around the cell of a position, inside the box, when it is in view. */
export function outlineCell(out: PassWriter, axis: AxisTransform, p: number, box: Box, color: number, width: number) {
  const x0 = axis.x(p);
  const x1 = axis.x(p + 1);
  if (!Number.isFinite(x0) || !Number.isFinite(x1)) return;
  const a = Math.max(axis.snap(x0), box.x);
  const b = Math.min(axis.snap(x1), box.x + box.width);
  if (b <= a) return;
  const half = width / 2;
  out.rect(a + half, box.y + half, Math.max(b - a - width, 0.5), box.height - width, NoColor, color, width);
}
