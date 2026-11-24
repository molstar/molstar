/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ResolvedTrack } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import type { Box } from './layer-type.js';
import { StripPlacement, type Strips } from './layout/strips.js';
import type { Placement } from './layout/types.js';

// A section's strips as the tracks layout reads them (_spec/README.md "Render", Section layouts): a strip's span and
// band column in its frame, and a pending section's placement (no strips).

/** The union of a strip's bands in its frame (`ctx.span`: what a layer with `across` outlines). */
export function stripSpan(strips: Strips, i: number): Box {
  const a = strips.axes[i];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (let j = 0; j < strips.bands; j++) {
    const k = i * strips.bands + j;
    if (strips.bandH[k] <= 0) continue;
    if (strips.bandY[k] < y0) y0 = strips.bandY[k];
    if (strips.bandY[k] + strips.bandH[k] > y1) y1 = strips.bandY[k] + strips.bandH[k];
  }
  return y1 > y0 ? { x: a.x0, y: y0, width: a.width, height: y1 - y0 } : { x: a.x0, y: 0, width: a.width, height: 0 };
}

/**
 * The bands of the tracks a layer's `across` names on a strip, the first to the last and those between (`ctx.span`: a
 * frame's outline), in its frame; the whole strip's when it names none shown.
 */
export function acrossSpan(
  strips: Strips,
  i: number,
  shown: readonly ResolvedTrack[],
  across: readonly unknown[],
): Box {
  let a = -1;
  let b = -1;
  for (let j = 0; j < shown.length; j++) {
    const id = shown[j].address.slice(shown[j].address.lastIndexOf('/') + 1);
    if (id !== across[0] && id !== across[across.length - 1]) continue;
    if (a < 0) a = j;
    b = j;
  }
  if (a < 0) return stripSpan(strips, i);
  const x = strips.axes[i];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (let j = a; j <= b; j++) {
    const k = i * strips.bands + j;
    if (strips.bandH[k] <= 0) continue;
    y0 = Math.min(y0, strips.bandY[k]);
    y1 = Math.max(y1, strips.bandY[k] + strips.bandH[k]);
  }
  return y1 > y0 ? { x: x.x0, y: y0, width: x.width, height: y1 - y0 } : stripSpan(strips, i);
}

/** A strip's band column in its frame: its bands above the page's bottom and the gaps between them; none when empty. */
export function stripColumn(strips: Strips, i: number, pageHeight: number): Box | undefined {
  const a = strips.axes[i];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (let j = 0; j < strips.bands; j++) {
    const k = i * strips.bands + j;
    if (strips.bandH[k] <= 0 || strips.bandY[k] + strips.dy[i] >= pageHeight) continue;
    if (strips.bandY[k] < y0) y0 = strips.bandY[k];
    if (strips.bandY[k] + strips.bandH[k] > y1) y1 = strips.bandY[k] + strips.bandH[k];
  }
  return y1 > y0 ? { x: a.x0, y: y0, width: a.width, height: y1 - y0 } : undefined;
}

/** A pending section (no axis yet): no strips, its tracks' height kept. */
export function pendingPlacement(
  top: number,
  box: Box,
  heights: readonly number[],
  pads: readonly number[],
  config: RenderConfig,
  leads?: readonly number[],
): Placement {
  const p = new StripPlacement();
  let h = 0;
  for (let j = 0; j < heights.length; j++)
    h += (j ? config.trackGap : 0) + (leads?.[j] ?? 0) + pads[2 * j] + heights[j] + pads[2 * j + 1];
  p.strips.reset(heights.length);
  p.contentHeight = h;
  p.box = { x: config.inset, y: top, width: box.width, height: h };
  return p;
}
