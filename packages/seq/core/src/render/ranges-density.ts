/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesDensity } from '@molstar/seqx-core/data';
import type { AxisTransform } from './axis.js';
import type { ImageData } from './draw-lists.js';
import { type LayerSource, toSection } from './layer-data.js';
import { mixColors, type SeqTokens } from './tokens.js';
import { itemPiece, mappedPieces, type RangesLay } from './ranges-lay.js';

// A ranges layer's heat strip (_spec/README.md "Render", Layer types): its items and a source's density bins counted
// a device column and painted low to high in the density tokens (`items.ts` reads the strip too).

/**
 * Adds the strip's `v` over the device columns data positions `[s, e)` land on in view: the positions are integers, the
 * count in the strip, `column` inlined (a number made and passed to a call that is not inlined is boxed: an allocation
 * a bin).
 */
export function addToStrip(strip: StripBuffers, src: LayerSource, axis: AxisTransform, s: number, e: number) {
  let pieces = itemPiece;
  if (src.path) pieces = toSection(src, s, e, mappedPieces);
  else {
    itemPiece[0] = s;
    itemPiece[1] = e;
  }
  const diff = strip.diff;
  const v = strip.v;
  const dpr = axis.dpr;
  const last = axis.columns - 1;
  for (let k = 0; k + 1 < pieces.length; k += 2) {
    axis.piecesIn(pieces[k], pieces[k + 1]);
    const ka = axis.ka;
    const kb = axis.kb;
    for (let a = ka; a <= kb; a++) {
      if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
      const f0 = Math.floor((axis.sx0 - axis.x0) * dpr);
      const f1 = Math.floor((axis.sx1 - 1e-6 - axis.x0) * dpr);
      const c0 = f0 < 0 ? 0 : f0 >= last ? last : f0;
      const c1 = Math.max(c0, f1 < 0 ? 0 : f1 >= last ? last : f1);
      diff[c0] += v;
      diff[c1 + 1] -= v;
    }
  }
}

/**
 * One heat strip: items (all, or with `overflowOnly` those past the lanes) and a source's density bins, counted a
 * device column, colored low to high. Its buffers are kept while the box's columns do not change.
 */
export function densityStrip(
  lay: RangesLay,
  src: LayerSource,
  axis: AxisTransform,
  densities: readonly { table: RangesDensity; level: number }[],
  tokens: SeqTokens,
  overflowOnly: boolean,
  share: number,
) {
  const strip = lay.strip;
  resetStrip(strip, axis.columns);
  strip.v = 1;
  for (let i = 0; i < lay.items; i++)
    if (!overflowOnly || lay.lanes[i] < 0)
      addToStrip(strip, src, axis, lay.starts[i], Math.max(lay.ends[i], lay.starts[i] + 1));
  for (let d = 0; d < densities.length; d++) {
    const { table, level } = densities[d];
    for (let i = 0; i < table.bin.length; i++) {
      strip.v = table.count[i];
      addToStrip(strip, src, axis, table.bin[i] * level, (table.bin[i] + 1) * level);
    }
  }
  lay.density = paintStrip(strip, axis.columns, tokens, share, axis);
}

/** A column's centre, held in an array (a fractional number passed to a call that is not inlined is boxed). */
const GapAt = new Float64Array(1);

/** A heat strip's buffers: counts by difference over the device columns, and the image they are painted into. */
export interface StripBuffers {
  image?: ImageData;
  diff: Float64Array;
  counts: Float64Array;
  v: number;
}

/** The strip emptied for `cols` device columns (its buffers kept while their number does not change). */
export function resetStrip(strip: StripBuffers, cols: number) {
  if (strip.diff.length !== cols + 1) {
    strip.diff = new Float64Array(cols + 1);
    strip.counts = new Float64Array(cols);
  } else strip.diff.fill(0);
  if (!strip.image || strip.image.width !== cols)
    strip.image = { width: cols, height: 1, data: new Uint8Array(4 * cols) };
  else strip.image.data.fill(0);
}

/**
 * The strip's counts (`addToStrip` made them) painted low to high in the density tokens: its image. `share` is the part
 * of the ramp the busiest column reaches (1: the whole ramp; a thin band uses the light end, `ThinStripShare`).
 */
export function paintStrip(
  strip: StripBuffers,
  cols: number,
  tokens: SeqTokens,
  share = 1,
  axis?: AxisTransform,
): ImageData {
  const { diff, counts } = strip;
  const data = strip.image!.data;
  let max = 0;
  let run = 0;
  // a group's gap (a text's) holds no position: its columns stay clear, as the cells beside it leave it (every drawn
  // column then has a position to hit, _spec/README.md "Render", Emphasis, The hover contract)
  const gaps = axis !== undefined && axis.gsize > 0;
  const dpr = axis?.dpr ?? 1;
  const x0 = axis?.x0 ?? 0;
  for (let c = 0; c < cols; c++) {
    run += diff[c];
    counts[c] = run;
    if (gaps) {
      GapAt[0] = x0 + (c + 0.5) / dpr;
      if (axis!.inGapOf(GapAt)) counts[c] = 0;
    }
    if (counts[c] > max) max = counts[c];
  }
  const low = tokens.colors['seq-density-low'];
  const high = tokens.colors['seq-density-high'];
  for (let c = 0; c < cols; c++) {
    if (counts[c] <= 0) continue;
    const color = mixColors(low, high, max > 0 ? (share * counts[c]) / max : 0);
    data[4 * c] = color >>> 24;
    data[4 * c + 1] = (color >>> 16) & 0xff;
    data[4 * c + 2] = (color >>> 8) & 0xff;
    data[4 * c + 3] = color & 0xff;
  }
  return strip.image!;
}
