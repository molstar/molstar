/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `matrix` renderer (seqspec:viewer#extension-points; seqspec:sources#levels-and-the-zoom-gate): positions × states
// as image tiles (`drawMatrixTiles`: at most 4,096 samples a tile, one image per tile of a shown run), a band per
// drawn state (its keyed order and filter, seqspec:format#keyed-sets), reduced to a row per device pixel where bands
// are thinner; hits in the draw's device frame (`matrixFrame`): the band a stretched tile shows in the device row, a
// reduced row's band whose color the tile drew (`matrixHitBand`), the hit's state the band's; the overlay outlines the
// hovered cell, and a data hover's row band is the hovered band.
import {
  type AxisView,
  drawMatrixTiles,
  type Hit,
  type HoverState,
  kernelEntry,
  MatrixBands,
  type MatrixState,
  MatrixTiles,
  matrixBandRows,
  matrixFrame,
  matrixHitBand,
  outlineMatrixCell,
} from '../draw/api.js';

const frame = new Int32Array(2);
/** The state under the point, from the device rows the draw fills; the position stays the hit. */
function hit(v: AxisView, s: MatrixState, _p: unknown, out: Hit): boolean {
  const b = out.block;
  if (b < 0 || s.states === 0 || !v.axis) return false;
  matrixFrame(v, frame);
  const y = (Math.floor(out.point[1] * v.dpr) | 0) - frame[0];
  const band = matrixHitBand(v, s, b, out.index - v.axis.base[b], y, frame[1]);
  if (band < 0) return false;
  out.state = s.stateOf[band];
  return true;
}

/** A data hover's row band: the hovered state's band alone; off one, the whole matrix. */
function rowHover(v: AxisView, s: MatrixState, _p: unknown, h: HoverState, out: Int32Array): boolean {
  const band = h.state >= 0 && h.state < s.bandOf.length ? s.bandOf[h.state] : -1;
  if (band >= 0) matrixBandRows(v, s, band, out);
  else {
    matrixFrame(v, out);
    out[1] += out[0];
  }
  return true;
}

/**
 * A band per state at a row's height times the row zoom, at most 18 px a band; the whole at least 24 px and at most
 * 160 px times the row zoom (480 at the default `maxRowScale`, 3); whole px, so the zoom does not box it.
 */
export const MatrixRenderer = kernelEntry<MatrixState>(
  'matrix',
  (s, _p, props, rowScale = 1) =>
    Math.round(Math.max(24, Math.min(160 * rowScale, s.states * Math.min(18, props.rowHeight * rowScale)))),
  {
    zooms: true,
    minHeight: 4,
    bands: MatrixBands,
    hit,
    draw: drawMatrixTiles,
    overlay: outlineMatrixCell,
    rowHover,
    cache: { create: () => new MatrixTiles(), dispose: (c) => (c as MatrixTiles).dispose() },
  },
);
