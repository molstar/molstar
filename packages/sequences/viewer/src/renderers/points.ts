/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `points` renderer (seqspec:viewer#extension-points): glyphs; underlines under letters when wrapped; density past
// 64 in view (below two device px a cell, and on the tiny path): `drawPoints`. Streamed: by region as `variants`,
// outlined in its regions (hit.ts hits them by position).
import {
  drawPoints,
  type Features,
  hitLineLanes,
  hitPoints,
  kernelEntry,
  outlinePoints,
  outlineStreamed,
} from '../draw/api.js';

export const PointsRenderer = kernelEntry<Features>('points', 20, {
  geometries: ['linear', 'wrapped'],
  hit: (v, f, _p, out) => hitLineLanes(v, f, out) || hitPoints(v, f, out),
  hitsPastEnd: true,
  draw: (v, f, p) => drawPoints(v, f, p, false, false),
  tiny: (v, f, p) => drawPoints(v, f, p, false, true),
  overlay: (v, f) =>
    outlineStreamed(v, false) || outlinePoints(v, f, false, v.height < (PointsRenderer.minHeight ?? 0)),
  sized: ['radius'],
});
