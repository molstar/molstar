/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `variants` renderer (seqspec:viewer#extension-points): markers colored by `alt`, stacked per position, a row per
// value when split; density as points. Linear only (wrapped takes its fallback): `drawPoints`, stacked. Streamed:
// outlined in its regions (ST3b; hit.ts hits them by position).
import { drawPoints, type Features, hitPoints, kernelEntry, outlinePoints, outlineStreamed } from "../draw/api.js";

export const VariantsRenderer = kernelEntry<Features>("variants", 32, {
  splitHeight: 8,
  // By drawn geometry, at its place in the stack (a boundary between positions, one past a block's end); else the column.
  hit: (v, f, _p, out) => hitPoints(v, f, out, true),
  hitsPastEnd: true,
  draw: (v, f, p) => drawPoints(v, f, p, true, false),
  tiny: (v, f, p) => drawPoints(v, f, p, true, true),
  overlay: (v, f) =>
    outlineStreamed(v, false) || outlinePoints(v, f, true, v.height < (VariantsRenderer.minHeight ?? 0)),
  sized: ["thickness"],
});
