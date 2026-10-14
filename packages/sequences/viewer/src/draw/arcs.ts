/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A pair's arc as the `pairs` draw and overlay stroke it (seqspec:viewer#extension-points `pairs`), part of the
// renderer API: the curve from `Arcs` (one geometry for the draw, the overlay and the hit), device px written through
// `verb` so no double crosses a call; an end out of the row or filtered out a tick toward its partner, an end past the
// data area a wedge.

import type { PairsState } from "../prepare.js";
import { ARC_PAST, Arcs } from "../renderer-helpers.js";
import { NO_COLOR, type PathBatch, PathMode, PathVerb, SkinToken } from "./backend.js";
import { keyAlpha, keyColor } from "./density.js";
import type { DrawView } from "./view.js";

/** A pair's color, or the quiet ink where its theme gives none. */
export const inkOf = (key: number) => (key === NO_COLOR ? SkinToken.inkQuiet : keyColor(key));

// Enum-like members read once (a per-frame path: see the heap guards).
const MOVE = PathVerb.move;
const LINE = PathVerb.line;
const QUAD = PathVerb.quad;
const STROKE = PathMode.stroke;
const FILL = PathMode.fill;
/**
 * The window's t range of the parabola being pushed. Held in a typed array: a literal 0 or 1 is a small integer, and
 * comparisons that saw only those (every arc in view) were optimized for them, then deoptimized at the first fraction.
 */
const win = new Float64Array(2);

/**
 * Pair k's arc (after `Arcs.frame`) into `arcs`, stroked `width` half device px, in its color key, its apex by
 * `height` (`Arcs.apex`): an end on another row of a stacked axis or filtered out is a tick toward its partner. The
 * curve is `Arcs.heightAt`'s, the parabola over its virtual ends (each taken at most a margin past the window's edge),
 * within the band at every zoom and moving continuously as a pan moves its ends; only its piece within the data area's
 * x is pushed. An end beyond the data area gets a wedge into `wedges` where the curve meets that edge. Each pair is
 * its own subpath. False: nothing of it is drawn here. The draw's and the overlay's (the same curve, stronger).
 */
export function pushArc(
  v: DrawView,
  s: PairsState,
  k: number,
  width: number,
  key: number,
  arcs: PathBatch,
  wedges: PathBatch,
  height = -1,
): boolean {
  const f = s.frame;
  const dpr = v.dpr;
  const both = Arcs.ends(v, s, k);
  const xa = s.xs[0];
  const xb = s.xs[1];
  if (!both && Number.isNaN(xa) && Number.isNaN(xb)) return false;
  if (both && (Math.max(xa, xb) < v.left - 1 || Math.min(xa, xb) > v.right + 1)) return false;
  let o = 0;
  if (!both) {
    const x = !Number.isNaN(xa) ? xa : xb;
    arcs.begin(inkOf(key), keyAlpha(key), STROKE, width);
    o = arcs.verb(MOVE);
    arcs.coords[o] = x * dpr;
    arcs.coords[o + 1] = f[0] * dpr;
    o = arcs.verb(LINE);
    arcs.coords[o] = x * dpr;
    arcs.coords[o + 1] = (f[0] + f[1] * Math.min(f[2], 10)) * dpr;
    return true;
  }
  const ring = v.ring;
  if (ring) {
    // On rings, a chord (seqspec:viewer#the-render-backend): a quadratic Bézier from end to end through a point on
    // their bisector, toward the center by their separation (opposite ends: through it), in the unrolled frame the ring
    // map takes it from.
    const dx = Math.abs(xb - xa);
    const sep = Math.min(dx, ring.w - dx) / (ring.w / 2);
    const rEnd = ring.r - (f[0] - ring.top) * ring.scale;
    arcs.begin(inkOf(key), keyAlpha(key), STROKE, width);
    o = arcs.verb(MOVE);
    arcs.coords[o] = xa * dpr;
    arcs.coords[o + 1] = f[0] * dpr;
    o = arcs.verb(QUAD);
    arcs.coords[o] = ((xa + xb) / 2 + (dx > ring.w / 2 ? ring.w / 2 : 0)) * dpr;
    arcs.coords[o + 1] = (ring.top + (ring.r - rEnd * (1 - sep)) / ring.scale) * dpr;
    arcs.coords[o + 2] = xb * dpr;
    arcs.coords[o + 3] = f[0] * dpr;
    return true;
  }
  Arcs.apex(s, height);
  // The parabola over the virtual ends (`Arcs.virtual`), only its piece within the data area's x (all of it for ends
  // in one cell): its ends on the curve, its control the curve's blossom at (s0, s1).
  const a = f[9];
  const y0 = f[0];
  const rise = f[1] * a;
  const va = s.xs[2];
  const span = s.xs[3] - va;
  win[0] = 0;
  win[1] = 1;
  if (span !== 0) {
    win[0] = Math.max(0, (v.left - ARC_PAST - va) / span);
    win[1] = Math.min(1, (v.right + ARC_PAST - va) / span);
  }
  const s0 = win[0];
  const s1 = win[1];
  if (!(s1 > s0) && span !== 0) return false;
  arcs.begin(inkOf(key), keyAlpha(key), STROKE, width);
  o = arcs.verb(MOVE);
  arcs.coords[o] = (va + span * s0) * dpr;
  arcs.coords[o + 1] = (y0 + rise * 4 * s0 * (1 - s0)) * dpr;
  o = arcs.verb(QUAD);
  arcs.coords[o] = (va + (span * (s0 + s1)) / 2) * dpr;
  arcs.coords[o + 1] = (y0 + rise * 2 * ((1 - s0) * s1 + s0 * (1 - s1))) * dpr;
  arcs.coords[o + 2] = (va + span * s1) * dpr;
  arcs.coords[o + 3] = (y0 + rise * 4 * s1 * (1 - s1)) * dpr;
  for (let side = -1; side <= 1; side += 2) {
    // A wedge at the edge an end lies beyond, pointing out of the view, where the curve meets it.
    const beyond = side < 0 ? Math.min(xa, xb) < v.left : Math.max(xa, xb) > v.right;
    if (!beyond) continue;
    const edge = side < 0 ? v.left : v.right;
    f[10] = edge;
    Arcs.heightAt(s);
    const h = f[11];
    if (h < 0) continue;
    const ex = edge * dpr;
    const y = y0 + f[1] * h;
    wedges.begin(inkOf(key), keyAlpha(key), FILL);
    o = wedges.verb(MOVE);
    wedges.coords[o] = ex;
    wedges.coords[o + 1] = y * dpr;
    o = wedges.verb(LINE);
    wedges.coords[o] = ex - side * 5 * dpr;
    wedges.coords[o + 1] = (y - 3) * dpr;
    o = wedges.verb(LINE);
    wedges.coords[o] = ex - side * 5 * dpr;
    wedges.coords[o + 1] = (y + 3) * dpr;
    wedges.close();
  }
  return true;
}
