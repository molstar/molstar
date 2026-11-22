/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ribbons drawn and hit (_spec/README.md "Render", Bridges): each a closed path of two cubics with vertical
// tangents, clipped to the band; the density strips and the pending hatch along the edges; the lit outlined; the hit
// closed form. Part of the bridge chunk.

import { type Location, location } from '@molstar/seqx-core/query';
import { grown, NoColor, Op, PathCommand, type PassName, type PassWriter } from './draw-lists.js';
import type { Box, Emphasis, Hit } from './layer-type.js';
import { aOf, bOf, Line, R, type RibbonsLay } from './pairs-ribbons-lay.js';

// ---- draw ----

/** A ribbon's closed path (two cubics with vertical tangents joined by its ends) into a pass, through `args`. */
export function ribbonPath(out: PassWriter, xy: Float64Array, o: number, close: boolean) {
  const a = out.args;
  out.pass.op(Op.path);
  const p = out.pass.paths;
  p.begin();
  p.commands = grown(p.commands, p.commandCount + 5);
  p.points = grown(p.points, 2 * (p.pointCount + 8));
  const c = p.commands;
  const ya = xy[o + 2];
  const yb = xy[o + 5];
  const y1 = ya + (yb - ya) / 3;
  const y2 = ya + (2 * (yb - ya)) / 3;
  const q = p.points;
  let k = 2 * p.pointCount;
  c[p.commandCount++] = PathCommand.move;
  c[p.commandCount++] = PathCommand.cubic;
  q[k++] = xy[o];
  q[k++] = ya;
  q[k++] = xy[o];
  q[k++] = y1;
  q[k++] = xy[o + 3];
  q[k++] = y2;
  q[k++] = xy[o + 3];
  q[k++] = yb;
  p.pointCount += 4;
  if (close) {
    c[p.commandCount++] = PathCommand.line;
    c[p.commandCount++] = PathCommand.cubic;
    c[p.commandCount++] = PathCommand.close;
    q[k++] = xy[o + 4];
    q[k++] = yb;
    q[k++] = xy[o + 4];
    q[k++] = y2;
    q[k++] = xy[o + 1];
    q[k++] = y1;
    q[k++] = xy[o + 1];
    q[k++] = ya;
    p.pointCount += 4;
  }
  const i = p.n++;
  p.cmdCount[i] = p.commandCount - p.cmdStart[i];
  p.fill[i] = a[4];
  p.stroke[i] = a[5];
  p.strokeWidth[i] = a[6];
  p.dash[i] = a[7];
}

/** A line ribbon's middle curve: its corners collapsed onto the ends' middles (into `Mid`). */
const Mid = new Float64Array(R);
function middle(rib: Float64Array, o: number): Float64Array {
  Mid[0] = Mid[1] = (rib[o] + rib[o + 1]) / 2;
  Mid[2] = rib[o + 2];
  Mid[3] = Mid[4] = (rib[o + 3] + rib[o + 4]) / 2;
  Mid[5] = rib[o + 5];
  return Mid;
}

function clipBand(out: PassWriter, b: Box) {
  const g = out.args;
  g[0] = b.x;
  g[1] = b.y;
  g[2] = b.width;
  g[3] = b.height;
  out.clipArgs();
}

/** A filled rect through `args`. */
function rect(out: PassWriter, x: number, y: number, w: number, h: number, fill: number) {
  const g = out.args;
  g[0] = x;
  g[1] = y;
  g[2] = w;
  g[3] = h;
  g[4] = fill;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** The pending hatch's stripes: 4 px wide every 8 px, 4 px into the band from the edge (the pending ink). */
const HatchDepth = 4;

export function drawData(lay: RibbonsLay, out: PassWriter) {
  const g = out.args;
  clipBand(out, lay.band);
  for (let i = 0; i < lay.n; i++) {
    const line = lay.flags[i] & Line;
    g[4] = line ? NoColor : lay.color[i];
    g[5] = line ? lay.color[i] | 0xff : NoColor;
    g[6] = line ? 1 : 0;
    g[7] = 0;
    if (line) ribbonPath(out, middle(lay.rib, R * i), 0, false);
    else ribbonPath(out, lay.rib, R * i, true);
  }
  const b = lay.band;
  if (lay.densityN && lay.densityMax) {
    const dens = lay.tokens!.colors['seq-density-high'];
    for (let side = 0; side < 2; side++) {
      const y = side === 0 ? lay.edgeA : lay.edgeB;
      const down = y <= b.y + b.height / 2;
      for (let c = 0; c < lay.densityN; c++) {
        const v = lay.density[side * lay.densityN + c];
        if (!v) continue;
        const fill = ((dens & 0xffffff00) | Math.round(64 + (191 * v) / lay.densityMax)) >>> 0;
        rect(out, b.x + c / lay.dpr, down ? y : y - 3, 1 / lay.dpr, 3, fill);
      }
    }
  }
  const h = lay.hatch;
  if (h.length) {
    const pending = lay.tokens!.colors['seq-pending'];
    for (let i = 0; i + 3 < h.length; i += 4) {
      const y = h[i + 3] ? h[i + 2] : h[i + 2] - HatchDepth;
      for (let x = h[i]; x < h[i + 1]; x += 8) rect(out, x, y, Math.min(4, h[i + 1] - x), HatchDepth, pending);
    }
  }
  out.unclip();
}

/** Whether ribbon `i`'s ends (its hull) meet a Location's ranges (a hovered position, the selection). */
function endIn(
  lay: RibbonsLay,
  i: number,
  parts: readonly { readonly system: string; readonly ranges: ArrayLike<number> }[],
): boolean {
  const prep = lay.reads[lay.src[i]];
  const k = lay.pair[i];
  const h = lay.hull;
  const o = 4 * i;
  const as = aOf(prep, k);
  const bs = bOf(prep, k);
  for (let p = 0; p < parts.length; p++) {
    const a = parts[p].system === as;
    const b = parts[p].system === bs;
    if (!a && !b) continue;
    const r = parts[p].ranges;
    for (let j = 0; j + 1 < r.length; j += 2)
      if ((a && h[o] < r[j + 1] && r[j] < h[o + 1]) || (b && h[o + 2] < r[j + 1] && r[j] < h[o + 3])) return true;
  }
  return false;
}

/** The hovered pair in each read (−1 none). */
const Hovered = new Int32Array(2);

/** Whether ribbon `i` is lit: hovered (its element, or an end at the hovered position), or selected. */
function lit(lay: RibbonsLay, i: number, pass: PassName, e: Emphasis): boolean {
  const prep = lay.reads[lay.src[i]];
  const k = lay.pair[i];
  if (pass === 'hover') {
    if (k === Hovered[lay.src[i]]) return true;
    const h = e.hover;
    if (!h || h.element) return false;
    if (h.position !== undefined) {
      const p = Math.floor(h.position);
      const o = 4 * i;
      if (h.system === aOf(prep, k) && p >= lay.hull[o] && p < lay.hull[o + 1]) return true;
      if (h.system === bOf(prep, k) && p >= lay.hull[o + 2] && p < lay.hull[o + 3]) return true;
    }
    const projected = e.projected?.hover;
    return projected !== undefined && endIn(lay, i, projected.parts);
  }
  const keys = e.elements?.keys;
  if (keys) {
    const key = prep.data!.keys[k];
    for (let j = 0; j < keys.length; j++) if (keys[j].owner === prep.owner && keys[j].key === key) return true;
  }
  return e.selection !== undefined && endIn(lay, i, e.selection.parts);
}

export function drawEmphasis(lay: RibbonsLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const h = e.hover;
  Hovered[0] = Hovered[1] = -1;
  if (pass === 'hover') {
    if (!h && !e.projected?.hover) return;
    if (h?.element) {
      let any = false;
      for (let r = 0; r < lay.reads.length && r < 2; r++) {
        const prep = lay.reads[r];
        Hovered[r] = h.element.owner === prep.owner ? (prep.index.get(h.element.key) ?? -1) : -1;
        any ||= Hovered[r] >= 0;
      }
      if (!any) return;
    }
  } else if (!e.elements?.keys.length && !e.selection) return;
  const tokens = lay.tokens!;
  const config = lay.config!;
  const ink = pass === 'hover' ? tokens.colors['seq-hover'] : tokens.colors['msui-accent'];
  const width = pass === 'hover' ? config.outline.hover + 1 : config.outline.selection + 1;
  const g = out.args;
  clipBand(out, lay.band);
  for (let i = 0; i < lay.n; i++) {
    if (!lit(lay, i, pass, e)) continue;
    const line = lay.flags[i] & Line;
    const hovered = lay.pair[i] === Hovered[lay.src[i]];
    // the hovered ribbon at full opacity with the hover's ink as its outline; the lit and selected outlined
    g[4] = pass === 'hover' && hovered && !line ? lay.color[i] | 0xff : NoColor;
    g[5] = ink;
    g[6] = width;
    g[7] = 0;
    if (line) ribbonPath(out, middle(lay.rib, R * i), 0, false);
    else ribbonPath(out, lay.rib, R * i, true);
  }
  out.unclip();
}

// ---- hit ----

/**
 * The ribbon under a page point, closed form (_spec/README.md "Render", Bridges): at the point's height `u` each
 * ribbon's sides are `x_a + (x_b − x_a)·S(u)`; the point is in it when its x lies between them (widened to
 * `pairs.hitPx` for a thin one), scanned last drawn first (the narrowest on top). Nothing allocated.
 */
export function hitRibbons(lay: RibbonsLay, at: Float64Array): Hit | undefined {
  const x = at[0];
  const y = at[1];
  const b = lay.band;
  if (!lay.prep || y < b.y || y >= b.y + b.height || x < b.x || x >= b.x + b.width) return undefined;
  const r = lay.config!.pairs.hitPx;
  const rib = lay.rib;
  for (let i = lay.n - 1; i >= 0; i--) {
    const o = R * i;
    const ya = rib[o + 2];
    const yb = rib[o + 5];
    const u = (y - ya) / (yb - ya);
    if (!(u >= 0 && u <= 1)) continue;
    const s = u * u * (3 - 2 * u);
    let lo = rib[o] + (rib[o + 3] - rib[o]) * s;
    let hi = rib[o + 1] + (rib[o + 4] - rib[o + 1]) * s;
    if (lay.flags[i] & Line) lo = hi = (lo + hi) / 2;
    if (lo > hi) {
      const t = lo;
      lo = hi;
      hi = t;
    }
    if (hi - lo < r) {
      const m = (lo + hi) / 2;
      lo = m - r / 2;
      hi = m + r / 2;
    }
    if (x >= lo && x <= hi) return ribbonHit(lay, i);
  }
  return undefined;
}

/** The two-part (or one-part) location of ends on systems `as`, `bs`. */
function endsLocation(as: string, bs: string, a0: number, a1: number, b0: number, b1: number): Location {
  return location({
    parts:
      as === bs
        ? [{ system: as, ranges: a0 <= b0 ? [a0, a1, b0, b1] : [b0, b1, a0, a1] }]
        : [
            { system: as, ranges: [a0, a1] },
            { system: bs, ranges: [b0, b1] },
          ],
  });
}

/**
 * A ribbon's hit, the drawn run: its pair's element (a merged ribbon's first block) and its ends' extent (a merged
 * ribbon's hull), its bridge; made once a pair or merged ribbon (`elements`, `extents`).
 */
function ribbonHit(lay: RibbonsLay, i: number): Hit {
  const prep = lay.reads[lay.src[i]];
  const d = prep.data!;
  const k = lay.pair[i];
  let element = prep.elements[k];
  if (!element) element = prep.elements[k] = { owner: prep.owner, key: d.keys[k] };
  const as = aOf(prep, k);
  const bs = bOf(prep, k);
  let extent: Location | undefined;
  if (lay.merged[i]) {
    const h = lay.hull;
    const o = 4 * i;
    extent = lay.extents[i] ??= endsLocation(as, bs, h[o], h[o + 1], h[o + 2], h[o + 3]);
  } else extent = prep.extents[k] ??= endsLocation(as, bs, d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k]);
  const h = lay.scratch;
  h.layer = lay.key;
  h.section = lay.bridgeId;
  h.system = as;
  h.element = element;
  h.extent = extent;
  h.position = undefined;
  h.ySystem = undefined;
  h.yPosition = undefined;
  h.target = undefined;
  h.bin = prep.level;
  h.bridge = lay.bridge;
  return h;
}
