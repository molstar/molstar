/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AxisTransform } from './axis.js';
import { grown, type PassWriter } from './draw-lists.js';
import type { LayoutContext } from './layer-type.js';
import type { RangesLay } from './ranges-lay.js';

// A ranges layer's connectors (_spec/README.md "Render", Layer types; axis.md, An item split by a path is joined by a
// connector): an item projected through a path into pieces with shown positions between them that the path does not
// reach (a domain through the codon map over a gene whose intron is shown) is joined by a hairline at its lane's
// mid-height, from one piece's end to the next's start, in the item's colour. Across a break (zero units wide) the
// two x are the same and nothing is drawn. The line is part of its item: a hit on it is the item's, and the hover
// outlines it with the item. The layout writes the connectors and the draw reads them: no number is made a frame.

/** The least reach, CSS px, of a connector: two ends this close are a break's (zero units wide) and join nothing. */
const MinReach = 1e-6;

/** The layout's settings of the layer's connectors: its param, the preset's width, ink and alpha. */
export function startConnectors(lay: RangesLay, ctx: LayoutContext, on: boolean): void {
  const c = ctx.config.ranges.connector;
  lay.cn = 0;
  lay.cwidth = c.width;
  lay.cink = c.ink === 'muted' ? ctx.tokens.colors['msui-text-muted'] : 0;
  lay.calpha = on ? c.alpha : 0;
}

/**
 * The connectors of one item whose pieces (`[start, end, …]` ascending, on the section's system) are laid in `lane`:
 * one between each two neighbours with positions shown between them, clipped to the data box.
 */
export function connect(
  lay: RangesLay,
  axis: AxisTransform,
  pieces: readonly number[],
  lane: number,
  item: number,
  left: number,
  right: number,
): void {
  if (lay.calpha <= 0) return;
  for (let k = 0; k + 3 < pieces.length; k += 2) {
    const e = pieces[k + 1];
    const s = pieces[k + 2];
    if (!(s > e)) continue;
    const xa = axis.x(e);
    const xb = axis.x(s);
    if (!(Math.abs(xb - xa) > MinReach)) continue;
    const x0 = Math.max(Math.min(xa, xb), left);
    const x1 = Math.min(Math.max(xa, xb), right);
    if (!(x1 > x0)) continue;
    const j = lay.cn++;
    if (j >= lay.cx0.length) {
      lay.cx0 = grown(lay.cx0, j + 1);
      lay.cx1 = grown(lay.cx1, j + 1);
      lay.clane = grown(lay.clane, j + 1);
      lay.citem = grown(lay.citem, j + 1);
      lay.cdir = grown(lay.cdir, j + 1);
      lay.canchor = grown(lay.canchor, j + 1);
    }
    lay.cx0[j] = x0;
    lay.cx1[j] = x1;
    lay.clane[j] = lane;
    lay.citem[j] = item;
    lay.cdir[j] = 0;
  }
}

/** The connector under `x` in `lane` (its lane's whole height hits it), or −1. */
export function connectorAt(lay: RangesLay, lane: number, x: number): number {
  for (let j = 0; j < lay.cn; j++) if (lay.clane[j] === lane && x >= lay.cx0[j] && x < lay.cx1[j]) return j;
  return -1;
}

/** Chevrons along a hairline: one every this many CSS px, this far from its ends. */
const ChevronStep = 12;
/** A chevron's three points, written per chevron (kept). */
const Chevron = new Float64Array(6);

/**
 * The data pass: each connector a line `cwidth` thick at its lane's middle, on the device grid, in its item's colour;
 * a gene line's with chevrons pointing along its strand (`cdir`), where the lane is tall enough to read them.
 */
export function drawConnectors(lay: RangesLay, out: PassWriter): void {
  const dpr = lay.axis.dpr;
  const w = Math.max(lay.cwidth, 1 / dpr);
  const half = Math.min(3, lay.laneHeight / 4);
  const a = out.args;
  for (let j = 0; j < lay.cn; j++) {
    const base = lay.cink !== 0 ? lay.cink : lay.fills[lay.citem[j]];
    const color = ((base & 0xffffff00) | Math.round((base & 255) * lay.calpha)) >>> 0;
    const y = Math.round((lay.laneY + lay.clane[j] * lay.laneHeight + (lay.laneHeight - w) / 2) * dpr) / dpr;
    out.fill(lay.cx0[j], y, lay.cx1[j] - lay.cx0[j], w, color);
    const dir = lay.cdir[j];
    if (dir === 0 || half < 1.5) continue;
    const cy = y + w / 2;
    // anchored to the hairline's content (its unclipped left end on the screen): a chevron every step, inside the line
    const tip = dir === 1 ? half : -half;
    const from = lay.canchor[j] + ChevronStep / 2;
    const skip = Math.max(0, Math.ceil((lay.cx0[j] + half - from) / ChevronStep));
    for (let x = from + skip * ChevronStep; x + half <= lay.cx1[j] - half; x += ChevronStep) {
      Chevron[0] = x - tip / 2;
      Chevron[1] = cy - half;
      Chevron[2] = x + tip / 2;
      Chevron[3] = cy;
      Chevron[4] = x - tip / 2;
      Chevron[5] = cy + half;
      a[4] = 0;
      a[5] = color;
      a[6] = 1;
      a[7] = 0;
      out.pathArgs(Chevron, 3, false);
    }
  }
}

/** The hover pass: the hovered item's connectors in the hover ink (read here), as thick as its outline. */
export function hoverConnectors(lay: RangesLay, out: PassWriter, key: string | number): void {
  inkConnectors(lay, out, key, lay.hoverWidth, 0);
}

/** The selection pass: the connectors of the items whose keys `keys` holds for this layer's owner, in the accent. */
export function selectConnectors(
  lay: RangesLay,
  out: PassWriter,
  keys: readonly { readonly owner: string; readonly key: string | number }[],
): void {
  for (let k = 0; k < keys.length; k++)
    if (keys[k].owner === lay.owner) inkConnectors(lay, out, keys[k].key, lay.cwidth, 1);
}

/** The connectors of one item in the hover ink (`which` 0) or the accent (1), `width` thick, on the device grid. */
function inkConnectors(lay: RangesLay, out: PassWriter, key: string | number, width: number, which: number): void {
  const dpr = lay.axis.dpr;
  const w = Math.max(width, 1 / dpr);
  const color = which === 0 ? lay.tokens.colors['seq-hover'] : lay.tokens.colors['msui-accent'];
  for (let j = 0; j < lay.cn; j++) {
    if (lay.keys[lay.citem[j]] !== key) continue;
    const y = Math.round((lay.laneY + lay.clane[j] * lay.laneHeight + (lay.laneHeight - w) / 2) * dpr) / dpr;
    out.fill(lay.cx0[j], y, lay.cx1[j] - lay.cx0[j], w, color);
  }
}
