/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LocationPart } from '@molstar/seqx-core/query';
import { NoColor, type PassName, type PassWriter } from './draw-lists.js';
import type { Emphasis } from './layer-type.js';
import { partOn, rectOver } from './marks.js';
import type { MatrixLay } from './matrix.js';

// The matrix layer's emphasis (_spec/README.md "Render", The matrix layer; ux-foundations §6, the hover contract): a
// member of its states system selected or hovered anywhere (the kinship grid's label, an alignment's row) lights its
// state's row across the band: the selection a tint with an outline, the hover an outline with a light tint, what is
// projected the same, lighter. Rows in view only, cut at the band: O(states in view).

/**
 * What the bands paint with: a fill color and its opacity, a stroke color and its opacity, a width (`Src`, set before
 * `bands`). The colors (above 2^30, no small integers) and the opacities go through typed arrays: a number passed to a
 * call that is not inlined is boxed, and an emphasis redraw allocates nothing once warm (as `grid-emphasis.ts`).
 */
const Src = new Float64Array(5);
const Ink = new Float64Array(3);

function ink() {
  const f = Src[0];
  const s = Src[2];
  Ink[0] = Src[1] > 0 ? ((f & 0xffffff00) | Math.round(Math.min(1, Src[1]) * (f & 0xff))) >>> 0 : NoColor;
  Ink[1] = Src[3] > 0 ? ((s & 0xffffff00) | Math.round(Math.min(1, Src[3]) * (s & 0xff))) >>> 0 : NoColor;
  Ink[2] = Src[4];
}

/** Paints states `[a, b)` in view as one band (contiguous states have no outline between them). */
function band(lay: MatrixLay, out: PassWriter, a: number, b: number) {
  const s0 = Math.max(a, lay.first);
  const s1 = Math.min(b, lay.last);
  if (s1 <= s0) return;
  const { box } = lay;
  const g = out.args;
  g[0] = box.x;
  g[1] = box.y + s0 * lay.rowScale - lay.bandTop;
  g[2] = box.width;
  g[3] = (s1 - s0) * lay.rowScale;
  g[4] = Ink[0];
  g[5] = Ink[1];
  g[6] = Ink[2];
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** A band over each run of `ranges` (contiguous ranges merged) in view, as the tint and outline in `Src`. */
function bands(lay: MatrixLay, out: PassWriter, ranges: ArrayLike<number>) {
  ink();
  const { box } = lay;
  const g = out.args;
  g[0] = box.x;
  g[1] = box.y;
  g[2] = box.width;
  g[3] = box.height;
  out.clipArgs();
  let from = -1;
  let to = -1;
  for (let k = 0; k + 1 < ranges.length; k += 2) {
    if (ranges[k] === to) {
      to = ranges[k + 1];
      continue;
    }
    if (from >= 0) band(lay, out, from, to);
    from = ranges[k];
    to = ranges[k + 1];
  }
  if (from >= 0) band(lay, out, from, to);
  out.unclip();
}

const One = new Int32Array(2);

function hoverRows(lay: MatrixLay, out: PassWriter, e: Emphasis) {
  const { config, tokens } = lay;
  if (!config || !tokens) return;
  const hover = tokens.colors['seq-hover'];
  const h = e.hover;
  // a hit of another layer naming a member of the states system (its own cell is the view's position box)
  if (h && h.layer !== lay.key) {
    if (h.system === lay.statesSystem && h.position !== undefined) {
      Src[0] = hover;
      Src[1] = config.emphasis.hoverTint;
      Src[2] = hover;
      Src[3] = 1;
      Src[4] = config.outline.hover;
      One[0] = h.position;
      One[1] = h.position + 1;
      bands(lay, out, One);
    }
    if (h.ySystem === lay.statesSystem && h.yPosition !== undefined && h.yPosition !== h.position) {
      Src[0] = hover;
      Src[1] = config.emphasis.hoverTint;
      Src[2] = hover;
      Src[3] = 1;
      Src[4] = config.outline.hover;
      One[0] = h.yPosition;
      One[1] = h.yPosition + 1;
      bands(lay, out, One);
    }
    return;
  }
  const p = e.projected;
  const part = h ? undefined : partOn(p?.hover, lay.statesSystem);
  if (part && p?.hoverSection !== lay.section) {
    Src[0] = 0;
    Src[1] = 0;
    Src[2] = hover;
    Src[3] = config.emphasis.projected;
    Src[4] = config.outline.hover;
    bands(lay, out, part.ranges);
  }
}

/**
 * Selection rows: its own part, else the projected one, lighter (never the sides of a rect selected on this section's
 * own cells: it lights its cells, as the grid's `rectOver`).
 */
function selectionRows(lay: MatrixLay, out: PassWriter, e: Emphasis) {
  const { config, tokens } = lay;
  if (!config || !tokens) return;
  const own = partOn(e.selection, lay.statesSystem);
  const part: LocationPart | undefined =
    own ??
    (rectOver(e.selection, lay.system, lay.statesSystem)
      ? undefined
      : partOn(e.projected?.selection, lay.statesSystem));
  if (!part) return;
  const k = own ? 1 : config.emphasis.projected;
  const accent = tokens.colors['msui-accent'];
  Src[0] = accent;
  Src[1] = config.emphasis.selectionTint * k;
  Src[2] = accent;
  Src[3] = config.emphasis.selectionOutline * k;
  Src[4] = config.outline.selection;
  bands(lay, out, part.ranges);
}

export function drawMatrixEmphasis(lay: MatrixLay, pass: PassName, out: PassWriter, e: Emphasis) {
  if (!lay.tokens || !lay.axis || !lay.states) return;
  if (pass === 'selection') selectionRows(lay, out, e);
  else if (pass === 'hover') hoverRows(lay, out, e);
}
