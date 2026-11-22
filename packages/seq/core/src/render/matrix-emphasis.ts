/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LocationPart } from '@molstar/seqx-core/query';
import { NoColor, type PassName, type PassWriter } from './draw-lists.js';
import type { Emphasis } from './layer-type.js';
import type { MatrixLay } from './matrix.js';
import { partOn } from './marks.js';

// The matrix layer's emphasis (_spec/README.md "Render", The matrix layer; ux-foundations §6, the hover contract): a
// member of its states system selected or hovered anywhere (the kinship grid's label, an alignment's row) lights its
// state's row across the band: the selection a tint with an outline, the hover an outline with a light tint, what is
// projected the same, lighter. Rows in view only, cut at the band: O(states in view).

/** A band over states `[s, s + 1)` for each range of `ranges` in view, as a tint and an outline. */
function bands(
  lay: MatrixLay,
  out: PassWriter,
  ranges: ArrayLike<number>,
  fill: number,
  stroke: number,
  width: number,
) {
  const { box } = lay;
  const g = out.args;
  g[0] = box.x;
  g[1] = box.y;
  g[2] = box.width;
  g[3] = box.height;
  out.clipArgs();
  for (let k = 0; k + 1 < ranges.length; k += 2) {
    const s0 = Math.max(ranges[k], lay.first);
    const s1 = Math.min(ranges[k + 1], lay.last);
    for (let s = s0; s < s1; s++) {
      const y = lay.box.y + s * lay.rowScale - lay.bandTop;
      g[0] = box.x;
      g[1] = y;
      g[2] = box.width;
      g[3] = lay.rowScale;
      g[4] = fill;
      g[5] = stroke;
      g[6] = width;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
  out.unclip();
}

/** A color at `a` (0 to 1) of its own opacity. */
function faded(color: number, a: number): number {
  return ((color & 0xffffff00) | Math.round(Math.min(1, a) * (color & 0xff))) >>> 0;
}

const One = new Int32Array(2);

function hoverRows(lay: MatrixLay, out: PassWriter, e: Emphasis) {
  const { config, tokens } = lay;
  if (!config || !tokens) return;
  const ink = tokens.colors['seq-hover'];
  const h = e.hover;
  // a hit of another layer naming a member of the states system (its own cell is the view's position box)
  if (h && h.layer !== lay.key) {
    if (h.system === lay.statesSystem && h.position !== undefined) {
      One[0] = h.position;
      One[1] = h.position + 1;
      bands(lay, out, One, faded(ink, 0.1), ink, config.outline.hover);
    }
    if (h.ySystem === lay.statesSystem && h.yPosition !== undefined && h.yPosition !== h.position) {
      One[0] = h.yPosition;
      One[1] = h.yPosition + 1;
      bands(lay, out, One, faded(ink, 0.1), ink, config.outline.hover);
    }
    return;
  }
  const p = e.projected;
  const part = h ? undefined : partOn(p?.hover, lay.statesSystem);
  if (part && p?.hoverSection !== lay.section)
    bands(lay, out, part.ranges, NoColor, faded(ink, config.emphasis.projected), config.outline.hover);
}

/** Selection rows: its own part, else the projected one, lighter. */
function selectionRows(lay: MatrixLay, out: PassWriter, e: Emphasis) {
  const { config, tokens } = lay;
  if (!config || !tokens) return;
  const own = partOn(e.selection, lay.statesSystem);
  const part: LocationPart | undefined = own ?? partOn(e.projected?.selection, lay.statesSystem);
  if (!part) return;
  const k = own ? 1 : config.emphasis.projected;
  const accent = tokens.colors['msui-accent'];
  bands(
    lay,
    out,
    part.ranges,
    faded(accent, config.emphasis.selectionTint * k),
    faded(accent, config.emphasis.selectionOutline * k),
    config.outline.selection,
  );
}

export function drawMatrixEmphasis(lay: MatrixLay, pass: PassName, out: PassWriter, e: Emphasis) {
  if (!lay.tokens || !lay.axis || !lay.states) return;
  if (pass === 'selection') selectionRows(lay, out, e);
  else if (pass === 'hover') hoverRows(lay, out, e);
}
