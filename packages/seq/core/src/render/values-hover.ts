/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { perPosition } from './config.js';
import { NoColor, type PassWriter, Round } from './draw-lists.js';
import type { Emphasis } from './layer-type.js';
import { hoverCell } from './marks.js';
import { Fill, hoverDot } from './values-marks.js';
import type { ValuesLay } from './values.js';

/**
 * The hover's outline on the value hovered: its bar's or its step's own extent, from the baseline to the value's end
 * (a position that has no value in view, the cell's column as every layer's).
 */
export function hoverValue(lay: ValuesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || h.position === undefined) return;
  const axis = lay.axis;
  if (h.system !== axis.system || !axis.cell(h.position)) return;
  const cx = (axis.sx0 + axis.sx1) / 2;
  let j = -1;
  for (let i = 0; i < lay.n; i++)
    if (lay.x0[i] <= cx && cx < lay.x1[i] && !Number.isNaN(lay.v[i])) {
      j = i;
      break;
    }
  if (j < 0) {
    hoverCell(lay, out, emphasis);
    return;
  }
  // `points` on a per-position step: a ring on the dot
  if (lay.encoding === 'points' && perPosition(lay.step)) {
    hoverDot(lay, out, j);
    return;
  }
  const { box } = lay;
  // `yFrame`'s and `ys`'s arithmetic written out: a hover allocates nothing (no number crosses a call boxed)
  const log = lay.yLog;
  const d0 = lay.domain ? lay.domain[0] : 0;
  const d1 = lay.domain ? lay.domain[1] : 1;
  const base = d0 <= 0 && d1 >= 0 ? 0 : d0;
  const v = lay.v[j];
  const l = Math.min(lay.lo[j], base);
  const lo = !log ? d0 : d0 > 0 ? Math.log10(d0) : Number.NEGATIVE_INFINITY;
  const hi = !log ? d1 : d1 > 0 ? Math.log10(d1) : Number.NEGATIVE_INFINITY;
  const sb = !log ? base : base > 0 ? Math.log10(base) : Number.NEGATIVE_INFINITY;
  const sv = !log ? v : v > 0 ? Math.log10(v) : Number.NEGATIVE_INFINITY;
  const sl = !log ? l : l > 0 ? Math.log10(l) : Number.NEGATIVE_INFINITY;
  const span = (hi - lo || 1) / Fill;
  const y1 = box.y + box.height;
  const bh = box.height;
  const yBase = y1 - ((Math.max(lo, Math.min(hi, sb)) - lo) / span) * bh;
  const yv = y1 - ((Math.max(lo, Math.min(hi, sv)) - lo) / span) * bh;
  const yl = y1 - ((Math.max(lo, Math.min(hi, sl)) - lo) / span) * bh;
  const top = Math.min(yv, yBase);
  const bottom = Math.max(yv, yl, yBase);
  // outside the extent by the offset, as every hover outline; a sliver of a bar still gets a visible frame
  const d = lay.hoverWidth / 2 + lay.hoverOffset;
  const c = lay.config;
  const cw = lay.x1[j] - lay.x0[j];
  const eh = Math.max(bottom - top, 1 / lay.axis.dpr);
  // the extent's squircle by the frame rule (as `hoverCell`), the outline's its plus `d`
  const r =
    cw < c.squircleMin || eh < c.squircleMin ? 0 : Math.min(eh * c.squircleRatio, c.squircleMax, cw / 2, eh / 2);
  const g = out.args;
  g[0] = lay.x0[j] - d;
  g[1] = top - d;
  g[2] = cw + 2 * d;
  g[3] = eh + 2 * d;
  g[4] = NoColor;
  g[5] = lay.tokens.colors['seq-hover'];
  g[6] = lay.hoverWidth;
  g[7] = r > 0 ? r + d : 0;
  g[8] = Round.both;
  out.rectArgs();
}
