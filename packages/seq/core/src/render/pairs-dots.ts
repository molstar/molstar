/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { AxisTransform } from './axis.js';
import { grown, NoColor, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { DotMin, pairColor, pairHit, type PairsLay, type PairsPrep, segmentDistance2 } from './pairs-lay.js';
import { type Hit, type LayoutContext } from './layer-type.js';

// The pairs layer in a grid's strip (_spec/README.md "Render", The pairs layer): a dot a pair of single positions, a
// diagonal for a stretch, at (a, b) and with `mirror` at (b, a) too; kinds at one cell side by side.

// ---- dots (a grid's strip) ----

export function layoutDots(lay: PairsLay, prep: PairsPrep, ctx: LayoutContext, a: AxisTransform, ay: AxisTransform) {
  const d = prep.data!;
  lay.mode = 'dots';
  const forward = prep.aSystem === a.system && prep.bSystem === ay.system;
  const back = prep.aSystem === ay.system && prep.bSystem === a.system;
  if (!forward && !(prep.mirror && back)) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.pairs.systems', 'The pairs are not on the grid’s two systems', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.pairs.systems'].hint,
      }),
    );
    return;
  }
  for (let k = 0; k < d.n; k++) {
    if (!prep.shown[k]) continue;
    // a stretch runs with its partner when both ends are on the same strand, against it otherwise (EC-206)
    const parallel = d.aStrand[k] === d.bStrand[k];
    if (forward) dot(lay, prep, a, ay, k, d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k], parallel);
    // a self-pair on a self grid is its own mirror: one dot
    const self = forward && d.aStart[k] === d.bStart[k] && d.aEnd[k] === d.bEnd[k];
    if (prep.mirror && back && !self)
      dot(lay, prep, a, ay, k, d.bStart[k], d.bEnd[k], d.aStart[k], d.aEnd[k], parallel);
  }
  sideBySide(lay);
}

/**
 * A pair's dot (single positions) or diagonal (a stretch: x's first with y's first when `parallel`, else with y's last),
 * kept by the triangle.
 */
function dot(
  lay: PairsLay,
  prep: PairsPrep,
  a: AxisTransform,
  ay: AxisTransform,
  k: number,
  xs: number,
  xe: number,
  ys: number,
  ye: number,
  parallel: boolean,
) {
  const t = prep.triangle;
  const cx = (xs + xe - 1) / 2;
  const cy = (ys + ye - 1) / 2;
  if ((t === 'upper' && !(cy < cx)) || (t === 'lower' && !(cy > cx))) return;
  const color = pairColor(lay, prep, k);
  if (lay.extents.any) lay.extents.add(prep.datum!.set(k));
  if (xe - xs === 1 && ye - ys === 1) {
    if (!a.cell(xs)) return;
    const ax0 = a.sx0;
    const ax1 = a.sx1;
    if (!ay.cell(ys)) return;
    const i = lay.nd++;
    lay.dots = grown(lay.dots, 4 * i + 4);
    lay.dotPair = grown(lay.dotPair, i + 1);
    lay.dotColor = grown(lay.dotColor, i + 1);
    // at least `DotMin` CSS px a side: a dot reads at any zoom
    const m = DotMin;
    const w = Math.max(Math.abs(ax1 - ax0), m);
    const h = Math.max(Math.abs(ay.sx1 - ay.sx0), m);
    lay.dots[4 * i] = (ax0 + ax1) / 2 - w / 2;
    lay.dots[4 * i + 1] = (ay.sx0 + ay.sx1) / 2 - h / 2;
    lay.dots[4 * i + 2] = w;
    lay.dots[4 * i + 3] = h;
    lay.dotPair[i] = k;
    lay.dotColor[i] = color;
    return;
  }
  // a stretch: from x's first to x's last, y's first to last (parallel) or last to first, through the cells' centres
  if (!a.cell(xs)) return;
  const x0 = (a.sx0 + a.sx1) / 2;
  if (!a.cell(xe - 1)) return;
  const x1 = (a.sx0 + a.sx1) / 2;
  if (!ay.cell(parallel ? ys : ye - 1)) return;
  const y0 = (ay.sx0 + ay.sx1) / 2;
  if (!ay.cell(parallel ? ye - 1 : ys)) return;
  const y1 = (ay.sx0 + ay.sx1) / 2;
  const i = lay.ns++;
  lay.segs = grown(lay.segs, 5 * i + 5);
  lay.segPair = grown(lay.segPair, i + 1);
  lay.segColor = grown(lay.segColor, i + 1);
  lay.segs[5 * i] = x0;
  lay.segs[5 * i + 1] = y0;
  lay.segs[5 * i + 2] = x1;
  lay.segs[5 * i + 3] = y1;
  lay.segs[5 * i + 4] = Math.max(1.5, Math.min(a.scale, ay.scale) * 0.6);
  lay.segPair[i] = k;
  lay.segColor[i] = color;
}

/** Dots at one cell share it: each a slice of its width, in order (EC-63). */
function sideBySide(lay: PairsLay) {
  if (lay.nd < 2) return;
  const d = lay.dots;
  const at = new Map<string, number[]>();
  for (let i = 0; i < lay.nd; i++) {
    const key = `${d[4 * i]}|${d[4 * i + 1]}`;
    let l = at.get(key);
    if (!l) at.set(key, (l = []));
    l.push(i);
  }
  for (const l of at.values()) {
    if (l.length < 2) continue;
    const w = d[4 * l[0] + 2] / l.length;
    for (let j = 0; j < l.length; j++) {
      d[4 * l[j]] += j * w;
      d[4 * l[j] + 2] = w;
    }
  }
}

export function drawDots(lay: PairsLay, out: PassWriter) {
  const g = out.args;
  // outlined in the text's ink (dark on light, light on dark): a dot reads over any color of the cells under it
  const line = lay.tokens!.colors['msui-text'];
  for (let i = 0; i < lay.nd; i++) {
    g[0] = lay.dots[4 * i];
    g[1] = lay.dots[4 * i + 1];
    g[2] = lay.dots[4 * i + 2];
    g[3] = lay.dots[4 * i + 3];
    g[4] = lay.dotColor[i];
    g[5] = line;
    g[6] = 1;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  const p = lay.seg;
  for (let i = 0; i < lay.ns; i++) {
    p[0] = lay.segs[5 * i];
    p[1] = lay.segs[5 * i + 1];
    p[2] = lay.segs[5 * i + 2];
    p[3] = lay.segs[5 * i + 3];
    g[4] = NoColor;
    g[5] = lay.segColor[i];
    g[6] = lay.segs[5 * i + 4];
    g[7] = 0;
    out.pathArgs(p, 2, false);
  }
}

/** The dot under a point (its box), else the diagonal within `pairs.hitPx` of it. */
export function hitDots(lay: PairsLay, at: Float64Array): Hit | undefined {
  const x = at[0];
  const y = at[1];
  const d = lay.dots;
  for (let i = lay.nd - 1; i >= 0; i--)
    if (x >= d[4 * i] && x < d[4 * i] + d[4 * i + 2] && y >= d[4 * i + 1] && y < d[4 * i + 1] + d[4 * i + 3])
      return pairHit(lay, lay.dotPair[i]);
  const r = lay.config!.pairs.hitPx;
  const s = lay.segs;
  for (let i = lay.ns - 1; i >= 0; i--) {
    const w = Math.max(r, s[5 * i + 4] / 2);
    if (segmentDistance2(s[5 * i], s[5 * i + 1], s[5 * i + 2], s[5 * i + 3], x, y) <= w * w)
      return pairHit(lay, lay.segPair[i]);
  }
  return undefined;
}
