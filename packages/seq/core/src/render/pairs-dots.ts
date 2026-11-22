/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { AxisTransform } from './axis.js';
import { grown, NoColor, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { aOf, bOf, DotMin, pairColor, pairHit, type PairsLay, type PairsPrep, segmentDistance2 } from './pairs-lay.js';
import { type Hit, type LayoutContext } from './layer-type.js';
import { withAlpha } from './tokens.js';
import { keeps, StripKeep, stripKeep, type Triangle } from './triangle.js';
import type { PairsData } from '@molstar/seqx-core/resolve';

// The pairs layer in a grid's strip (_spec/README.md "Render", The pairs layer): a dot a pair of single positions, a
// diagonal for a stretch, at (a, b) and with `mirror` at (b, a) too; kinds at one cell side by side.

// ---- dots (a grid's strip) ----

/**
 * A strip's pairs (module-private, by lay): each pair on its two systems (`dir` 1) or, with `mirror`, the other way (2,
 * both: 3), and the bins' greatest weight among them; made when the data, the strip's systems or `mirror` change, so a
 * warm frame scans no pair it does not draw.
 */
interface StripPairs {
  data: PairsData | undefined;
  shown: Uint8Array | undefined;
  a: string;
  b: string;
  mirror: boolean;
  n: number;
  pair: Int32Array;
  dir: Uint8Array;
  most: number;
}

const OnStrip = new WeakMap<PairsLay, StripPairs>();

function stripPairs(lay: PairsLay, prep: PairsPrep, a: AxisTransform, ay: AxisTransform): StripPairs {
  let sp = OnStrip.get(lay);
  if (!sp) {
    sp = {
      data: undefined,
      shown: undefined,
      a: '',
      b: '',
      mirror: false,
      n: 0,
      pair: new Int32Array(16),
      dir: new Uint8Array(16),
      most: 0,
    };
    OnStrip.set(lay, sp);
  }
  const d = prep.data!;
  if (sp.data === d && sp.shown === prep.shown && sp.a === a.system && sp.b === ay.system && sp.mirror === prep.mirror)
    return sp;
  sp.data = d;
  sp.shown = prep.shown;
  sp.a = a.system;
  sp.b = ay.system;
  sp.mirror = prep.mirror;
  sp.n = 0;
  sp.most = 0;
  const units = prep.endUnits !== undefined;
  const w = prep.weight;
  for (let k = 0; k < d.n; k++) {
    if (!prep.shown[k]) continue;
    const as = units ? aOf(prep, k) : prep.aSystem;
    const bs = units ? bOf(prep, k) : prep.bSystem;
    const dir =
      (as === a.system && bs === ay.system ? 1 : 0) | (prep.mirror && as === ay.system && bs === a.system ? 2 : 0);
    if (!dir) continue;
    const i = sp.n++;
    sp.pair = grown(sp.pair, i + 1);
    sp.dir = grown(sp.dir, i + 1);
    sp.pair[i] = k;
    sp.dir[i] = dir;
    if (w && w[k] > sp.most) sp.most = w[k];
  }
  return sp;
}

export function layoutDots(lay: PairsLay, prep: PairsPrep, ctx: LayoutContext, a: AxisTransform, ay: AxisTransform) {
  const d = prep.data!;
  lay.mode = 'dots';
  // each pair's ends on their own systems (an assembly's sequences, a complex's units): over blocks a strip is a block
  // pair, and draws the pairs of its two systems alone (a dot plot of two genomes, EC-209)
  const units = prep.endUnits !== undefined;
  const forward = prep.aSystem === a.system && prep.bSystem === ay.system;
  const back = prep.aSystem === ay.system && prep.bSystem === a.system;
  if (!units && !forward && !(prep.mirror && back)) {
    if (!ctx.section.blocks)
      ctx.report(
        issueOf(RenderIssues, 'seq.render.pairs.systems', 'The pairs are not on the grid’s two systems', {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.pairs.systems'].hint,
        }),
      );
    return;
  }
  // the triangle on this strip (`stripKeep`): a cell by its display ranks, never across two systems or two sections
  const keep = stripKeep(prep.triangle, ctx);
  if (keep === StripKeep.None) return;
  Keep.t = keep === StripKeep.Cells ? prep.triangle : 'both';
  Keep.sx = a.dir;
  Keep.sy = ay.dir;
  const sp = stripPairs(lay, prep, a, ay);
  const w = prep.weight;
  for (let q = 0; q < sp.n; q++) {
    const k = sp.pair[q];
    const fw = (sp.dir[q] & 1) !== 0;
    const bw = (sp.dir[q] & 2) !== 0;
    // a self-pair on a self grid is its own mirror: one mark
    const self = fw && d.aStart[k] === d.bStart[k] && d.aEnd[k] === d.bEnd[k];
    if (w) {
      // bins (a source's pairs zoomed out): a rect a bin pair at its bins' size, its opacity by count, the square root
      // over the bins on this strip (as a bridge's, `ribbons.binOpacity`); with `mirror` at (b, a) too
      if (fw) binRect(lay, prep, ctx, a, ay, k, sp.most, d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k]);
      if (bw && !self) binRect(lay, prep, ctx, a, ay, k, sp.most, d.bStart[k], d.bEnd[k], d.aStart[k], d.aEnd[k]);
      continue;
    }
    // a stretch runs with its partner when both ends are on the same strand, against it otherwise (EC-206)
    const parallel = d.aStrand[k] === d.bStrand[k];
    if (fw) dot(lay, prep, a, ay, k, d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k], parallel);
    if (bw && !self) dot(lay, prep, a, ay, k, d.bStart[k], d.bEnd[k], d.aStart[k], d.aEnd[k], parallel);
  }
  sideBySide(lay);
}

/** The strip's triangle and the axes' directions (display ranks: a reversed axis's positions negated), set a layout. */
const Keep = { t: 'both' as Triangle, sx: 1, sy: 1 };

/** Whether a mark centred at positions (cx, cy) is kept by the strip's triangle, above or below the drawn diagonal. */
function kept(cx: number, cy: number): boolean {
  return keeps(Keep.t, Keep.sx * cx, Keep.sy * cy);
}

/** `[s, e)`'s extent in view on an axis into `out[o]`, `out[o + 1]` (its first piece in view, clipped); false if none. */
function extentIn(ax: AxisTransform, s: number, e: number, out: Float64Array, o: number): boolean {
  ax.piecesIn(s, e);
  for (let k = ax.ka; k <= ax.kb; k++)
    if (ax.spanIn(k, s, e)) {
      out[o] = ax.sx0;
      out[o + 1] = ax.sx1;
      return true;
    }
  return false;
}

const BinBox = new Float64Array(4);

/** A bin pair as a rect of its two bins (at least `DotMin` CSS px a side), its opacity by its count's share. */
function binRect(
  lay: PairsLay,
  prep: PairsPrep,
  ctx: LayoutContext,
  a: AxisTransform,
  ay: AxisTransform,
  k: number,
  most: number,
  xs: number,
  xe: number,
  ys: number,
  ye: number,
) {
  if (!kept((xs + xe - 1) / 2, (ys + ye - 1) / 2)) return;
  if (!extentIn(a, xs, xe, BinBox, 0) || !extentIn(ay, ys, ye, BinBox, 2)) return;
  const color = pairColor(lay, prep, k);
  if (lay.extents.any) lay.extents.add(prep.datum!.set(k));
  const { min, max } = ctx.config.ribbons.binOpacity;
  const share = most > 0 ? Math.sqrt(prep.weight![k] / most) : 1;
  const i = lay.nd++;
  lay.dots = grown(lay.dots, 4 * i + 4);
  lay.dotPair = grown(lay.dotPair, i + 1);
  lay.dotColor = grown(lay.dotColor, i + 1);
  const w = Math.max(BinBox[1] - BinBox[0], DotMin);
  const h = Math.max(BinBox[3] - BinBox[2], DotMin);
  lay.dots[4 * i] = (BinBox[0] + BinBox[1]) / 2 - w / 2;
  lay.dots[4 * i + 1] = (BinBox[2] + BinBox[3]) / 2 - h / 2;
  lay.dots[4 * i + 2] = w;
  lay.dots[4 * i + 3] = h;
  lay.dotPair[i] = k;
  lay.dotColor[i] = withAlpha(color, min + (max - min) * share);
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
  if (!kept((xs + xe - 1) / 2, (ys + ye - 1) / 2)) return;
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
