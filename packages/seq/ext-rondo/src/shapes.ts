/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LayoutContext } from '@molstar/seqx-core/extension';
import type { RondoConfig } from './config.js';
import { addArc, addCubic, grow, type RibbonFrame, ribbonCorners } from './geometry.js';
import type { RondoLay, RondoPrep } from './lay.js';

// The ribbons' shapes (_spec/README.md, Stems, Ribbons, Stubs): the laid-out stems, ribbons, stubs and fans each as a
// polygon in draw order, the hits they stand for, and the grid over the hole a hit reads (one cover: the closed hit
// equals the drawn run).

/** A polygon's kind. */
export const Kind = { trumpet: 0, ribbon: 1, stub: 2, fan: 3 } as const;

/** The labels taken so far: their box and tip (x0, y0, x1, y1, x, y a label), reused. */
let Boxes = new Float64Array(64);
let ByCount = new Int32Array(16);

/**
 * The stubs' labels (`chr7 · 37` at the tip, centred): the strongest first, each only where its box (padded 2 px) meets
 * no label taken and its tip is `stubLabelPx` from every taken tip; the rest unlabelled (their partner in the card).
 */
function stubLabels(lay: RondoLay, cfg: RondoConfig, ctx: LayoutContext, tipD: number) {
  const f = lay.frame;
  const font = lay.font;
  const size = Number.parseFloat(font) || 12;
  if (ByCount.length < lay.nb) ByCount = new Int32Array(lay.nb);
  let n = 0;
  for (let s = 0; s < lay.nb; s++) if (lay.stubPoly[s] >= 0) ByCount[n++] = s;
  const order = ByCount.subarray(0, n).sort((a, b) => lay.stubCount[b] - lay.stubCount[a] || a - b);
  let taken = 0;
  for (const s of order) {
    const text = `${lay.names[lay.stubPartner[s]]} · ${Math.round(lay.stubCount[s]).toLocaleString('en-US')}`;
    f.point((lay.stubT0[s] + lay.stubT1[s]) / 2, tipD + 4);
    const x = f.px;
    const y = f.py;
    const w = ctx.text.measure(font, text);
    // the box about the baseline (`y + 4`): the cap height above, a descent below
    const x0 = x - w / 2 - 2;
    const x1 = x + w / 2 + 2;
    const y0 = y + 4 - size - 2;
    const y1 = y + 4 + size * 0.3 + 2;
    let free = true;
    for (let j = 0; j < taken && free; j++) {
      const o = 6 * j;
      const overlaps = x0 < Boxes[o + 2] && x1 > Boxes[o] && y0 < Boxes[o + 3] && y1 > Boxes[o + 1];
      if (overlaps || Math.hypot(Boxes[o + 4] - x, Boxes[o + 5] - y) < cfg.stubLabelPx) free = false;
    }
    if (!free) continue;
    if (Boxes.length < 6 * (taken + 1)) Boxes = grow(Boxes, 12 * (taken + 1));
    const o = 6 * taken;
    Boxes[o] = x0;
    Boxes[o + 1] = y0;
    Boxes[o + 2] = x1;
    Boxes[o + 3] = y1;
    Boxes[o + 4] = x;
    Boxes[o + 5] = y;
    taken++;
    lay.stubLabel[s] = text;
    lay.stubLabelX[s] = x;
    lay.stubLabelY[s] = y;
    lay.stubLabelW[s] = w;
  }
}

/** A polygon made: its kind, what it stands for, its colour. */
function made(lay: RondoLay, kind: number, ref: number, color: number): number {
  const i = lay.polys.end();
  if (i >= lay.polyKind.length) {
    lay.polyKind = grow(lay.polyKind, i + 1);
    lay.polyRef = grow(lay.polyRef, i + 1);
    lay.polyColor = grow(lay.polyColor, i + 1);
  }
  lay.polyKind[i] = kind;
  lay.polyRef[i] = ref;
  lay.polyColor[i] = color;
  return i;
}

const alpha = (c: number, a: number) => ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;

/** A system's hue (its index in the systems' order: the domain taken before `include`). */
export function hueOf(lay: RondoLay, sys: number): number {
  return lay.hues[sys % lay.hues.length];
}

const C4 = new Float64Array(4);

/** The polygons in draw order: trumpets, ribbons widest first (the narrow on top, hittable), stubs, fans. */
export function polygons(lay: RondoLay, prep: RondoPrep, cfg: RondoConfig, ctx: LayoutContext) {
  const f = lay.frame;
  const P = lay.polys;
  const seg = Math.max(2, cfg.flatten);
  const dn = cfg.stemDepth * f.room;
  const op = lay.fillOpacity;
  // trumpets (a stem's mouth to its neck: two cubics flaring), stems: `none` draws none
  for (let s = 0; s < lay.ns; s++) {
    if (!lay.stemPlaced[s] || lay.stemCount[s] <= 0) continue;
    const t0 = lay.stemT0[s];
    const t1 = lay.stemT1[s];
    const mid = (t0 + t1) / 2;
    const n0 = mid - lay.stemNeck[s];
    const n1 = mid + lay.stemNeck[s];
    if (prep.stems === 'none') {
      lay.stemPoly[s] = -1;
      continue;
    }
    P.begin();
    f.point(t0, 0);
    P.add(f.px, f.py);
    addArc(f, P, t0, t1, 0);
    addCubic(f, P, seg, t1, 0, t1, dn / 2, n1, dn / 2, n1, dn);
    addArc(f, P, n1, n0, dn);
    addCubic(f, P, seg, n0, dn, n0, dn / 2, t0, dn / 2, t0, 0);
    lay.stemPoly[s] = made(lay, Kind.trumpet, s, alpha(hueOf(lay, lay.stemSystem[s]), op));
  }
  // ribbons widest first
  const ribs: number[] = [];
  for (let r = 0; r < lay.nr; r++) if (lay.ribDrawn[r]) ribs.push(r);
  ribs.sort((x, y) => lay.ribCount[y] - lay.ribCount[x] || x - y);
  const bundle = prep.bundle;
  for (const r of ribs) {
    const a = lay.ribA[r];
    const b = lay.ribB[r];
    // uncrossed: A's low edge to B's high, A's high to B's low; crossed: low to low, high to high
    ribbonCorners(
      lay.ribA0[r],
      lay.ribA1[r],
      lay.stemDir[a],
      lay.ribB0[r],
      lay.ribB1[r],
      lay.stemDir[b],
      lay.ribStrand[r],
      C4,
    );
    const dc = controlDepth(f, dn, bundle, C4[0], C4[1]);
    const dc2 = controlDepth(f, dn, bundle, C4[2], C4[3]);
    P.begin();
    f.point(C4[0], dn);
    P.add(f.px, f.py);
    addCubic(f, P, seg, C4[0], dn, C4[0], dc, C4[1], dc, C4[1], dn);
    f.point(C4[2], dn);
    P.add(f.px, f.py);
    addCubic(f, P, seg, C4[2], dn, C4[2], dc2, C4[3], dc2, C4[3], dn);
    lay.ribPoly[r] = made(lay, Kind.ribbon, r, alpha(hueOf(lay, lay.stemSystem[a]), op));
  }
  // stubs: a taper from the neck toward the centre, in the partner's hue, labelled at its tip where it has room
  const tipD = dn + cfg.stubShare * f.room;
  for (let s = 0; s < lay.nb; s++) {
    lay.stubPoly[s] = -1;
    lay.stubLabel[s] = '';
    if (!lay.stubDrawn[s] || prep.stubs === 'none') continue;
    const t0 = lay.stubT0[s];
    const t1 = lay.stubT1[s];
    const mid = (t0 + t1) / 2;
    P.begin();
    f.point(t0, dn);
    P.add(f.px, f.py);
    f.point(t1, dn);
    P.add(f.px, f.py);
    f.point(mid, tipD);
    P.add(f.px, f.py);
    lay.stubPoly[s] = made(lay, Kind.stub, s, alpha(hueOf(lay, lay.stubPartner[s]), op));
  }
  stubLabels(lay, cfg, ctx, tipD);
  // fans: a stem's weak links, one taper of their summed width, last in its neck
  for (let s = 0; s < lay.ns; s++) {
    lay.stemFanPoly[s] = -1;
    if (!lay.stemPlaced[s] || lay.stemFan[s] <= 0) continue;
    const t0 = lay.stemCut[s];
    const t1 = t0 + Math.max(lay.k * lay.stemFan[s], 0.5 / lay.dpr) * f.perPx(dn);
    P.begin();
    f.point(t0, dn);
    P.add(f.px, f.py);
    f.point(t1, dn);
    P.add(f.px, f.py);
    f.point((t0 + t1) / 2, dn + (cfg.stubShare / 2) * f.room);
    P.add(f.px, f.py);
    lay.stemFanPoly[s] = made(lay, Kind.fan, s, alpha(hueOf(lay, lay.stemSystem[s]), op * 0.6));
  }
  // the hole (a ring's) or the band (a row's): the wash's polygon, never hit
  P.begin();
  if (f.polar) addArc(f, P, 0, 2 * Math.PI, dn);
  else {
    P.add(f.x0, f.r);
    P.add(f.x1, f.r);
    P.add(f.x1, f.r + f.room);
    P.add(f.x0, f.r + f.room);
  }
  lay.hole = made(lay, -1, -1, 0);
  hitsOf(lay);
}

/** A ribbon's controls' depth: `(1 − bundle)` of the neck's radius from the centre on a ring; on a row by its span. */
function controlDepth(f: RibbonFrame, dn: number, bundle: number, ta: number, tb: number): number {
  if (f.polar) return f.r - (1 - bundle) * (f.r - dn);
  const reach = Math.min(1, Math.abs(tb - ta) / Math.max(1, f.x1 - f.x0));
  return dn + (f.room - dn) * Math.max(0.15, reach) * (0.5 + bundle / 2);
}

/** Each polygon's hit and its key's index (made at layout: a hit allocates nothing). */
function hitsOf(lay: RondoLay) {
  const owner = lay.owner;
  const section = lay.sectionId;
  const P = lay.polys;
  lay.order = grow(lay.order, P.n);
  for (let i = 0; i < P.n; i++) {
    lay.order[i] = i;
    const kind = lay.polyKind[i];
    const ref = lay.polyRef[i];
    let key = '';
    let stem = -1;
    let count = 0;
    if (kind === Kind.trumpet || kind === Kind.fan) {
      stem = ref;
      key = `stem:${lay.systems[lay.stemSystem[ref]]}:${lay.stemBin[ref]}`;
      count = kind === Kind.fan ? lay.stemFan[ref] : lay.stemCount[ref];
      lay.stemKey[ref] = key;
    } else if (kind === Kind.ribbon) {
      stem = lay.ribA[ref];
      const b = lay.ribB[ref];
      key = `ribbon:${lay.systems[lay.stemSystem[stem]]}:${lay.stemBin[stem]}|${lay.systems[lay.stemSystem[b]]}:${lay.stemBin[b]}|${lay.ribStrand[ref]}`;
      count = lay.ribCount[ref];
      lay.ribKey[ref] = key;
    } else if (kind === Kind.stub) {
      stem = lay.stubStem[ref];
      key = `stub:${lay.systems[lay.stemSystem[stem]]}:${lay.stemBin[stem]}|${lay.systems[lay.stubPartner[ref]]}`;
      count = lay.stubCount[ref];
      lay.stubKey[ref] = key;
    } else {
      lay.hits.push(undefined as never);
      continue;
    }
    // a fan is hit as its stem (its key the stem's): the first made (the trumpet) is the key's
    if (!lay.byKey.has(key)) lay.byKey.set(key, i);
    const sys = lay.systems[lay.stemSystem[stem]];
    lay.hits.push({
      layer: lay.key,
      section,
      system: sys,
      element: { owner, key },
      position: (lay.stemStart[stem] + lay.stemEnd[stem]) / 2,
      bin: lay.stemLevel,
      count,
    });
  }
}

/** A uniform grid over the polygons' boxes, each cell listing the polygons over it (`rondo.hitCellPx`). */
export function hitGrid(lay: RondoLay, cellPx: number) {
  const P = lay.polys;
  let x0 = Infinity;
  let y0 = Infinity;
  let x1 = -Infinity;
  let y1 = -Infinity;
  for (let i = 0; i < P.n; i++) {
    if (i === lay.hole) continue;
    x0 = Math.min(x0, P.box[4 * i]);
    y0 = Math.min(y0, P.box[4 * i + 1]);
    x1 = Math.max(x1, P.box[4 * i + 2]);
    y1 = Math.max(y1, P.box[4 * i + 3]);
  }
  const c = Math.max(1, cellPx);
  lay.cell = c;
  if (!(x1 >= x0)) {
    lay.gridCols = lay.gridRows = 0;
    return;
  }
  lay.gridX = x0;
  lay.gridY = y0;
  lay.gridCols = Math.max(1, Math.ceil((x1 - x0) / c) + 1);
  lay.gridRows = Math.max(1, Math.ceil((y1 - y0) / c) + 1);
  const cells = lay.gridCols * lay.gridRows;
  lay.gridStart = grow(lay.gridStart, cells + 1);
  lay.gridStart.fill(0, 0, cells + 1);
  const span = (i: number, add: (cell: number) => void) => {
    const cx0 = Math.floor((P.box[4 * i] - x0) / c);
    const cy0 = Math.floor((P.box[4 * i + 1] - y0) / c);
    const cx1 = Math.floor((P.box[4 * i + 2] - x0) / c);
    const cy1 = Math.floor((P.box[4 * i + 3] - y0) / c);
    for (let y = cy0; y <= cy1; y++) for (let x = cx0; x <= cx1; x++) add(y * lay.gridCols + x);
  };
  for (let i = 0; i < P.n; i++) if (i !== lay.hole) span(i, (k) => lay.gridStart[k + 1]++);
  for (let k = 0; k < cells; k++) lay.gridStart[k + 1] += lay.gridStart[k];
  lay.gridItems = grow(lay.gridItems, lay.gridStart[cells]);
  const fill = new Int32Array(cells);
  // topmost first: the last drawn first in each cell
  for (let i = P.n - 1; i >= 0; i--)
    if (i !== lay.hole) span(i, (k) => (lay.gridItems[lay.gridStart[k] + fill[k]++] = i));
}

/** The polygon under a page point, topmost first, or −1 (through the grid: a few candidates; nothing allocated). */
export function polygonAt(lay: RondoLay, x: number, y: number): number {
  if (!lay.gridCols) return -1;
  const cx = Math.floor((x - lay.gridX) / lay.cell);
  const cy = Math.floor((y - lay.gridY) / lay.cell);
  if (cx < 0 || cy < 0 || cx >= lay.gridCols || cy >= lay.gridRows) return -1;
  const k = cy * lay.gridCols + cx;
  for (let j = lay.gridStart[k]; j < lay.gridStart[k + 1]; j++) {
    const i = lay.gridItems[j];
    if (lay.polys.contains(i, x, y)) return i;
  }
  return -1;
}
