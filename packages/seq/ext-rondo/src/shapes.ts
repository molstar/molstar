/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type LayoutContext, TakenLabels } from '@molstar/seqx-core/extension';
import { ColorMode } from './colour.js';
import type { RondoConfig } from './config.js';
import { addArc, addCubic, grow, type RibbonFrame, ribbonCorners } from './geometry.js';
import type { RondoLay, RondoPrep } from './lay.js';

// The ribbons' shapes (_spec/README.md, Stems, Ribbons, Stubs): the laid-out stems, ribbons, stubs and fans each as a
// polygon in draw order, the hits they stand for, and the grid over the hole a hit reads (one cover: the closed hit
// equals the drawn run).

/** A polygon's kind. */
export const Kind = { trumpet: 0, ribbon: 1, stub: 2, fan: 3 } as const;

/** The labels taken when a layer is laid out alone (a section's layers share its own: `LayoutContext.taken`). */
const Alone = new TakenLabels();
let ByCount = new Int32Array(16);

/**
 * The stubs' labels (`chr7 · 37` at the tip, centred): the strongest first, each only where its box (padded 2 px) meets
 * no label taken and its tip is `stubLabelPx` from every taken tip; the rest unlabelled (their partner in the card).
 * The labels taken are the section's, across its layers (two kinds of links over one ring: their stubs at one tip).
 */
function stubLabels(lay: RondoLay, cfg: RondoConfig, ctx: LayoutContext, tipD: number) {
  const f = lay.frame;
  const font = lay.font;
  const size = Number.parseFloat(font) || 12;
  if (ByCount.length < lay.nb) ByCount = new Int32Array(lay.nb);
  let n = 0;
  // a folded stub unlabelled (every stem has one: its text would repeat round the ring); its count in its card, its
  // partners unfolded on hover
  for (let s = 0; s < lay.nb; s++) if (lay.stubPoly[s] >= 0 && lay.stubPartner[s] >= 0) ByCount[n++] = s;
  const order = ByCount.subarray(0, n).sort((a, b) => lay.stubCount[b] - lay.stubCount[a] || a - b);
  const taken = ctx.taken ?? Alone;
  if (!ctx.taken) taken.reset();
  for (const s of order) {
    const partner = lay.stubPartner[s];
    const text = `${lay.names[partner]} · ${Math.round(lay.stubLinks[s]).toLocaleString('en-US')}`;
    f.point((lay.stubT0[s] + lay.stubT1[s]) / 2, tipD + 4);
    const x = f.px;
    const y = f.py;
    const w = ctx.text.measure(font, text);
    // the box about the baseline (`y + 4`): the cap height above, a descent below
    const x0 = x - w / 2 - 2;
    const x1 = x + w / 2 + 2;
    const y0 = y + 4 - size - 2;
    const y1 = y + 4 + size * 0.3 + 2;
    if (!taken.free(x0, y0, x1, y1, x, y, cfg.stubLabelPx)) continue;
    taken.take(x0, y0, x1, y1, x, y);
    lay.stubLabel[s] = text;
    lay.stubLabelX[s] = x;
    lay.stubLabelY[s] = y;
    lay.stubLabelW[s] = w;
  }
}

/**
 * A thin polygon's slop: a ribbon or stub under `2 × hitSlopPx` wide is hit within `hitSlopPx` of its centre line, so
 * past its outline by what its half width leaves (a hairline ribbon a pointer can stay on).
 */
function slop(lay: RondoLay, i: number, widthPx: number, hitSlopPx: number) {
  lay.polySlop[i] = Math.max(0, hitSlopPx - widthPx / 2);
}

/** A polygon made: its kind, what it stands for, its colour. */
function made(lay: RondoLay, kind: number, ref: number, color: number): number {
  const i = lay.polys.end();
  if (i >= lay.polyKind.length) {
    lay.polyKind = grow(lay.polyKind, i + 1);
    lay.polyRef = grow(lay.polyRef, i + 1);
    lay.polyColor = grow(lay.polyColor, i + 1);
    lay.polySlop = grow(lay.polySlop, i + 1);
  }
  lay.polySlop[i] = 0;
  lay.polyKind[i] = kind;
  lay.polyRef[i] = ref;
  lay.polyColor[i] = color;
  return i;
}

const alpha = (c: number, a: number) => ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;

/**
 * A system's hue (its index in the systems' order: the domain taken before `include`); the layer's one colour when it
 * has one (a fixed colour, a field's): stems, stubs and fans take it.
 */
export function hueOf(lay: RondoLay, sys: number): number {
  return lay.one >= 0 ? lay.one : lay.hues[sys % lay.hues.length];
}

const C4 = new Float64Array(4);

/** The polygons in draw order: trumpets, ribbons widest first (the narrow on top, hittable), stubs, fans. */
export function polygons(lay: RondoLay, prep: RondoPrep, cfg: RondoConfig, ctx: LayoutContext) {
  const f = lay.frame;
  const P = lay.polys;
  const seg = Math.max(2, cfg.flatten);
  const dn = cfg.stemDepth * f.room;
  // a fixed colour (a kind over others) opaque, so it reads as one colour over any hue; else `ribbons.fillOpacity`
  const op = lay.mode === ColorMode.fixed ? lay.opacity : lay.fillOpacity * lay.opacity;
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
    // on one chromosome or a pair the ribbons and stubs carry the information: the trumpets lighter there
    const stemOp = lay.focus ? op * cfg.focusStems : op;
    lay.stemPoly[s] = made(lay, Kind.trumpet, s, alpha(hueOf(lay, lay.stemSystem[s]), stemOp));
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
    const color = lay.byRow ? lay.ribColor[r] : hueOf(lay, lay.stemSystem[a]);
    lay.ribPoly[r] = made(lay, Kind.ribbon, r, alpha(color, op));
    slop(lay, lay.ribPoly[r], lay.k * lay.ribCount[r], cfg.hitSlopPx);
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
    // a folded stub (partner −1) in the muted text's colour (no one partner's hue), or the layer's one colour
    const hue = lay.stubPartner[s] < 0 ? (lay.one >= 0 ? lay.one : lay.muted) : hueOf(lay, lay.stubPartner[s]);
    lay.stubPoly[s] = made(lay, Kind.stub, s, alpha(hue, op));
    slop(lay, lay.stubPoly[s], lay.k * lay.stubCount[s], cfg.hitSlopPx);
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
      // a field colour by link: the stem pair's ribbons told apart by colour
      if (lay.byRow) key += `|${(lay.ribColor[ref] >>> 8).toString(16).padStart(6, '0')}`;
      count = lay.ribCount[ref];
      lay.ribKey[ref] = key;
    } else if (kind === Kind.stub) {
      stem = lay.stubStem[ref];
      const partner = lay.stubPartner[ref];
      key = `stub:${lay.systems[lay.stemSystem[stem]]}:${lay.stemBin[stem]}|${partner < 0 ? 'other' : lay.systems[partner]}`;
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
    const e = lay.polySlop[i];
    x0 = Math.min(x0, P.box[4 * i] - e);
    y0 = Math.min(y0, P.box[4 * i + 1] - e);
    x1 = Math.max(x1, P.box[4 * i + 2] + e);
    y1 = Math.max(y1, P.box[4 * i + 3] + e);
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
    const e = lay.polySlop[i];
    const cx0 = Math.floor((P.box[4 * i] - e - x0) / c);
    const cy0 = Math.floor((P.box[4 * i + 1] - e - y0) / c);
    const cx1 = Math.floor((P.box[4 * i + 2] + e - x0) / c);
    const cy1 = Math.floor((P.box[4 * i + 3] + e - y0) / c);
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

/** A point in and its distance out (`Polygons.distance`). */
const Dist = new Float64Array(3);

/**
 * The polygon under a page point (`at`), or −1 (through the grid: a few candidates; nothing allocated), topmost drawn first:
 * a thin ribbon or stub covers its outline and `rondo.hitSlopPx` about its centre line, so a hairline drawn over a
 * wide ribbon is hit along it; where thin covers overlap, the nearest centre line; the first wide cover below ends it.
 */
export function polygonAt(lay: RondoLay, at: ArrayLike<number>): number {
  if (!lay.gridCols) return -1;
  // the point read here, never passed: a double passed to a call not inlined is boxed
  const x = at[0];
  const y = at[1];
  const cx = Math.floor((x - lay.gridX) / lay.cell);
  const cy = Math.floor((y - lay.gridY) / lay.cell);
  if (cx < 0 || cy < 0 || cx >= lay.gridCols || cy >= lay.gridRows) return -1;
  const k = cy * lay.gridCols + cx;
  const P = lay.polys;
  const b = P.box;
  let near = -1;
  let best = Infinity;
  for (let j = lay.gridStart[k]; j < lay.gridStart[k + 1]; j++) {
    const i = lay.gridItems[j];
    const e = lay.polySlop[i];
    if (e <= 0) {
      if (P.contains(i, x, y)) return near >= 0 ? near : i;
      continue;
    }
    if (x < b[4 * i] - e || x > b[4 * i + 2] + e || y < b[4 * i + 1] - e || y > b[4 * i + 3] + e) continue;
    // the centre line's distance less the slop (inside: the outline's distance taken off)
    Dist[0] = x;
    Dist[1] = y;
    P.distance(i, Dist);
    const d = P.contains(i, x, y) ? -Dist[2] - e : Dist[2] - e;
    if (d <= 0 && d < best) {
      best = d;
      near = i;
    }
  }
  return near;
}
