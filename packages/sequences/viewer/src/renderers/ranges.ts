/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `ranges` renderer (seqspec:viewer#extension-points): packed lanes by item over its whole extent, a hairline
// across each gap between an item's fragments, brackets, labels cut with an ellipsis, boundaries with glyphs; wrapped,
// underlines on the lanes of the line; more items in view than device columns (one mode per layer per frame, with
// hysteresis): density per column and color (the shared `drawDensity`, coverage), never one color's bars over
// another's. Every item counts: one narrower than a pixel is a device pixel wide. `tiny` merges the lanes into one
// band. A label is centered on its piece's whole extent and moves with it, sticking to the data area's edge only where
// it would leave the
// piece's visible part (seqspec:sources#levels-and-the-zoom-gate, the sticky label). Over its track's letters a
// boundary keeps off their text: a small triangle pinned to the row's top edge and a hairline down the boundary, or the
// hairline alone (seqspec:viewer#the-render-backend). The overlay outlines the hovered item's pieces
// (seqspec:interaction#hover).
import { Annotations } from "@molstar/seq-core/model";
import {
  type Axis,
  AxisGeometry,
  type AxisView,
  bar,
  barBand,
  boundaryEdge,
  type DrawView,
  densityBand,
  densityColumns,
  densityMode,
  densityShown,
  drawLineLanes,
  drawStreamedItems,
  type Features,
  featureAt,
  featureGlyph,
  featureHoverKey,
  featureInk,
  featureKey,
  featureSize,
  fitIn,
  fitted,
  GlyphBatch,
  type Hit,
  HoverOutline,
  hitItem,
  hitLayer,
  hitLineLanes,
  hitStreamed,
  hitWindow,
  hoverEmphasis,
  hoveredAt,
  hoveredItem,
  itemsInView,
  kernelEntry,
  keyAlpha,
  Lod,
  laneHeight,
  laneInset,
  lanePad,
  laneRows,
  NO_SIZE,
  outlineLineLanes,
  outlineStreamed,
  owned,
  PathBatch,
  PathMode,
  PathVerb,
  placeRun,
  RectBatch,
  SkinToken,
  segWindow,
  shownHull,
  sizePx,
  sizeTheme,
  TextAlign,
  TextBaseline,
  TextBatch,
  type TextManagerApi,
  ThemeLocation,
  visibleRange,
} from "../draw/api.js";

type Params = Readonly<Record<string, unknown>>;

const rects = new RectBatch(256);
const glyphs = new GlyphBatch(32);
const names = new TextBatch(32);
const heads = new PathBatch(32);
const range = new Int32Array(2);
const hull = new Int32Array(2);
const win = new Int32Array(2);
const piece = new Float64Array(64);
const fit = new Float64Array(2);
/** The data area's left and right edges (device px), set by the draw. */
const edge = new Int32Array(2);
/** A labeled piece's whole extent (device px), and the shown runs it is found in. */
const ext = new Int32Array(2);
/** A label's half width and center (device px). */
const at = new Float64Array(2);
let shown = new Int32Array(16);

/** Per prepared state, its items' directions by `direction_field` (1, -1; 0: none). Edge-rate: on first draw. */
const directions = new WeakMap<Features, { field: string; d: Int8Array }>();
const FORWARD = new Set(["+", "1", "forward", "fwd"]);
const REVERSE = new Set(["-", "-1", "reverse", "rev"]);
function directionsOf(v: DrawView, f: Features, p: Params): Int8Array | null {
  const field = p.direction_field;
  const set = v.axis?.set;
  if (typeof field !== "string" || !set) return null;
  const known = directions.get(f);
  if (known && known.field === field) return known.d;
  const d = new Int8Array(f.annotation.length);
  for (let k = 0; k < d.length; k++) {
    const a = set.annotations[f.annotation[k]];
    const x = a && f.row[k] >= 0 ? String(Annotations.rowValue(a, set, f.row[k], field) ?? "").toLowerCase() : "";
    d[k] = FORWARD.has(x) ? 1 : REVERSE.has(x) ? -1 : 0;
  }
  directions.set(f, { field, d });
  return d;
}
/** The items' directions while `items` draws (`direction_field`); null: no arrowheads. */
let dirs: Int8Array | null = null;
/** An arrowhead (a filled triangle, device px) from its base x to its tip, over [y0, y1). */
function arrowhead(base: number, tip: number, y0: number, y1: number, ink: number, alpha: number) {
  heads.begin(ink, alpha, PathMode.fill);
  let o = heads.verb(PathVerb.move);
  heads.coords[o] = base;
  heads.coords[o + 1] = y0;
  o = heads.verb(PathVerb.line);
  heads.coords[o] = tip;
  heads.coords[o + 1] = (y0 + y1) >> 1;
  o = heads.verb(PathVerb.line);
  heads.coords[o] = base;
  heads.coords[o + 1] = y1;
  heads.close();
}

// A boundary over its track's letters (`Axis.lettered`, at the letters LOD): its form and geometry, per block.
const GLYPH = 0;
const PIN = 1;
const LINE = 2;
// Read once: a per-frame read of a module object's member allocated under the test runner's module wrapper.
const LINEAR = AxisGeometry.linear;
const LETTERS = Lod.letters;
const FILL = PathMode.fill;
const MOVE = PathVerb.move;
const LINE_TO = PathVerb.line;
/** The hit's reach each side of a pinned boundary (CSS px): the triangle's half width at most, and a few px. */
const PIN_REACH = 4.5;
/**
 * The block's boundary form (`pin[0]`: `PIN`, `LINE`, else `GLYPH`) and, pinned, in device px: [1] the triangle's height
 * and [2] its half width, [3] the letters' ink top and [4] its bottom (a descender's room), [5] the room each side of
 * a boundary between two letters' advances, [6] the row's top and [7] its bottom, [8] the hairline's width.
 */
const pin = new Float64Array(9);
/**
 * The letters' font `pinFor` last read (id, DPR), its capital's ink height (device px) and its advance (CSS px); the
 * DPR's hairline width, triangle height and half width at most (device px); the text manager's version it was read
 * at, and the manager (`inkText`): a font id names a font of one manager only, and a version drops its metrics (VL-09).
 * Held weakly: a disposed viewer's manager is not kept alive by the last pin drawn (LC-09).
 */
const inkOf = new Float64Array([-1, 0, 0, 0, 0, 0, 0, -1]);
let inkText: WeakRef<TextManagerApi> | null = null;
/** Whether layer `layer` draws over its track's letters at their LOD (draw and hit alike). */
function lettered(v: AxisView, layer: number): boolean {
  const ax = v.axis;
  return (
    ax !== null &&
    v.geometry === LINEAR &&
    v.lod === LETTERS &&
    v.height >= 10 &&
    layer >= 0 &&
    layer < ax.lettered.length &&
    ax.lettered[layer] === 1
  );
}
/**
 * `pin` for block b: over letters, the room above their text (the ink of the one-letter font their row takes at most,
 * as letters.ts's `textPx(1)` and `letterLines` place it) holds a triangle down to a device px above the ink; less than
 * 3 CSS px of it (a triangle too small to read), the hairline alone. The form.
 */
function pinFor(v: DrawView, b: number): number {
  pin[0] = GLYPH;
  if (!lettered(v, v.layer)) return GLYPH;
  // Integer and branch arithmetic: a Math builtin's double boxed under the test runner's module wrapper.
  const dpr = v.dpr;
  const h0 = v.height;
  const pad = h0 > 10 ? 2 : 0;
  let f = v.top * dpr + 0.5;
  let top = f | 0;
  if (top > f) top--;
  f = (v.top + pad) * dpr + 0.5;
  let y0 = f | 0;
  if (y0 > f) y0--;
  f = (v.top + h0 - pad) * dpr + 0.5;
  let y1 = f | 0;
  if (y1 > f) y1--;
  f = (v.top + h0) * dpr + 0.5;
  let bottom = f | 0;
  if (bottom > f) bottom--;
  const cell = v.cell * v.span[b];
  const p = v.props;
  const cap = p.maxLetterPx > 0 ? p.maxLetterPx : p.labelFontPx;
  // letters.ts's `textPx(1)`: the cell's 1.25, the row's 0.72, the cap, at least 6. The rule is `letterFontPx`'s
  // (letter-text.ts), inline here: through a call its doubles would be boxed under the test runner (LC-13).
  let px = cell * 1.25 < h0 * 0.72 ? cell * 1.25 : h0 * 0.72;
  px = px | 0;
  if (px > cap) px = cap;
  if (px < 6) px = 6;
  const font = v.text.font("mono", px | 0);
  // The font's ink and advance and the DPR's sizes, read when the font or DPR changes (a call's double is boxed).
  const tm = v.text;
  if (font !== inkOf[0] || dpr !== inkOf[1] || tm !== inkText?.deref() || tm.version.value !== inkOf[7]) {
    inkText = new WeakRef(tm);
    inkOf[7] = tm.version.value;
    inkOf[0] = font;
    inkOf[1] = dpr;
    inkOf[2] = v.text.ink(font, Math.round(dpr * 64))[0];
    inkOf[3] = v.text.cell(font, null);
    inkOf[4] = Math.max(1, Math.round(dpr));
    inkOf[5] = Math.round(5 * dpr);
    inkOf[6] = Math.round(4 * dpr);
  }
  const ink = inkOf[2];
  // `letterLines`' capital baseline: half way between two px, the higher.
  f = (y0 + y1 + ink) / 2 - 0.5;
  let base = f | 0;
  if (base < f) base++;
  const gap = inkOf[4];
  f = ink / 4;
  let below = f | 0;
  if (below < f) below++;
  pin[3] = base - ink;
  pin[4] = base + below;
  pin[5] = ((cell - inkOf[3]) * dpr) / 2;
  pin[6] = top;
  pin[7] = bottom;
  pin[8] = gap;
  let h = pin[3] - top - gap;
  if (h > inkOf[5]) h = inkOf[5];
  h = h | 0;
  let half = (h * 0.9 + 0.5) | 0;
  if (half > inkOf[6]) half = inkOf[6];
  pin[1] = h;
  pin[2] = half;
  pin[0] = h >= 3 * dpr ? PIN : LINE;
  return pin[0];
}
/**
 * A pinned boundary on block b's shown edge e, as `pin` says: the triangle's base [`tri[0]`, `tri[1]`) on the row's top
 * edge and its tip at `tri[2]`, the boundary's x. At the block's ends the triangle's half outside the data area is left
 * out (a half triangle against the edge), so it is never cut.
 */
const tri = new Float64Array(3);
function pinPlace(v: DrawView, b: number, e: number): void {
  let x = v.xd(b, e as never);
  if (x < edge[0]) x = edge[0];
  if (x > edge[1]) x = edge[1];
  const half = pin[0] === PIN ? pin[2] : 0;
  const end = e === 0 || e === (v.axis as Axis).shown[b];
  tri[0] = end && x - half < edge[0] ? x : x - half;
  tri[1] = end && x + half > edge[1] ? x : x + half;
  tri[2] = x;
}
/** A pinned boundary: its hairline down the row, and its triangle where `pin` has room for one. */
function pinned(v: DrawView, b: number, e: number, ink: number, alpha: number) {
  const w = pin[8];
  pinPlace(v, b, e);
  const x = tri[2];
  // The hairline, a device column on the boundary (inside the data area); where the letters leave it no room, above
  // and below their ink.
  let xl = ((x + 0.5) | 0) - (w >> 1);
  if (xl < edge[0]) xl = edge[0];
  if (xl > edge[1] - w) xl = edge[1] - w;
  const y0 = pin[6];
  const y1 = pin[7];
  rects.begin(ink, alpha);
  if (pin[5] >= w) rects.push(xl, y0, w, y1 - y0);
  else {
    const a = pin[3] - w;
    const c = pin[4] + w;
    if (a > y0) rects.push(xl, y0, w, a - y0);
    if (y1 > c) rects.push(xl, c, w, y1 - c);
  }
  if (pin[0] !== PIN) return;
  heads.begin(ink, alpha, FILL);
  let o = heads.verb(MOVE);
  heads.coords[o] = tri[0];
  heads.coords[o + 1] = y0;
  o = heads.verb(LINE_TO);
  heads.coords[o] = tri[1];
  heads.coords[o + 1] = y0;
  o = heads.verb(LINE_TO);
  heads.coords[o] = x;
  heads.coords[o + 1] = y0 + pin[1];
  heads.close();
}
/**
 * A mark's center x on block b's shown edge e, in `gx`: [0] the x, [1] its half width, [2] and [3] the data area's
 * edges; at the block's ends, [0] moved inside them. Through scratch: a double passed to a call is boxed.
 */
const gx = new Float64Array(4);
function inside(v: AxisView, b: number, e: number): void {
  if (e !== 0 && e !== (v.axis as Axis).shown[b]) return;
  if (gx[0] < gx[2] + gx[1]) gx[0] = gx[2] + gx[1];
  else if (gx[0] > gx[3] - gx[1]) gx[0] = gx[3] - gx[1];
}

/** An item: one rect per piece in view, a bracket's three, or a boundary's glyph; its label where it fits. */
function item(v: DrawView, f: Features, k: number, b: number, y0: number, y1: number, brackets: boolean, font: number) {
  const dpr = v.dpr;
  const key = featureKey(v, f, k);
  const ink = featureInk(key);
  const alpha = keyAlpha(key);
  if (f.stop[k] <= f.start[k]) {
    // A boundary: the glyph (a tick by default) on the edge before its position (a reversed piece's: its right).
    const e = boundaryEdge(v, b, f.start[k]);
    if (e < 0 || e < range[0] || e > range[1]) return;
    // Over letters, pinned clear of their text.
    if (pin[0] !== GLYPH) {
      pinned(v, b, e, ink, alpha);
      return;
    }
    const id = featureGlyph(v, f, k, 7);
    const r = Math.min(sizePx(featureSize(v, f, k, sizeTheme(v, 1, 2)), 5, 1, 12) * dpr, (y1 - y0) / 2);
    glyphs.begin(ink, alpha);
    const i = glyphs.add(id < 0 ? 7 : id);
    gx[0] = v.xd(b, e as never);
    gx[1] = r;
    gx[2] = edge[0];
    gx[3] = edge[1];
    inside(v, b, e);
    glyphs.x[i] = gx[0];
    glyphs.y[i] = (y0 + y1) / 2;
    glyphs.radius[i] = r;
    return;
  }
  // The thickness channel: a bar of a share of its lane, centered (brackets keep their own).
  if (!brackets) {
    barBand(featureSize(v, f, k, sizeTheme(v, 0, 2)), y0, y1);
    y0 = bar[0];
    y1 = bar[1];
  }
  const n = placeRun(v, b, f.start[k], f.stop[k], piece);
  // A direction (`direction_field`): an arrowhead at the item's end, its tip on the end's outer edge on screen (a
  // reversed piece's left), cut from the piece that reaches it, at most 3/4 of the bar's height long.
  const d = dirs === null || brackets ? 0 : dirs[k];
  const qe = d === 0 ? -1 : v.shownOf(b, d > 0 ? f.stop[k] - 1 : f.start[k]);
  const sd = qe < 0 ? 0 : d * v.dir[b];
  const tip = sd > 0 ? v.xd(b, qe + 1) : sd < 0 ? v.xd(b, qe) : 0;
  const hl = Math.max(2, ((y1 - y0) * 3) >> 2);
  rects.begin(ink, alpha);
  for (let q = 0; q < n; q++) {
    let x0 = piece[2 * q] | 0;
    const span = (piece[2 * q + 1] | 0) - x0;
    // A device px apart from the next, unless the piece is cut by the data area's edge (it ends there).
    let w = Math.max(1, span - (span > 3 * dpr && x0 + span < edge[1] ? 1 : 0));
    // Only on a piece half again as long as the head (a narrower one stays a bar: no spikes when zoomed out).
    if (sd > 0 && 2 * w >= 3 * hl && Math.abs(x0 + span - tip) <= 1) {
      w -= hl;
      arrowhead(x0 + w, x0 + w + hl, y0, y1, ink, alpha);
    } else if (sd < 0 && 2 * w >= 3 * hl && Math.abs(x0 - tip) <= 1) {
      arrowhead(x0 + hl, x0, y0, y1, ink, alpha);
      x0 += hl;
      w -= hl;
    }
    if (brackets) {
      const t = Math.min(Math.round(4 * dpr) | 0, y1 - y0);
      rects.push(x0, y1 - 1, w, 1);
      rects.push(x0, y1 - t, 1, t);
      rects.push(x0 + w - 1, y1 - t, 1, t);
    } else rects.push(x0, y0, w, y1 - y0);
  }
  const label = f.label[k];
  if (font < 0 || label === null || !extent(v, b, f.start[k], f.stop[k])) return;
  // Centered on the piece's whole extent (its bar: a device px off a wide one's end), so it moves with the content; it
  // sticks 4 px inside the data area's edge only where it would leave the piece's visible part. Cut with an ellipsis to
  // the visible part, down to three graphemes (Nightingale's).
  const a0 = ext[0];
  const a1 = ext[1] - (ext[1] - a0 > 3 * dpr ? 1 : 0);
  const v0 = a0 > edge[0] ? a0 : edge[0];
  const v1 = a1 < edge[1] ? a1 : edge[1];
  const p = v.text.label(label, font);
  fit[0] = (v1 - v0) / dpr - 8;
  const m = fitIn(p, fit);
  if (m === 0 || m < Math.min(3, p.ends.length)) return;
  // Half its width, then its center (device px), through scratch: no double phi for Maglev to box (RN1).
  at[0] = p.ends[m - 1];
  if (m < p.ends.length) at[0] += fit[1];
  at[0] *= dpr / 2;
  at[1] = (a0 + a1) / 2;
  const pad = 4 * dpr;
  names.begin(brackets ? SkinToken.ink : SkinToken.surface);
  const t = names.add(fitted(p, m));
  names.x[t] = at[1];
  if (at[1] - at[0] < v0 + pad) names.x[t] = v0 + pad + at[0];
  else if (at[1] + at[0] > v1 - pad) names.x[t] = v1 - pad - at[0];
  names.y[t] = Math.round((y0 + y1) / 2) + 0.5;
}

/**
 * Into `ext`: the whole extent (device px, linear) of the first shown run of [start, stop) that reaches into the data
 * area, never cut by the view (nor by a copy's strip: what a full draw places); false: none does.
 */
function extent(v: DrawView, b: number, start: number, stop: number): boolean {
  let n = v.runToShown(b, start, stop, shown);
  if (2 * n > shown.length) {
    shown = new Int32Array(4 * n);
    n = v.runToShown(b, start, stop, shown);
  }
  for (let r = 0; r < n; r++) {
    const x0 = v.xd(b, shown[2 * r]);
    const x1 = v.xd(b, shown[2 * r + 1]);
    if (x1 <= edge[0] || x0 >= edge[1]) continue;
    ext[0] = x0;
    ext[1] = x1;
    return true;
  }
  return false;
}

/** A hairline across the gap between an item's fragments, at the lane's middle. */
function join(v: DrawView, f: Features, k: number, b: number, y: number) {
  const j = f.joinFrom[k];
  if (j < 0 || (f.on[k] >= 0 && f.on[k] !== v.row) || !owned(k)) return;
  const n = placeRun(v, b, f.stop[j], f.start[k], piece);
  const key = featureKey(v, f, k);
  rects.begin(featureInk(key), keyAlpha(key));
  for (let q = 0; q < n; q++)
    rects.push(piece[2 * q] | 0, y, Math.max(1, (piece[2 * q + 1] | 0) - (piece[2 * q] | 0)), 1);
}

/**
 * Every item in view, lane by lane (`flat`: the tiny path, one band, no brackets or labels; `lanes`: a streamed layer's,
 * its region's table drawn, only the items it owns).
 */
function items(v: DrawView, f: Features, p: Params, flat: boolean, lanes = f.lanes) {
  const ax = v.axis;
  if (!ax) return;
  const ys = laneRows(v, f, flat, lanes);
  edge[0] = Math.round(v.left * v.dpr) | 0;
  edge[1] = Math.round(v.right * v.dpr) | 0;
  const brackets = !flat && p.display === "brackets";
  dirs = flat ? null : directionsOf(v, f, p);
  // Labels need a lane of 10 px and over; one font for the draw. On rings they are the section's, outside (ring-labels.ts).
  const font = !flat && !v.ring && p.labels === true && ys[1] - ys[0] >= 10 * v.dpr ? labelFont(v, ys[1] - ys[0]) : -1;
  for (let b = 0; b < ax.count && (b + 1) * f.lanes < f.seg.length; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    shownHull(v, b, range[0], range[1], hull);
    if (flat) pin[0] = GLYPH;
    else pinFor(v, b);
    for (let lane = 0; lane < f.lanes && lane < lanes; lane++) {
      const g = b * f.lanes + lane;
      const y0 = ys[2 * lane];
      const y1 = ys[2 * lane + 1];
      // One past the hull: a boundary after the last position in view (at the piece's end) is drawn on its edge.
      segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1] + 1, win);
      // The item after the window may join from an item inside it.
      for (let k = Math.max(win[0], f.seg[g] + 1); k <= win[1] && k < f.seg[g + 1]; k++)
        join(v, f, k, b, (y0 + y1) >> 1);
      for (let k = win[0]; k < win[1]; k++) {
        if (
          (f.stop[k] <= f.start[k] ? f.start[k] < hull[0] : f.stop[k] <= hull[0] || f.start[k] >= hull[1]) ||
          (f.on[k] >= 0 && f.on[k] !== v.row) ||
          !owned(k)
        )
          continue;
        item(v, f, k, b, y0, y1, brackets, font);
      }
    }
  }
}
const region = (v: DrawView, f: Features, p: unknown, lanes: number) => items(v, f, p as Params, false, lanes);
const regionFlat = (v: DrawView, f: Features, p: unknown, lanes: number) => items(v, f, p as Params, true, lanes);

/** The labels' font for a lane of `h` device px, set on the labels' batch; its id. */
function labelFont(v: DrawView, h: number): number {
  const font = v.text.font("ui", Math.max(7, Math.min(11, Math.floor((h * 0.72) / v.dpr))));
  names.style(font, TextAlign.center, TextBaseline.middle);
  fit[1] = v.text.label("…", font).width;
  return font;
}

function draw(v: DrawView, f: Features, p: Params, flat: boolean) {
  rects.reset();
  glyphs.reset();
  names.reset();
  heads.reset();
  ThemeLocation.reset(v.loc);
  if (v.geometry === AxisGeometry.wrapped) drawLineLanes(v, f, rects);
  // Streamed (ST3a): region by region (or density), pending regions as placeholders.
  else if (drawStreamedItems(v, p, false, false, flat ? regionFlat : region, rects)) {
  } else {
    // Items while no more are in view than device columns, else density (one mode for every block, with hysteresis).
    const columns = densityColumns(v);
    if (v.reductions && densityMode(v, itemsInView(v, f, columns + 1), columns)) densityBand(v, f, p, false, rects);
    else items(v, f, p, flat);
  }
  if (rects.n > 0) v.backend.rects(rects);
  if (heads.n > 0) v.backend.paths(heads);
  if (glyphs.n > 0) {
    glyphs.font = v.text.font("ui", 12, 600);
    v.backend.glyphs(glyphs);
  }
  if (names.n > 0) v.backend.text(names);
}

/**
 * The lane under y, then by drawn geometry: a boundary's glyph holding the point (the nearest), else the feature
 * covering the position there, the one drawn on top (`featureAt`), when the point is on its band (a thickness theme's
 * centered share of the lane). A gap between lanes or features, or off a band, is the position's hover; past the
 * block's end (no position), only a glyph counts.
 */
function hit(v: AxisView, f: Features, p: Params, out: Hit): boolean {
  // Wrapped, its underlines on the line's lanes.
  if (hitLineLanes(v, f, out)) return true;
  // On a rows track's row, its features are found by lane (hit.ts's item lookup); streamed, in its regions (ST3b).
  if (v.geometry !== AxisGeometry.linear || out.block < 0 || v.row !== -1) return false;
  // Past the block's end (no position, `hitsPastEnd`): a boundary's glyph only.
  const past = out.index < 0;
  if (!past && hitStreamed(v, out, true)) return true;
  // Over letters, a boundary is pinned: within a few px of it over the row, in any lane.
  const over = lettered(v, hitLayer.layer);
  if (over) {
    const p = pinAt(v, f, out);
    if (p >= 0) {
      hitItem(out, (v.axis as Axis).set, f.annotation[p], f.row[p]);
      return true;
    }
  }
  const h = laneHeight(v.height, f.lanes);
  const at = (out.point[1] - v.top - lanePad(v.height)) / h;
  const lane = Math.floor(at);
  const inset = laneInset(h) / h;
  let k =
    at - lane >= inset && at - lane < 1 - inset && lane < f.lanes ? (over ? -1 : boundaryAt(v, f, lane, h, out)) : -2;
  if (k === -1 && !past) k = featureAt(f, out.block, lane, out.index);
  // Over letters a boundary is its pin alone: the position after it is its residue's.
  if (over && k >= 0 && f.stop[k] <= f.start[k]) k = -1;
  if (past && k < 0) return false;
  const flat = v.height < (RangesRenderer.minHeight ?? 0);
  if (k >= 0 && f.stop[k] > f.start[k] && !flat && p.display !== "brackets") {
    // The bar's band as the draw places it (`laneRows`, then `barBand`; device px).
    const top = v.top + lanePad(v.height) + lane * h;
    const y0 = Math.round((top + inset * h) * v.dpr) | 0;
    const y1 = Math.max(y0 + 1, Math.round((top + h - inset * h) * v.dpr) | 0);
    barBand(featureSize(hitLayer, f, k, sizeTheme(hitLayer, 0, 2)), y0, y1);
    const y = out.point[1] * v.dpr;
    if (y < bar[0] || y >= bar[1]) k = -1;
  }
  if (k >= 0) hitItem(out, (v.axis as Axis).set, f.annotation[k], f.row[k]);
  return true;
}
/** The boundary in `lane` whose glyph (as `item` draws it) holds the hit's point, the nearest; -1: none. */
function boundaryAt(v: AxisView, f: Features, lane: number, h: number, out: Hit): number {
  const b = out.block;
  if (!hitWindow(v, b, out.point[0], 12, near)) return -1;
  const g = b * f.lanes + lane;
  const cy = v.top + lanePad(v.height) + (lane + 0.5) * h;
  const t = sizeTheme(hitLayer, 1, 2);
  let best = -1;
  let least = Number.POSITIVE_INFINITY;
  segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], near[0], near[1] + 1, win);
  for (let k = win[0]; k < win[1]; k++) {
    const e = f.stop[k] <= f.start[k] ? boundaryEdge(v, b, f.start[k]) : -1;
    if (e < 0) continue;
    const r = Math.min(sizePx(featureSize(hitLayer, f, k, t), 5, 1, 12), h / 2 - laneInset(h));
    gx[0] = v.x(b, e as never);
    gx[1] = r;
    gx[2] = v.left;
    gx[3] = v.right;
    inside(v, b, e);
    const dx = out.point[0] - gx[0];
    const dy = out.point[1] - cy;
    const d = dx * dx + dy * dy;
    if (d <= (r + 0.5) * (r + 0.5) && d < least) {
      least = d;
      best = k;
    }
  }
  return best;
}
const near = new Int32Array(2);
/** The pinned boundary (any lane) within `PIN_REACH` of the hit's x on the row, the nearest; -1: none. */
function pinAt(v: AxisView, f: Features, out: Hit): number {
  const b = out.block;
  if (out.point[1] < v.top || out.point[1] >= v.top + v.height || !hitWindow(v, b, out.point[0], PIN_REACH, near))
    return -1;
  let best = -1;
  let least = PIN_REACH + 0.5;
  for (let lane = 0; lane < f.lanes; lane++) {
    const g = b * f.lanes + lane;
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], near[0], near[1] + 1, win);
    for (let k = win[0]; k < win[1]; k++) {
      if (f.stop[k] > f.start[k] || (f.on[k] >= 0 && f.on[k] !== v.row)) continue;
      const e = boundaryEdge(v, b, f.start[k]);
      if (e < 0) continue;
      gx[0] = v.x(b, e as never);
      gx[1] = 0;
      gx[2] = v.left;
      gx[3] = v.right;
      inside(v, b, e);
      const d = out.point[0] > gx[0] ? out.point[0] - gx[0] : gx[0] - out.point[0];
      if (d <= least) {
        least = d;
        best = k;
      }
    }
  }
  return best;
}

const hovered = new Int32Array(2);
/**
 * The hovered item outlined inside its own pixels: each piece's bar as the draw places it, a bracket's box, a
 * boundary's glyph outlined inset; flat below the minimum height, as `tiny` draws; wrapped, its underlines. Nothing in
 * density.
 */
function overlay(v: DrawView, f: Features, p: Params) {
  if (!hoverEmphasis(v) || !v.axis || outlineStreamed(v, v.height >= (RangesRenderer.minHeight ?? 0))) return;
  if (v.geometry === AxisGeometry.wrapped) {
    outlineLineLanes(v, f);
    return;
  }
  const x = v.geometry === AxisGeometry.linear ? hoveredItem(v, f, hovered) : null;
  const flat = v.height < (RangesRenderer.minHeight ?? 0);
  if (x === null || (!flat && v.reductions && densityShown(v))) return;
  const dpr = v.dpr;
  const ys = laneRows(v, f, flat);
  const right = Math.round(v.right * dpr) | 0;
  edge[0] = Math.round(v.left * dpr) | 0;
  edge[1] = right;
  const brackets = !flat && p.display === "brackets";
  HoverOutline.reset();
  for (let j = hovered[0]; j < hovered[1]; j++) {
    const k = hoveredAt(v, f, x, j);
    if (k < 0) continue;
    const b = f.block[k];
    const y0 = ys[2 * f.lane[k]];
    const y1 = ys[2 * f.lane[k] + 1];
    const key = featureHoverKey(v, f, k);
    if (f.stop[k] <= f.start[k]) {
      // A boundary's glyph, as the draw places it, in a cell of the lane's height; pinned, its triangle or hairline.
      const e = boundaryEdge(v, b, f.start[k]);
      if (e < 0) continue;
      const x = v.xd(b, e as never);
      if (pinFor(v, b) !== GLYPH) {
        pinPlace(v, b, e);
        const t0 = (tri[0] + 0.5) | 0;
        if (pin[0] === PIN) HoverOutline.box(v, key, t0, pin[6], ((tri[1] + 0.5) | 0) - t0, pin[1]);
        else HoverOutline.box(v, key, ((tri[2] + 0.5) | 0) - (pin[8] >> 1), pin[6], pin[8], pin[7] - pin[6]);
        continue;
      }
      const id = featureGlyph(v, f, k, 7);
      const r = Math.min(sizePx(featureSize(v, f, k, sizeTheme(v, 1, 2)), 5, 1, 12) * dpr, (y1 - y0) / 2) | 0;
      gx[0] = x;
      gx[1] = r;
      gx[2] = edge[0];
      gx[3] = edge[1];
      inside(v, b, e);
      HoverOutline.glyph(v, key, id < 0 ? 7 : id, gx[0], (y0 + y1) >> 1, r, r);
      continue;
    }
    const n = placeRun(v, b, f.start[k], f.stop[k], piece);
    for (let q = 0; q < n; q++) {
      const x0 = piece[2 * q] | 0;
      const span = (piece[2 * q + 1] | 0) - x0;
      const w = Math.max(1, span - (span > 3 * dpr && x0 + span < right ? 1 : 0));
      barBand(brackets ? NO_SIZE : featureSize(v, f, k, sizeTheme(v, 0, 2)), y0, y1);
      const t = brackets ? Math.min(Math.round(4 * dpr) | 0, y1 - y0) : bar[1] - bar[0];
      HoverOutline.box(v, key, x0, (brackets ? y1 : bar[1]) - t, w, t);
    }
  }
  HoverOutline.hand(v);
}

export const RangesRenderer = kernelEntry<Features>("ranges", 22, {
  geometries: ["linear", "wrapped"],
  hit: (v, f, p, out) => hit(v, f, p as Params, out),
  hitsPastEnd: true,
  draw: (v, f, p) => draw(v, f, p as Params, false),
  tiny: (v, f, p) => draw(v, f, p as Params, true),
  overlay: (v, f, p) => overlay(v, f, p as Params),
  sized: ["thickness", "radius"],
  // With `labels`, on rings: each labeled item, for the section's labels outside the rings.
  ringLabels: (f, p, out) => {
    if ((p as Params).labels !== true) return;
    for (let k = 0; k < f.count; k++) {
      const label = f.label[k];
      if (label) out.add(f.block[k], f.start[k], f.stop[k], label, f.annotation[k], f.row[k]);
    }
  },
});
