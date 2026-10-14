/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `pairs` renderer (seqspec:viewer#extension-points): arcs within and between blocks, ends with glyphs; the pairs
// in view, the ends in view and an arc's hit all come from the end index and interval tree prepare built
// (`pairsInView`, `endsIn`, `pairsAt`), shared with the overlay (`litPairs`, `pushArc`). Arcs are one `PathBatch`
// sorted by color (a counting sort of the pairs in view), so a backend strokes each color once; below a device px a
// cell, and below the row's minimum (`tiny`), the ends are density per device column. The overlay strokes the arcs of
// pairs with a position of either end lit (`litPairs`, through the end index: never every pair), the same curves,
// stronger.
import { Annotations } from "@molstar/seq-core/model";
import {
  Arcs,
  AxisGeometry,
  type AxisView,
  CellIO,
  cellLayout,
  DensityKind,
  DOT_R,
  type DrawView,
  densityAdd,
  densityDraw,
  densitySlot,
  densityStart,
  endsIn,
  GlyphBatch,
  type Hit,
  HOVER_WIDTH,
  hitItem,
  hitWindow,
  inkOf,
  kernelEntry,
  keyAlpha,
  keyOf,
  litPairs,
  MarkBit,
  NO_COLOR,
  NO_SIZE,
  type PairsState,
  PathBatch,
  packAlpha,
  pairsAt,
  pairsInView,
  pushArc,
  RectBatch,
  SIZE_ONE,
  SkinToken,
  sizeTheme,
  TextAlign,
  TextBaseline,
  TextBatch,
  ThemeLocation,
  visibleRange,
} from "../draw/api.ts";

type Params = Readonly<Record<string, unknown>>;

/** A pair end's glyph radius: the viewer's dot. */
const R = DOT_R;
const io = new Float64Array(CellIO.SIZE);
const elementAt = new Int32Array(4);
/** The "+n" of a split cell's hidden kinds, n to 99: strings held, never built per draw. */
const MORE = Array.from({ length: 100 }, (_, n) => `+${n}`);
const MARKER_PX = 9;
/** On rings, the most pairs whose ends are drawn as glyphs below a device px a cell. */
const RING_GLYPHS = 256;

/**
 * How the cell of end `e` splits into one glyph per kind (`cellLayout`): the glyphs shown, 0 when it is not split; `io`
 * then holds the pitch, the marker and the cell's center x. Its room is the cell's width on screen.
 */
function splitAt(v: AxisView, s: PairsState, e: number): number {
  const c = s.endCell[s.endAt[e]];
  if (s.cellKinds[c] < 2) return 0;
  const b = s.block[e];
  const k0 = v.shownOf(b, s.start[e]);
  const k1 = v.shownOf(b, s.stop[e] - 1);
  if (k0 < 0 || k1 < 0) return 0;
  // A reversed piece shows an end's last residue first: the room and the center are the same either way.
  const lo = Math.min(k0, k1);
  const hi = Math.max(k0, k1);
  // The room is the cells' (not the snapped edges' difference, which a mirrored draw rounds the other way).
  io[CellIO.width] = v.cell * (hi - lo + 1);
  io[CellIO.center] = (v.x(b, (hi + 1) as never) + v.x(b, lo as never)) / 2;
  io[CellIO.widest] = 2 * R;
  return cellLayout(s, c, io);
}
/**
 * The center x of slot `i` of a split cell (`splitAt` has run for end `e`) into `io[at]`: from the
 * cell's left in the order of the kinds, from its right in a reversed piece (so a reversed draw mirrors the forward
 * one); slot -1 is the marker's.
 */
function slotX(v: AxisView, s: PairsState, e: number, shown: number, i: number) {
  const pitch = io[CellIO.pitch];
  const marker = io[CellIO.marker];
  const left = io[CellIO.center] - (pitch * shown + marker) / 2;
  const rev = v.dir[s.block[e]] < 0;
  io[CellIO.at] =
    i < 0
      ? rev
        ? left + marker / 2
        : left + pitch * shown + marker / 2
      : left + (rev ? marker + pitch * (shown - 1 - i + 0.5) : pitch * (i + 0.5));
}

let found = new Int32Array(64);
/**
 * By drawn geometry: the end whose glyph (its radius theme's size) holds the point, the nearest, else the end at the
 * hit's position; else the topmost drawn arc within its half width (at least 4 px) of the point: the overlay's last lit
 * arc, else the draw's (the last of the largest color rank, by first key), as `draw` and `overlay` stack them. The pair
 * as the item.
 */
function hit(v: AxisView, s: PairsState, p: Params, out: Hit): boolean {
  if (v.geometry !== AxisGeometry.linear || out.block < 0 || s.count === 0) return false;
  Arcs.frame(v, s, p);
  const f = s.frame;
  const c = colors.get(s);
  if (p.display !== "arcs" && Math.abs(out.point[1] - f[0]) <= (c?.radius ? 8 : 5)) {
    let at = endNear(v, s, c, out, f[0]);
    const near = at >= 0;
    const q = keyOf(out.block, out.index);
    for (let lo = 0, hi = s.endKey.length - 1; !near && lo <= hi; ) {
      const mid = (lo + hi) >> 1;
      if (s.endKey[mid] <= q) {
        at = mid;
        lo = mid + 1;
      } else hi = mid - 1;
    }
    const e = at >= 0 ? s.endOf[at] : -1;
    if (e >= 0 && s.block[e] === out.block && (near || s.stop[e] > out.index)) {
      Arcs.end(v, s, e >> 1, e & 1);
      const shown = splitAt(v, s, e);
      if (shown > 0) {
        // A cell with several kinds: the glyph under the pointer, its kind in the hit's `part` (the overlay lights
        // that kind's pairs at this end). The marker is no target; past it, the arcs below.
        let slot = -1;
        for (let i = 0; i < shown; i++) {
          slotX(v, s, e, shown, i);
          if (Math.abs(out.point[0] - io[CellIO.at]) <= io[CellIO.pitch] / 2) slot = i;
        }
        if (slot < 0) return p.display === "ends" || chords.has(s) ? false : hitArcs(v, s, p, out, f);
        const c = s.endCell[at];
        const rep = s.endOf[s.kindRep[s.cellKindAt[c] + slot]];
        out.side = rep & 1;
        out.part = s.kindList[s.cellKindAt[c] + slot];
        return hitItem(out, s.set, s.annotation[rep >> 1], s.row[rep >> 1]);
      }
      out.side = e & 1;
      return hitItem(out, s.set, s.annotation[e >> 1], s.row[e >> 1]);
    }
  }
  // On rings the arcs are chords through the hub, which hits do not reach: the ends only.
  if (p.display === "ends" || chords.has(s)) return false;
  return hitArcs(v, s, p, out, f);
}

/** The end (its place in key order) whose glyph holds the point, the nearest; -1: none. `y`: the ends' line. */
function endNear(v: AxisView, s: PairsState, c: PairColors | undefined, out: Hit, y: number): number {
  const b = out.block;
  if (!hitWindow(v, b, out.point[0], c?.radius ? 8 : R, nearby)) return -1;
  // The ends in that window, by their keys.
  let lo = 0;
  for (let hi = s.endKey.length; lo < hi; ) {
    const mid = (lo + hi) >> 1;
    if (s.endKey[mid] < keyOf(b, nearby[0])) lo = mid + 1;
    else hi = mid;
  }
  let best = -1;
  let least = Number.POSITIVE_INFINITY;
  for (let j = lo; j < s.endKey.length && s.endKey[j] < keyOf(b, nearby[1]); j++) {
    const e = s.endOf[j];
    Arcs.end(v, s, e >> 1, e & 1);
    const dx = out.point[0] - s.xs[e & 1];
    const dy = out.point[1] - y;
    const size = c?.radius ? c.radius[e >> 1] : NO_SIZE;
    const r = size < 0 ? R : Math.min(8, Math.max(1, size / SIZE_ONE));
    const d = dx * dx + dy * dy;
    if (d <= (r + 0.5) * (r + 0.5) && d < least) {
      least = d;
      best = j;
    }
  }
  return best;
}
const nearby = new Int32Array(2);

/** The topmost drawn arc within its half width (at least 4 px) of the point (`f`: the frame, after `Arcs.frame`). */
function hitArcs(v: AxisView, s: PairsState, p: Params, out: Hit, f: Float64Array): boolean {
  // The arcs over the position, from the interval tree (a stab): never a scan of the earlier pairs.
  if (found.length < s.count) found = new Int32Array(s.count);
  const n = pairsAt(s, out.block, out.index, found);
  let best = -1;
  let top = -1;
  const hover = p.arcs === "hover";
  // The draw's last color ranks (only drawn arcs hit, so the draw has made them); the overlay's arcs over all of them.
  const c = colors.get(s);
  const over = c ? c.ranks : 0;
  for (let j = 0; j < n; j++) {
    const k = found[j];
    const z = s.drawn[k] > 0 ? over + s.drawn[k] : hover ? -1 : c ? c.rank[k] : 0;
    // Found by first key, as the draw stacks a color's arcs: a later one of the same place is on top.
    if (z < top || z < 0 || !Arcs.ends(v, s, k)) continue;
    const xa = s.xs[0];
    const xb = s.xs[1];
    if (xa === xb) continue;
    // The draw's curve (`Arcs.heightAt`).
    Arcs.apex(s, c?.height ? c.height[k] : NO_SIZE);
    f[10] = out.point[0];
    Arcs.heightAt(s);
    const h = f[11];
    if (h < 0) continue;
    const w = c?.width ? c.width[k] / SIZE_ONE / 2 : 0;
    if (Math.abs(out.point[1] - f[0] - f[1] * h) <= (w > 4 ? w : 4)) {
      // Under a width theme the draw stacks a color's arcs widest first: of two at one place, the thinner is on top
      // (VL-24), by this layer's own widths (the draw's `half` is whichever pairs layer drew last: LC-01).
      if (z === top && best >= 0 && c?.width && halfWidth(c, k, v.dpr, 0) > halfWidth(c, best, v.dpr, 0)) continue;
      best = k;
      top = z;
    }
  }
  return best >= 0 && hitItem(out, s.set, s.annotation[best], s.row[best]);
}

/** Per (prepared state, color theme): each pair's color key (alpha folded in, as density's) and its color's rank. */
interface PairColors {
  theme: unknown;
  /** The width, height and radius themes it was made under, and per pair their sizes (`NO_SIZE`: the viewer's); null: none. */
  sizes: readonly [unknown, unknown, unknown];
  width: Int32Array | null;
  height: Int32Array | null;
  radius: Int32Array | null;
  key: Int32Array;
  rank: Int32Array;
  ranks: number;
}
const colors = new WeakMap<PairsState, PairColors>();
/** Pair k's location: its first placed end's, its annotation and row. */
function locate(v: DrawView, s: PairsState, k: number) {
  const ax = v.axis;
  const e = s.block[2 * k] >= 0 ? 2 * k : 2 * k + 1;
  if (ax) ThemeLocation.at(v.loc, ax, v.layer, s.block[e], s.start[e]);
  v.loc.annotation = s.annotation[k];
  v.loc.row = s.row[k];
}
/** Edge-rate, apart from the draw. */
function newColors(v: DrawView, s: PairsState): PairColors {
  const t = v.color;
  const [wt, ht, rt] = [sizeTheme(v, 0, 3), sizeTheme(v, 1, 3), sizeTheme(v, 2, 3)];
  const key = new Int32Array(s.count);
  const rank = new Int32Array(s.count);
  const width = wt ? new Int32Array(s.count) : null;
  const height = ht ? new Int32Array(s.count) : null;
  const radius = rt ? new Int32Array(s.count) : null;
  const seen = new Map<number, number>();
  for (let k = 0; k < s.count; k++) {
    locate(v, s, k);
    if (width && wt) width[k] = wt.value(v.loc);
    if (height && ht) height[k] = ht.value(v.loc);
    if (radius && rt) radius[k] = rt.value(v.loc);
    const c = t ? t.value(v.loc) : NO_COLOR;
    const a = c >= 0 && t?.alpha ? t.alpha(v.loc) : 255;
    key[k] = packAlpha(c, a);
    if (!seen.has(key[k])) seen.set(key[k], seen.size);
    rank[k] = seen.get(key[k]) as number;
  }
  const x = { theme: t, sizes: [wt, ht, rt] as const, width, height, radius, key, rank, ranks: seen.size };
  colors.set(s, x);
  return x;
}
const colorsOf = (v: DrawView, s: PairsState): PairColors => {
  const x = colors.get(s);
  return x &&
    x.theme === v.color &&
    x.sizes[0] === sizeTheme(v, 0, 3) &&
    x.sizes[1] === sizeTheme(v, 1, 3) &&
    x.sizes[2] === sizeTheme(v, 2, 3)
    ? x
    : newColors(v, s);
};
const arcs = new PathBatch(64);
const wedges = new PathBatch(16);
const glyphs = new GlyphBatch(64);
const bars = new RectBatch(256);
const markers = new TextBatch(16);
let list = new Int32Array(64);
let sorted = new Int32Array(64);
let byRank = new Int32Array(64);
const range = new Int32Array(2);
const span = new Int32Array(2);

/**
 * The glyphs of the split cell holding place `q` (end `e`; `splitAt` and `Arcs.end` have run for it): one per kind, each
 * a pair of that kind (the first at the cell) drawn as its color and glyph theme say, side by side from the cell's
 * left, then the "+n" of the kinds that did not fit.
 */
function drawSplit(v: DrawView, s: PairsState, c: PairColors, q: number, e: number, shown: number) {
  const dpr = v.dpr;
  const gt = v.glyph;
  const cell = s.endCell[q];
  const pitch = io[CellIO.pitch];
  const radius = Math.max(1.5, pitch / 2 - 0.5);
  const hidden = io[CellIO.hidden];
  const y = s.frame[0];
  for (let i = 0; i < shown; i++) {
    const k = s.endOf[s.kindRep[s.cellKindAt[cell] + i]] >> 1;
    locate(v, s, k);
    const id = gt ? gt.value(v.loc) : -1;
    glyphs.begin(inkOf(c.key[k]), keyAlpha(c.key[k]));
    const g = glyphs.add(id < 0 ? 0 : id);
    slotX(v, s, e, shown, i);
    glyphs.x[g] = io[CellIO.at] * dpr;
    glyphs.y[g] = y * dpr;
    glyphs.radius[g] = radius * dpr;
  }
  if (hidden > 0) {
    markers.begin(SkinToken.inkQuiet);
    const t = markers.add(MORE[Math.min(hidden, 99)]);
    slotX(v, s, e, shown, -1);
    markers.x[t] = io[CellIO.at] * dpr;
    markers.y[t] = y * dpr;
  }
}

/**
 * Each end in view (the end index over every drawn block's shown range): its glyph at its residue's center on the
 * ends' line, or with `dense` counted in its device column (`densityStart` before).
 */
function drawEnds(v: DrawView, s: PairsState, c: PairColors, dense: boolean) {
  const ax = v.axis;
  const dpr = v.dpr;
  const gt = v.glyph;
  for (let b = 0; ax && b < ax.count; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    const i = ax.indexOf(b, range[0] as never);
    const j = ax.indexOf(b, (range[1] - 1) as never);
    const across = ax.wrapAt(b) > range[0] && ax.wrapAt(b) < range[1];
    const lo = across ? ax.base[b] : Math.min(i, j);
    endsIn(s, b, lo, across ? ax.base[b] + ax.length[b] : Math.max(i, j) + 1, span);
    for (let q = span[0]; q < span[1]; q++) {
      const e = s.endOf[q];
      if (s.block[e] !== b || s.stop[e] <= lo) continue;
      const k = e >> 1;
      Arcs.end(v, s, k, e & 1);
      const x = s.xs[e & 1];
      if (Number.isNaN(x) || x < v.left - R || x > v.right + R) continue;
      const key = c.key[k];
      if (dense) {
        densityAdd(Math.round(x * dpr) | 0, densitySlot(key), 1);
        continue;
      }
      const shown = s.cellKinds[s.endCell[q]] > 1 ? splitAt(v, s, e) : 0;
      if (shown > 0) {
        // Several kinds at this cell: one glyph each, side by side, drawn once, at the cell's first place in view.
        if (q === Math.max(s.cellLo[s.endCell[q]], span[0])) drawSplit(v, s, c, q, e, shown);
        continue;
      }
      locate(v, s, k);
      const id = gt ? gt.value(v.loc) : -1;
      glyphs.begin(inkOf(key), keyAlpha(key));
      const g = glyphs.add(id < 0 ? 0 : id);
      glyphs.x[g] = x * dpr;
      glyphs.y[g] = s.frame[0] * dpr;
      const size = c.radius ? c.radius[k] : NO_SIZE;
      glyphs.radius[g] = (size < 0 ? R : Math.min(8, Math.max(1, size / SIZE_ONE))) * dpr;
    }
  }
}

/** Pair k's stroke width in half device px: its width theme's, rounded (the batch's quantum), else `dflt`. */
function halfWidth(c: PairColors, k: number, dpr: number, dflt: number): number {
  const size = c.width ? c.width[k] : NO_SIZE;
  return size < 0 ? dflt : Math.max(1, Math.min(255, Math.round((size * dpr) / (SIZE_ONE / 2))));
}
let half = new Uint8Array(64);
let byWidth = new Int32Array(64);
const widths = new Int32Array(257);
/**
 * The pairs in view into `sorted`, by their color's rank (a counting sort, stable), and within a color, with a width
 * theme, the widest first (a first counting sort by width, so thin arcs lie over wide ones and a backend strokes each
 * (color, width) once): their count; `half[k]` holds each one's width.
 */
function sortByColor(v: DrawView, s: PairsState, c: PairColors): number {
  if (list.length < s.count) {
    list = new Int32Array(s.count);
    sorted = new Int32Array(s.count);
    byWidth = new Int32Array(s.count);
    half = new Uint8Array(s.count);
  }
  if (byRank.length < c.ranks + 1) byRank = new Int32Array(c.ranks + 1);
  const n = pairsInView(v, s, list);
  let from = list;
  if (c.width) {
    const dflt = Math.max(1, Math.round(2 * v.dpr));
    widths.fill(0);
    for (let j = 0; j < n; j++) half[list[j]] = halfWidth(c, list[j], v.dpr, dflt);
    for (let j = 0; j < n; j++) widths[256 - half[list[j]]]++;
    for (let w = 1; w <= 256; w++) widths[w] += widths[w - 1];
    for (let j = 0; j < n; j++) byWidth[widths[255 - half[list[j]]]++] = list[j];
    from = byWidth;
  }
  byRank.fill(0, 0, c.ranks + 1);
  for (let j = 0; j < n; j++) byRank[c.rank[from[j]] + 1]++;
  for (let r = 1; r <= c.ranks; r++) byRank[r] += byRank[r - 1];
  for (let j = 0; j < n; j++) sorted[byRank[c.rank[from[j]]]++] = from[j];
  return n;
}

function hand(v: DrawView) {
  if (bars.n > 0) v.backend.rects(bars);
  if (glyphs.n > 0) {
    glyphs.font = v.text.font("ui", 12, 600);
    v.backend.glyphs(glyphs);
  }
  if (markers.n > 0) {
    markers.style(v.text.font("ui", MARKER_PX, 600), TextAlign.center, TextBaseline.middle);
    v.backend.text(markers);
  }
  if (arcs.n > 0) v.backend.paths(arcs);
  if (wedges.n > 0) v.backend.paths(wedges);
}

/** Ends as glyphs (density below a device px a cell), then every arc in view with `arcs: "always"`. */
/** The states last drawn on rings, as chords (`pushArc`): their arcs do not hit. */
const chords = new WeakSet<PairsState>();
function draw(v: DrawView, s: PairsState, p: Params) {
  if (v.ring) chords.add(s);
  else chords.delete(s);
  if (v.geometry !== AxisGeometry.linear || !v.axis || s.count === 0) return;
  const c = colorsOf(v, s);
  arcs.reset();
  wedges.reset();
  glyphs.reset();
  markers.reset();
  bars.reset();
  Arcs.frame(v, s, p);
  // On rings a chord's ends stay glyphs while there are few pairs (a plasmid map's), so each chord is seen to land.
  const dense = v.cell * v.dpr < 1 && !(v.ring && s.count <= RING_GLYPHS);
  if (p.display !== "arcs") {
    if (dense) densityStart(v);
    drawEnds(v, s, c, dense);
    const y = s.frame[0];
    const y0 = Math.round((y - R) * v.dpr) | 0;
    if (dense) densityDraw(bars, y0, Math.round((y + R) * v.dpr) | 0, DensityKind.stacked, SkinToken.inkQuiet);
  }
  if (p.display !== "ends" && p.arcs !== "hover") {
    const n = sortByColor(v, s, c);
    const w = Math.max(1, Math.round(2 * v.dpr));
    for (let j = 0; j < n; j++) {
      const k = sorted[j];
      pushArc(v, s, k, c.width ? half[k] : w, c.key[k], arcs, wedges, c.height ? c.height[k] : NO_SIZE);
    }
    hand(v);
    if (p.labels === true && !v.ring) drawLabels(v, s, c, p, n);
    return;
  }
  hand(v);
}

/** Per state, its pairs' labels (`labels: true`): a row's `label`, else its width theme's field's value. Edge-rate. */
const pairLabels = new WeakMap<PairsState, { field: string | null; text: readonly (string | null)[] }>();
function labelsOf(s: PairsState, p: Params): readonly (string | null)[] {
  const w = p.width as { theme?: string; params?: { field?: unknown } } | null | undefined;
  const field = w?.theme === "field_scale" && typeof w.params?.field === "string" ? w.params.field : null;
  const known = pairLabels.get(s);
  if (known && known.field === field) return known.text;
  const text = Array.from({ length: s.count }, (_, k) => {
    const a = s.set.annotations[s.annotation[k]];
    const x = a
      ? (Annotations.rowValue(a, s.set, s.row[k], "label") ??
        (field ? Annotations.rowValue(a, s.set, s.row[k], field) : null))
      : null;
    return x == null || x === "" ? null : String(x);
  });
  pairLabels.set(s, { field, text });
  return text;
}
/** A pair label's size, CSS px. */
const LABEL_PX = 10;
const names = new TextBatch(32);
/**
 * Each drawn arc's label (`labels: true`; seqspec:format#layer-types `pairs`) inside it under its apex, where it fits:
 * the arc taller than the text and wider than it. After the arcs (`sorted`'s first `n`), so they never cover it.
 */
function drawLabels(v: DrawView, s: PairsState, c: PairColors, p: Params, n: number) {
  const text = labelsOf(s, p);
  const font = v.text.font("ui", LABEL_PX, 600);
  const f = s.frame;
  const dpr = v.dpr;
  names.reset();
  names.style(font, TextAlign.center, TextBaseline.middle);
  for (let j = 0; j < n; j++) {
    const k = sorted[j];
    const t = text[k];
    if (t === null || !Arcs.ends(v, s, k)) continue;
    Arcs.apex(s, c.height ? c.height[k] : NO_SIZE);
    const mid = (s.xs[0] + s.xs[1]) / 2;
    // Only an apex in view, of an arc over its real ends (none taken in at a margin: `Arcs.virtual`). Below the apex by
    // half the text and the stroke's half width; its bottom `lift` px over the ends' line, where the arc (a parabola)
    // is `span · sqrt(1 − lift / apex)` wide: the text fits there, inside the stroke, 2 px a side.
    const span = s.xs[3] - s.xs[2];
    if (mid < v.left || mid > v.right || span !== Math.abs(s.xs[1] - s.xs[0])) continue;
    const stroke = c.width ? half[k] / (4 * dpr) : 1;
    const below = LABEL_PX / 2 + 2 + stroke;
    const lift = f[9] - below - LABEL_PX / 2;
    if (lift < 1) continue;
    const room = span * Math.sqrt(1 - lift / f[9]) - 2 * stroke - 4;
    if (v.text.label(t, font).width > room) continue;
    names.begin(SkinToken.ink);
    const i = names.add(t);
    names.x[i] = mid * dpr;
    names.y[i] = Math.round((f[0] + f[1] * (f[9] - below)) * dpr);
  }
  if (names.n > 0) v.backend.text(names);
}

/** Below the row's minimum: the ends in view as density over the band, from the same end index. */
function tiny(v: DrawView, s: PairsState, _p: Params) {
  if (v.geometry !== AxisGeometry.linear || !v.axis || s.count === 0) return;
  const c = colorsOf(v, s);
  bars.reset();
  glyphs.reset();
  markers.reset();
  arcs.reset();
  wedges.reset();
  densityStart(v);
  drawEnds(v, s, c, true);
  const y0 = Math.round(v.top * v.dpr) | 0;
  densityDraw(bars, y0, Math.round((v.top + v.height) * v.dpr) | 0, DensityKind.stacked, SkinToken.inkQuiet);
  hand(v);
}

const LIT = MarkBit.selection | MarkBit.highlight | MarkBit.accent | MarkBit.hover;
/** Per state, the pairs the overlay lit last (room for every pair, its count last): `drawn` is cleared through it. */
const lit = new WeakMap<PairsState, Int32Array>();
const litOf = (s: PairsState): Int32Array => lit.get(s) ?? newLit(s);
function newLit(s: PairsState): Int32Array {
  const x = new Int32Array(s.count + 1);
  lit.set(s, x);
  return x;
}

/** Per state, its pairs by item: the item ids sorted, then the pair of each (edge rate: built on its first item hover). */
const byItem = new WeakMap<PairsState, Int32Array>();
function newByItem(s: PairsState): Int32Array {
  const id = (k: number) => s.set.annotations[s.annotation[k]]?.item[s.row[k]] ?? -1;
  const ks = Array.from({ length: s.count }, (_, k) => k).sort((a, b) => id(a) - id(b));
  const x = new Int32Array(2 * s.count);
  for (let j = 0; j < s.count; j++) {
    x[j] = id(ks[j]);
    x[s.count + j] = ks[j];
  }
  byItem.set(s, x);
  return x;
}
/** The pair of item `item` of annotation `annotation` (-1: of any); -1: none. Ids repeat across annotations only. */
function pairOf(s: PairsState, annotation: number, item: number): number {
  const x = byItem.get(s) ?? newByItem(s);
  let a = 0;
  let b = s.count;
  while (a < b) {
    const m = (a + b) >> 1;
    if (x[m] < item) a = m + 1;
    else b = m;
  }
  for (; a < s.count && x[a] === item; a++)
    if (annotation < 0 || s.annotation[x[s.count + a]] === annotation) return x[s.count + a];
  return -1;
}

/** The cell the hover is on (`elementAt`: the hit's block and index), as the end of pair `k` that holds it; -1: none. */
function cellOfHover(s: PairsState, k: number): number {
  for (let side = 0; side < 2; side++) {
    const e = 2 * k + side;
    if (s.block[e] === elementAt[0] && s.start[e] <= elementAt[1] && elementAt[1] < s.stop[e] && s.endAt[e] >= 0)
      return s.endCell[s.endAt[e]];
  }
  return -1;
}
/** Whether pair `k` is of kind `kind` and has an end in cell `cell`. */
function inGroup(s: PairsState, k: number, cell: number, kind: number): boolean {
  if (s.kind[k] !== kind) return false;
  const a = s.endAt[2 * k];
  const b = s.endAt[2 * k + 1];
  return (a >= 0 && s.endCell[a] === cell) || (b >= 0 && s.endCell[b] === cell);
}

/**
 * The arcs of lit ends on the overlay (a hover, a highlight, the selection), a hovered pair of the layer alone for the
 * hover (`marks.item` of `marks.itemAnnotation`, last, on top; not the pairs sharing its ends): the same curves as the draw, 1.5 px with
 * `arcs: "hover"` (the only arcs) else 2.5 px, each recorded in `s.drawn` by its place (with `arcs: "hover"` only
 * they hit; the last drawn is on top).
 */
function overlay(v: DrawView, s: PairsState, p: Params) {
  if (v.geometry !== AxisGeometry.linear || !v.axis || s.count === 0 || p.display === "ends") return;
  const hover = p.arcs === "hover";
  const ids = litOf(s);
  for (let j = 0; j < ids[s.count]; j++) s.drawn[ids[j]] = 0;
  const item = v.marks.item(v.layer);
  const hovered = item >= 0 ? pairOf(s, v.marks.itemAnnotation(v.layer), item) : -1;
  let n = litPairs(v, s, v.marks, hovered >= 0 ? LIT & ~MarkBit.hover : LIT, ids);
  if (hovered >= 0) {
    // A glyph of a split cell (the hit's `part` is its kind): every pair of that kind at that end, else the pair alone.
    const cell = v.marks.element(v.layer, elementAt) && elementAt[3] >= 0 ? cellOfHover(s, hovered) : -1;
    const kind = elementAt[3];
    let m = 0;
    for (let j = 0; j < n; j++)
      if (ids[j] !== hovered && !(cell >= 0 && inGroup(s, ids[j], cell, kind))) ids[m++] = ids[j];
    if (cell >= 0) {
      for (let q = s.cellLo[cell]; q < s.cellLo[cell + 1]; q++) {
        const e = s.endOf[q];
        const k = e >> 1;
        // A pair with both ends in the cell once.
        if (s.kind[k] === kind && !((e & 1) === 1 && s.endCell[s.endAt[e - 1]] === cell)) ids[m++] = k;
      }
      // The hovered pair last, on top.
      for (let j = 0; j < m; j++) {
        if (ids[j] !== hovered) continue;
        ids[j] = ids[m - 1];
        ids[m - 1] = hovered;
        break;
      }
    } else ids[m++] = hovered;
    n = m;
  }
  ids[s.count] = n;
  if (n === 0) return;
  const c = colorsOf(v, s);
  arcs.reset();
  wedges.reset();
  bars.reset();
  glyphs.reset();
  markers.reset();
  ThemeLocation.reset(v.loc);
  Arcs.frame(v, s, p);
  // Half device px: `HOVER_WIDTH` (the hover outline's, seqspec:interaction#hover) with `arcs: "hover"`, else 2.5 px.
  const w = Math.max(1, Math.round((hover ? 2 * HOVER_WIDTH : 5) * v.dpr)) | 0;
  for (let j = 0; j < n; j++) {
    // A wide arc is lit over its own width.
    const own = c.width ? halfWidth(c, ids[j], v.dpr, 0) : 0;
    const h = c.height ? c.height[ids[j]] : NO_SIZE;
    if (pushArc(v, s, ids[j], Math.max(w, own + 2), c.key[ids[j]], arcs, wedges, h)) s.drawn[ids[j]] = j + 1;
  }
  hand(v);
}

/** Linear only: the format's `pairs` is a tracks type; in the text layout its fallback, `points`, draws the ends. */
export const PairsRenderer = kernelEntry<PairsState>("pairs", 48, {
  miniature: false,
  hit: (v, s, p, out) => hit(v, s, p as Params, out),
  draw: (v, s, p) => draw(v, s, p as Params),
  tiny: (v, s, p) => tiny(v, s, p as Params),
  overlay: (v, s, p) => overlay(v, s, p as Params),
  sized: ["width", "height", "radius"],
});
