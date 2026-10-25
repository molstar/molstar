/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `secondary` renderer (seqspec:viewer#extension-points): helices as cylinders (their ends rounded where the helix
// ends, square where a line, the view or an intron cuts it), strands as arrows and turns over a coil line; below 3 px
// a feature (and on the tiny path) a bar; with `labels`, each feature's label over its shape (seqspec:format#layer-types);
// wrapped, underlines on the lanes of the line. Placed through `placeRun`, in device px. The overlay outlines the
// hovered element inside its own pixels (seqspec:interaction#hover).
import {
  AxisGeometry,
  type DrawView,
  drawLineLanes,
  type Features,
  featureHoverKey,
  featureInk,
  featureKey,
  featureSize,
  fitIn,
  fitted,
  HoverOutline,
  hitLineLanes,
  hoverEmphasis,
  hoveredAt,
  hoveredItem,
  hoverFill,
  hoverInk,
  hoverWidth,
  kernelEntry,
  keyAlpha,
  outlineLineLanes,
  PathBatch,
  PathMode,
  PathVerb,
  placeRun,
  RectBatch,
  runHull,
  SecondaryKind,
  SkinToken,
  segWindow,
  TextAlign,
  TextBaseline,
  TextBatch,
  shownHull,
  sizeTheme,
  ThemeLocation,
  visibleRange,
} from '../draw/api.js';

// Read once: a module member read per frame allocates under the test runner's module wrapper.
const STRAND = SecondaryKind.strand;
const HELIX = SecondaryKind.helix;
const OTHER = SecondaryKind.other;
const WRAPPED = AxisGeometry.wrapped;
/** The labels (`labels`), one font a draw. */
const names = new TextBatch(16);
const fitArgs = new Float64Array(2);
/** A label's extent (device px): the feature's pieces drawn. */
const ext = new Int32Array(2);
const helixHull = new Int32Array(2);
const INK = SkinToken.ink;
const rects = new RectBatch(64);
const shapes = new PathBatch(32);
const range = new Int32Array(2);
const hull = new Int32Array(2);
const win = new Int32Array(2);
const piece = new Float64Array(64);
/**
 * The opaque features drawn at the current start (a block's lane is sorted by start): index, color key, half height.
 * The copies of a merged block often carry one feature each alike (3J3Q: 18,965 rows, 253 distinct), and an opaque shape
 * over its equal changes no pixel, so it is skipped; past `SAME` variants at one start, the rest are drawn.
 */
const SAME = 16;
const same = new Int32Array(3 * SAME);
let sameCount = 0;
let sameStart = -1;

/** Whether feature k (key, half height hk) is one already drawn at its start; else noted, while there is room. */
function drawnAlike(f: Features, k: number, key: number, hk: number, last: Uint8Array): boolean {
  if (f.start[k] !== sameStart) {
    sameStart = f.start[k];
    sameCount = 0;
  }
  for (let i = 0; i < sameCount; i++) {
    const j = same[3 * i];
    if (
      same[3 * i + 1] === key &&
      same[3 * i + 2] === hk &&
      f.stop[j] === f.stop[k] &&
      f.kind[j] === f.kind[k] &&
      last[j] === last[k] &&
      f.label[j] === f.label[k]
    )
      return true;
  }
  if (sameCount < SAME) {
    same[3 * sameCount] = k;
    same[3 * sameCount + 1] = key;
    same[3 * sameCount + 2] = hk;
    sameCount++;
  }
  return false;
}

/** A vertex of the shape being built (device px, integers: no double crosses the call). */
function vertex(first: boolean, x: number, y: number) {
  const o = shapes.verb(first ? PathVerb.move : PathVerb.line);
  shapes.coords[o] = x;
  shapes.coords[o + 1] = y;
}
function curve(cx: number, cy: number, x: number, y: number) {
  const o = shapes.verb(PathVerb.quad);
  shapes.coords[o] = cx;
  shapes.coords[o + 1] = cy;
  shapes.coords[o + 2] = x;
  shapes.coords[o + 3] = y;
}

/**
 * One piece [x0, x1) of a feature of `kind`, on the coil line at cy, half-height h (all device px); a strand's head on
 * the piece's right (`head` 1), its left (-1) or none (0).
 */
function shape(kind: number, x0: number, x1: number, cy: number, h: number, ink: number, alpha: number, head = 1) {
  // A helix's ends: `head` bit 1 its left rounded, bit 2 its right (`helixEnds`).
  if (kind === HELIX) cylinder(x0, x1, cy, h, ink, alpha, head);
  else if (kind === SecondaryKind.strand) {
    shapes.begin(ink, alpha, PathMode.fill);
    const hl = head === 0 ? 0 : Math.min(2 * h, (x1 - x0) >> 1);
    const body = Math.round(h * 0.55) | 0;
    // Mirrored for a head on the left: from the tip's side.
    const a = head < 0 ? x1 : x0;
    const z = head < 0 ? x0 : x1;
    const n = head < 0 ? z + hl : z - hl;
    vertex(true, a, cy - body);
    vertex(false, n, cy - body);
    if (hl > 0) {
      vertex(false, n, cy - h);
      vertex(false, z, cy);
      vertex(false, n, cy + h);
    }
    vertex(false, n, cy + body);
    vertex(false, a, cy + body);
    shapes.close();
  } else if (kind === SecondaryKind.turn) {
    shapes.begin(ink, alpha, PathMode.fill);
    const top = cy - (Math.round(2.4 * h) | 0);
    vertex(true, x0, cy);
    curve((x0 + x1) >> 1, top, x1, cy);
    vertex(false, x1, cy + 1);
    curve((x0 + x1) >> 1, top + 2, x0, cy + 1);
  } else {
    const e = Math.max(1, Math.round(h * 0.45) | 0);
    rects.push(x0, cy - e, x1 - x0, 2 * e);
  }
}

/**
 * A helix's piece [x0, x1) as a cylinder: its body, each end rounded (a quarter ellipse above and below the axis, as
 * wide as 0.6 of its half height, at most half the piece) where `ends` says the helix ends there (bit 1 left, 2 right),
 * else square. Device px, integers.
 */
function cylinder(x0: number, x1: number, cy: number, h: number, ink: number, alpha: number, ends: number) {
  const rx = Math.min(Math.max(1, Math.round(h * 0.6)) | 0, (x1 - x0) >> 1);
  const l = ends & 1 ? rx : 0;
  const r = ends & 2 ? rx : 0;
  shapes.begin(ink, alpha, PathMode.fill);
  vertex(true, x0 + l, cy - h);
  vertex(false, x1 - r, cy - h);
  if (r > 0) {
    curve(x1, cy - h, x1, cy);
    curve(x1, cy + h, x1 - r, cy + h);
  } else vertex(false, x1, cy + h);
  vertex(false, x0 + l, cy + h);
  if (l > 0) {
    curve(x0, cy + h, x0, cy);
    curve(x0, cy - h, x0 + l, cy - h);
  } else vertex(false, x0, cy - h);
  shapes.close();
}

/**
 * Which ends of piece [x0, x1) of helix k (block b) are the helix's own (bit 1 left, 2 right): where the piece's edge is
 * its first or last shown cell's (a reversed piece too: shown order is screen order); a line's, the view's or an
 * intron's cut stays square. Ints in and out.
 */
function helixEnds(v: DrawView, f: Features, k: number, b: number, x0: number, x1: number): number {
  if (!runHull(v, b, f.start[k], Math.max(f.stop[k], f.start[k] + 1), helixHull)) return 0;
  const lo = helixHull[0];
  const hi = helixHull[1];
  const right = v.geometry === WRAPPED ? v.xe(b, hi - 1) : v.xd(b, hi);
  return (Math.abs(x0 - v.xd(b, lo)) <= 1 ? 1 : 0) | (Math.abs(x1 - right) <= 1 ? 2 : 0);
}

/** Feature k's label over its pieces' extent (`ext`), in the label row [y0, y1) (device px), where it fits. */
function label(v: DrawView, f: Features, k: number, y: number) {
  const text = f.label[k];
  if (text === null || ext[1] <= ext[0]) return;
  const p = v.text.label(text, names.font);
  fitArgs[0] = (ext[1] - ext[0]) / v.dpr;
  const m = fitIn(p, fitArgs);
  // Whole, or cut down to three graphemes at least (a short label, α1, is never cut to nothing).
  if (m === 0 || (m < p.ends.length && m < 3)) return;
  const t = names.add(fitted(p, m));
  names.x[t] = (ext[0] + ext[1]) >> 1;
  names.y[t] = y;
}

function draw(v: DrawView, f: Features, tiny: boolean, labeled: boolean) {
  const ax = v.axis;
  rects.reset();
  shapes.reset();
  names.reset();
  ThemeLocation.reset(v.loc);
  if (v.geometry === AxisGeometry.wrapped && !v.trackForm) drawLineLanes(v, f, rects);
  else if (ax) {
    const dpr = v.dpr;
    const labels = shapeBand(v, labeled, tiny);
    const cy = band[0];
    const h = band[1];
    const top = v.top;
    const height = v.height;
    if (labels) {
      names.style(
        v.text.font('ui', Math.max(7, Math.min(11, Math.floor(height * 0.4)))),
        TextAlign.center,
        TextBaseline.bottom,
      );
      names.begin(INK);
      fitArgs[1] = v.text.label('…', names.font).width;
    }
    const labelY = Math.round((top + height * 0.45) * dpr) | 0;
    const thickness = sizeTheme(v, 0, 1);
    const last = lastFragments(v, f);
    // The coil line across each block's view, under the features.
    rects.begin(SkinToken.inkFaint);
    for (let b = 0; b < ax.count; b++) {
      if (v.visHi[b] <= v.visLo[b]) continue;
      // Wrapped (a track drawn in full), to the line's last cell: past it a group's gap or the next line.
      const end = v.geometry === AxisGeometry.wrapped ? v.xe(b, v.visHi[b] - 1) : v.xd(b, v.visHi[b]);
      rects.push(v.xd(b, v.visLo[b]), cy, end - v.xd(b, v.visLo[b]), 1);
    }
    for (let b = 0; b < ax.count && (b + 1) * f.lanes < f.seg.length; b++) {
      visibleRange(v, b, range);
      if (range[0] >= range[1]) continue;
      shownHull(v, b, range[0], range[1], hull);
      for (let lane = 0; lane < f.lanes; lane++) {
        const g = b * f.lanes + lane;
        segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
        sameStart = -1;
        for (let k = win[0]; k < win[1]; k++) {
          const stop = Math.max(f.stop[k], f.start[k] + 1);
          if (stop <= hull[0] || (f.on[k] >= 0 && f.on[k] !== v.row)) continue;
          const key = featureKey(v, f, k);
          const alpha = keyAlpha(key);
          const hk = thick(h, featureSize(v, f, k, thickness));
          if (alpha === 255 && drawnAlike(f, k, key, hk, last)) continue;
          const n = placeRun(v, b, f.start[k], stop, piece);
          const ink = featureInk(key);
          rects.begin(ink, alpha);
          ext[0] = 0x7fffffff;
          ext[1] = -0x7fffffff;
          for (let q = 0; q < n; q++) {
            const x0 = piece[2 * q] | 0;
            const x1 = Math.max(x0 + 1, piece[2 * q + 1] | 0);
            const kind = tiny || x1 - x0 < 3 * dpr ? OTHER : f.kind[k];
            const head =
              kind === STRAND ? headAt(v, f, k, b, x0, x1) : kind === HELIX ? helixEnds(v, f, k, b, x0, x1) : 0;
            shape(kind, x0, x1, cy, hk, ink, alpha, head);
            if (x0 < ext[0]) ext[0] = x0;
            if (x1 > ext[1]) ext[1] = x1;
          }
          if (labels) label(v, f, k, labelY);
        }
      }
    }
  }
  if (rects.n > 0) v.backend.rects(rects);
  if (shapes.n > 0) v.backend.paths(shapes);
  if (names.n > 0) v.backend.text(names);
}

const lasts = new WeakMap<Features, Uint8Array>();
/**
 * Per feature, 1 where it is its item's last fragment on its block, the C-terminal one: of its item's rows (rows sharing
 * an id are one item, LC-04; a range across a cyclic origin is its two runs, head then tail) the last, and of that row's
 * fragments the one reaching furthest (an exon of a strand split by an intron is a fragment each). Edge rate: once per
 * prepared state.
 */
function lastFragments(v: DrawView, f: Features): Uint8Array {
  let m = lasts.get(f);
  if (m) return m;
  m = new Uint8Array(f.count);
  const set = v.axis?.set;
  const last = new Map<string, number>();
  for (let k = 0; k < f.count; k++) {
    const item = set?.annotations[f.annotation[k]]?.item[f.row[k]] ?? f.row[k];
    const key = `${f.block[k]} ${f.annotation[k]} ${item}`;
    const j = last.get(key);
    if (j === undefined || f.row[k] > f.row[j] || (f.row[k] === f.row[j] && f.stop[k] > f.stop[j])) last.set(key, k);
  }
  for (const k of last.values()) m[k] = 1;
  lasts.set(f, m);
  return m;
}
const headRuns = new Int32Array(16);
/** Block b's shown position of the highest index of [start, stop) shown (LC-05: the last one may be filtered out); -1. */
function lastShown(v: DrawView, b: number, start: number, stop: number): number {
  const n = Math.min(headRuns.length >> 1, v.runToShown(b, start, stop, headRuns));
  let best = -1;
  let at = -1;
  for (let r = 0; r < n; r++)
    for (let e = 0; e < 2; e++) {
      const q = e === 0 ? headRuns[2 * r] : headRuns[2 * r + 1] - 1;
      const i = v.indexOf(b, q as never);
      if (i > best) {
        best = i;
        at = q;
      }
    }
  return at;
}
/**
 * Which side of piece [x0, x1) of feature k (block b) holds a strand's head: its C-terminal end's outer edge on screen
 * (a reversed piece's left), the last of its positions shown, on the piece that reaches it, of its item's last fragment;
 * 0: none there (a piece or a fragment before it, or none of it shown).
 * Ints in and out.
 */
function headAt(v: DrawView, f: Features, k: number, b: number, x0: number, x1: number): number {
  if (lastFragments(v, f)[k] === 0) return 0;
  const qe = lastShown(v, b, f.start[k], Math.max(f.stop[k], f.start[k] + 1));
  if (qe < 0) return 0;
  if (v.dir[b] > 0) return Math.abs(x1 - v.xd(b, qe + 1)) <= 1 ? 1 : 0;
  return Math.abs(x0 - v.xd(b, qe)) <= 1 ? -1 : 0;
}

/** The coil line and a shape's half height (device px). */
const centerY = (v: DrawView) => Math.round((v.top + v.height / 2) * v.dpr) | 0;
const halfHeight = (v: DrawView) => Math.max(1, Math.round(Math.max(1, v.height * 0.31) * v.dpr)) | 0;
const band = new Int32Array(2);
/**
 * The shapes' coil line and half height (device px) into `band`, the draw's and the hover outline's: with labels, the
 * shapes in the band's lower 55%, the labels over them (seqspec:format#layer-types). Whether the labels show.
 */
function shapeBand(v: DrawView, labeled: boolean, tiny: boolean): boolean {
  const labels = labeled && !tiny && v.height >= 14;
  const top = v.top;
  const height = v.height;
  if (labels) {
    v.top = top + height * 0.45;
    v.height = height * 0.55;
  }
  band[0] = centerY(v);
  band[1] = halfHeight(v);
  v.top = top;
  v.height = height;
  return labels;
}
/** The thickness channel: half height `h` (device px) at a share of the shape's (a size theme's int; < 0: whole). */
const thick = (h: number, size: number) => (size < 0 ? h : Math.max(1, (h * Math.min(size, 256) + 128) >> 8));

const hovered = new Int32Array(2);
/**
 * The hovered element inside its own pixels: a helix's or a bar's box outlined inside; a strand's arrow in the hover
 * ink with the arrow `hoverWidth` smaller in its fill over it; a turn's thin arch in the ink; wrapped, its underline.
 */
function overlay(v: DrawView, f: Features, labeled: boolean) {
  if (!hoverEmphasis(v) || !v.axis) return;
  if (v.geometry === AxisGeometry.wrapped && !v.trackForm) {
    outlineLineLanes(v, f);
    return;
  }
  const x = v.geometry === AxisGeometry.linear ? hoveredItem(v, f, hovered) : null;
  if (x === null) return;
  const tiny = v.height < (SecondaryRenderer.minHeight ?? 0);
  const dpr = v.dpr;
  const d = hoverWidth(dpr);
  shapeBand(v, labeled, tiny);
  const cy = band[0];
  const h0 = band[1];
  const thickness = sizeTheme(v, 0, 1);
  rects.reset();
  shapes.reset();
  HoverOutline.reset();
  for (let j = hovered[0]; j < hovered[1]; j++) {
    const k = hoveredAt(v, f, x, j);
    if (k < 0) continue;
    const n = placeRun(v, f.block[k], f.start[k], Math.max(f.stop[k], f.start[k] + 1), piece);
    const key = featureHoverKey(v, f, k);
    const h = thick(h0, featureSize(v, f, k, thickness));
    for (let q = 0; q < n; q++) {
      const x0 = piece[2 * q] | 0;
      const x1 = Math.max(x0 + 1, piece[2 * q + 1] | 0);
      const kind = tiny || x1 - x0 < 3 * dpr ? OTHER : f.kind[k];
      if (kind === HELIX) HoverOutline.box(v, key, x0, cy - h, x1 - x0, 2 * h);
      else if (kind === SecondaryKind.other) {
        const e = Math.max(1, Math.round(h * 0.45) | 0);
        HoverOutline.box(v, key, x0, cy - e, x1 - x0, 2 * e);
      } else {
        const head = kind === STRAND ? headAt(v, f, k, f.block[k], x0, x1) : 0;
        shape(kind, x0, x1, cy, h, hoverInk(v.skin, key), 255, head);
        if (kind === STRAND && x1 - x0 > 4 * d && h > 2 * d)
          shape(kind, x0 + d, x1 - d, cy, h - d, hoverFill(v.skin, key), 255, head);
      }
    }
  }
  if (shapes.n > 0) v.backend.paths(shapes);
  HoverOutline.hand(v);
}

export const SecondaryRenderer = kernelEntry<Features>('secondary', 18, {
  geometries: ['linear', 'wrapped'],
  draw: (v, f, p) => draw(v, f, false, (p as { labels?: boolean }).labels === true),
  tiny: (v, f) => draw(v, f, true, false),
  overlay: (v, f, p) => overlay(v, f, (p as { labels?: boolean }).labels === true),
  // Wrapped, its underlines on the line's lanes; linear, the position's item (hit.ts's lookup).
  hit: (v, f, _p, out) => hitLineLanes(v, f, out),
  sized: ['thickness'],
});
