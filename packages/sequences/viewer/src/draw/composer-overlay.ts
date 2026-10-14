/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's overlay and navigator passes (draw/composer.ts): what is lit, redrawn whole (row marks and the
// hovered row under the position marks, streamed layers' progress, renderers' overlays, the cursor's ring, the ring's
// window), then the rulerRail and the navigator's window; and the navigator's row. Per frame nothing is allocated
// (RN1's draw pattern).

import { type Axis, AxisGeometry, type SectionAxisView } from "../axis.ts";
import { SeqViewerConfig } from "../config.ts";
import type { ViewModel } from "../graph.ts";
import { HitRegion, HoverKind } from "../hit.ts";
import type { HoverState } from "../hover.ts";
import { cursorRect } from "../keyboard.ts";
import { EntryKind, type SectionLayout, SectionLayouts } from "../layout.ts";
import { MarkBit } from "../marks.ts";
import type { NavigatorEntry } from "../navigators.ts";
import { RW_CYCLIC, RW_HI, RW_LO, RW_X0, RW_X1, ringBounds, ringWindowFrame } from "../ring-window.ts";
import type { ViewportState } from "../viewport.ts";
import { gripsLit, litPart, railGrown, WindowKind, windowFrame, windowIdle } from "../windows.ts";
import { Pass, type RenderBackend, SkinToken } from "./backend.ts";
import {
  banded,
  boxOf,
  type ComposerFrame,
  fail,
  gutterWidth,
  index,
  name,
  next,
  open,
  place,
  prime,
  release,
  setClip,
  stateOf,
  tint,
  view,
} from "./composer-frame.ts";
import { pending, positionMarks, rowMarks } from "./composer-marks.ts";
import { navDraw, navDrawWindow, WINDOW_EDGE_INK, WINDOW_FILL } from "./nav.ts";
import { Moving, progress1D, progressBegin, progressEnd, progressFlush, progressRows } from "./progress.ts";

// The ring window's inks, read once (a read per frame allocated under the test runner's module wrapper).
const RW_RIM = SkinToken.inkQuiet;
const RW_HALO = SkinToken.surface;
const RW_GRIP = SkinToken.grip;
const RW_GRIP_LIT = SkinToken.gripHover;
const RW_KIND = WindowKind.ring;

// The other enum-like members, read once alike.
const CIRCULAR = AxisGeometry.circular;
const LINEAR = AxisGeometry.linear;
const TRACK = EntryKind.track;
const COLUMN = HitRegion.column;
const DATA = HitRegion.data;
const ALIGNED = HoverKind.aligned;
const HEADER = HoverKind.header;
const ROW = HoverKind.row;
const BIT_HOVER = MarkBit.hover;
const FADE = Moving.fade;
const PASS_NAV = Pass.nav;
const PASS_NAV_OVERLAY = Pass.navOverlay;
const PASS_OVERLAY = Pass.overlay;
const SHIMMER_FPS = SeqViewerConfig.Stream.ShimmerFps;
const REDUCED_MOTION = SeqViewerConfig.Viewport.ReducedMotion;
const TOKEN_FOCUS_RING = SkinToken.focusRing;
const NAVIGATOR = WindowKind.navigator;

/** The navigator's axis in its row (`navDraw`: its entry's, else from its parts); redrawn on a size, layout or skin change. */
export function nav(
  f: ComposerFrame,
  be: RenderBackend,
  l: SectionLayout,
  av: SectionAxisView,
  s: ViewportState,
  version: number,
) {
  const a = f.navAt;
  const e = l.navHeight > 0 ? navOf(f, be, l, av, s) : null;
  if (!e || e.parts.onRuler) return;
  // The map's image (its key: layout, regions, fonts, model, skin) may change with no layout or size change; its size
  // is the band's, so the row is set first (the last pass left the view as it was: a resize and a fresh mount differed).
  f.nview.top = 0;
  f.nview.height = l.navHeight;
  const map = f.mini.image(f.nview);
  if (a[0] === version && a[1] === l.width && a[2] === av.dpr && map === f.nview.map) return;
  a[0] = version;
  a[1] = l.width;
  a[2] = av.dpr;
  f.nview.map = map;
  navBand(f, PASS_NAV, l, av);
  navDraw(e, f.nview);
  be.endPass();
}
/**
 * What is lit, redrawn whole: row marks and the hovered row under the position marks, renderers' overlays, the
 * cursor's ring; then the rulerRail and the navigator's window.
 */
export function overlay(
  f: ComposerFrame,
  be: RenderBackend,
  l: SectionLayout,
  av: SectionAxisView,
  s: ViewportState,
  version: number,
) {
  const W = Math.ceil(l.width * av.dpr) | 0;
  const H = Math.ceil(l.height * av.dpr) | 0;
  be.beginPass(PASS_OVERLAY, W, H, av.dpr);
  be.clear(0, 0, W, H);
  const main = f.v;
  const m = f.marks;
  f.v = f.ov;
  try {
    const v = f.v;
    view(f, v, be, av);
    v.lanes.set(main.lanes);
    if (m) v.marks = m;
    setClip(f, 0, W, 0, H);
    f.fills.reset();
    f.edges.reset();
    f.clade = -1;
    progress(f, l, s);
    if (m) rowMarks(f, l, m, version);
    rowBand(f, l, s);
    flushMarks(f);
    if (f.clade >= 0) cladeHover(f, l, f.clade);
    f.bottom = l.navHeight;
    for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) if (boxOf(f, l, e)) f.bottom = Math.max(f.bottom, f.box[3]);
    if (m) positionMarks(f, l, s, m);
    pending(f, l, s);
    flushMarks(f);
    overlays(f, l, s);
    cursorRing(f);
    if (l.geometry === CIRCULAR) ringWindow(f, l);
  } finally {
    f.v = main;
  }
  const e = navOf(f, be, l, av, s);
  const n = f.nview;
  const re = l.railEntry;
  if (e?.parts.onRuler && re >= 0) {
    // Along the pinned ruler's bottom border; its grips as near the rail as the ruler allows, never over the rows.
    const bottom = Math.round((l.entryTop[re] + l.entryHeight[re]) * av.dpr) / av.dpr;
    n.top = bottom - 8;
    n.height = 8;
    n.gripTop = l.entryTop[re];
    n.gripBottom = bottom;
    navDraw(e, n);
    navDrawWindow(e, n);
  }
  be.endPass();
  if (!e || e.parts.onRuler || l.navHeight <= 0) return;
  navBand(f, PASS_NAV_OVERLAY, l, av);
  navDrawWindow(e, n);
  be.endPass();
}
function flushMarks(f: ComposerFrame) {
  if (f.fills.n > 0) f.v.backend.rects(f.fills);
  if (f.edges.n > 0) f.v.backend.rects(f.edges);
  f.fills.reset();
  f.edges.reset();
}
/**
 * The hovered row across the header and data columns, and the data area at half strength (a hovered row never reads
 * as selected); a hovered tree node: the rows it stands for, its clade drawn after (`clade`). A hover on a track's
 * data (linear): the `hoverRow` behavior's band at its strengths, the track's box (a split row's, a rows track's
 * row) or its layer's `rowHover`; drawn from the hover as it is, which it never changes.
 */
function rowBand(f: ComposerFrame, l: SectionLayout, s: ViewportState) {
  // The marks' hover (an export's own, empty unless it draws the hover), never the context's directly.
  const h = f.channels?.hovered ?? f.ctx.hover.live.current;
  const v = f.v;
  const r = l.rows;
  if (h.section !== index(f) || h.track < 0) return;
  const own = h.kind === ALIGNED || h.kind === ROW || h.kind === HEADER;
  const onData = l.geometry === LINEAR && (h.region === DATA || h.region === COLUMN);
  const of = own || !onData ? null : f.ctx.ui.rowHover.value;
  const p = of === null ? null : of(h.section);
  if (!own && p?.band !== "entry") return;
  const stacked = !own && v.axis?.stacked === true;
  let y0 = Number.NaN;
  let y1 = Number.NaN;
  if (r.layer >= 0 && h.track === r.track) {
    const top0 = l.scrolledTop + l.entryTop[r.entry] - s.scrollY;
    const sp = f.span;
    sp[0] = sp[1] = SectionLayouts.displayOfLane(r, h.lane);
    const hd = f.renderers[r.layer]?.header;
    // A tree node's rows and clade for a row's or header's hover only (a wedge row's data hover lights none).
    const node = own ? h.treeNode : -1;
    if (node >= 0 && hd?.span) hd.span((f.model as ViewModel).prepared[r.layer].state, node, sp);
    if (sp[0] < 0) return;
    y0 = Math.max(l.scrolledTop, top0 + sp[0] * r.height);
    y1 = Math.min(l.height, top0 + (sp[1] + 1) * r.height);
    if (node >= 0) f.clade = node;
  } else
    for (let e = next(l, s, -1); e >= 0 && Number.isNaN(y0); e = next(l, s, e))
      if (
        l.entryTrack[e] === h.track &&
        (h.split === -1 || l.entryRow[e] === h.split) &&
        boxOf(f, l, e) &&
        (!stacked || (h.block >= 0 ? l.entryBlock[e] === h.block : pointerIn(f)))
      ) {
        if (!own && !layerBand(f, l, e, h)) return;
        y0 = f.box[2];
        y1 = f.box[3];
      }
  if (!(y1 > y0)) return;
  const Y0 = Math.round(y0 * v.dpr) | 0;
  const Y1 = Math.round(y1 * v.dpr) | 0;
  // The strengths as alphas, inline (a double passed to a call is boxed).
  const head = p === null ? 255 : Math.max(0, Math.min(255, Math.round(p.header * 255))) | 0;
  const data = p === null ? 128 : Math.max(0, Math.min(255, Math.round(p.data * 255))) | 0;
  tint(f, BIT_HOVER, head);
  f.fills.push(0, Y0, f.ext[2], Y1 - Y0);
  tint(f, BIT_HOVER, data);
  f.fills.push(f.ext[2], Y0, f.ext[3] - f.ext[2], Y1 - Y0);
}
/**
 * Layer k's own band on entry e (its renderer's `rowHover`), within the entry's visible box (`box[2]`, `box[3]`);
 * false: none. Without the member, the box as it is.
 */
function layerBand(f: ComposerFrame, l: SectionLayout, e: number, h: HoverState): boolean {
  const k = h.layer;
  const r = k >= 0 ? f.renderers[k] : undefined;
  if (!r?.rowHover) return true;
  if (f.failed.has(k)) return false;
  const v = f.v;
  const b = f.box;
  const out = f.bandRows;
  v.top = b[0];
  v.height = b[4];
  place(f, l, e);
  prime(f, k, l, e);
  try {
    const m = f.model as ViewModel;
    if (!r.rowHover(v, stateOf(f, k, l, e), m.resolved.layers.params[k], h, out)) return false;
    b[2] = Math.max(b[2], out[0] / v.dpr);
    b[3] = Math.min(b[3], out[1] / v.dpr);
    return true;
  } catch (err) {
    fail(f, k, err);
    return false;
  } finally {
    release(f);
  }
}
/** Stacked, off a block (past its end): whether the pointer is over the entry in `box` (its block's lane). */
function pointerIn(f: ComposerFrame): boolean {
  const p = f.ctx.pointer.current;
  const b = f.box;
  return p.section === index(f) && p.y >= b[0] && p.y < b[0] + b[1];
}
/** A hovered tree node's clade over the tree (its `header.hover`), clipped to the tree's box over the scrolled rows. */
function cladeHover(f: ComposerFrame, l: SectionLayout, node: number) {
  const r = l.rows;
  const k = r.layer;
  const hd = f.renderers[k]?.header;
  const v = f.v;
  const s = (f.model as ViewModel).prepared[k].state;
  const p = (f.model as ViewModel).resolved.layers.params[k];
  const hw = Math.min(l.headerWidth, l.dataLeft - gutterWidth(f, l));
  const inset = hd?.hover ? Math.min(hd.inset(s, p), hw * 0.5) : 0;
  if (!hd?.hover || inset <= 0 || !(r.height > 0) || l.height <= l.scrolledTop || f.failed.has(k)) return;
  const top0 = l.scrolledTop + l.entryTop[r.entry] - v.scrollY;
  const first = Math.max(0, Math.floor((l.scrolledTop - top0) / r.height));
  const last = Math.min(r.count - 1, Math.floor((l.height - top0) / r.height));
  const T = Math.round(l.scrolledTop * v.dpr) | 0;
  if (first > last) return;
  v.backend.pushClip(0, T, Math.round(inset * v.dpr) | 0, (Math.round(l.height * v.dpr) - T) | 0);
  prime(f, k, l, r.entry);
  v.top = top0;
  v.height = r.height;
  v.left = 0;
  try {
    hd.hover(v, s, p, node, first, last);
  } catch (err) {
    fail(f, k, err);
  } finally {
    release(f);
    v.left = l.dataLeft;
    v.backend.popClip();
  }
}
/**
 * Streamed layers' regions in view (linear): skeletons, error bands and fade-ins over the main pass
 * (seqspec:sources#levels-and-the-zoom-gate; draw/progress.ts), per entry under its clip; `animating` when one
 * shimmers or fades.
 */
function progress(f: ComposerFrame, l: SectionLayout, s: ViewportState) {
  const rv = (f.model as ViewModel).resolved;
  const v = f.v;
  f.animating = false;
  if (l.geometry !== LINEAR) return;
  const reduced = f.ctx.config.get(REDUCED_MOTION) === true;
  f.now[0] = f.ctx.scheduler.now();
  const fps = f.ctx.config.get(SHIMMER_FPS);
  progressBegin(v, f.now, f.live && !reduced, fps > 0);
  for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) {
    const t = l.entryTrack[e];
    if (l.entryKind[e] !== TRACK || t < 0 || !boxOf(f, l, e)) continue;
    let opened = false;
    for (let k = l.entryLayer[e]; k < rv.layers.count && rv.layers.track[k] === t; k++) {
      const st = f.ctx.streams.get(f.refs[k]);
      const rows = st !== null && st.source.systems === null;
      if (st === null || (rows && (k !== l.rows.layer || st.source.remote === null || !banded(f, l, e)))) continue;
      if (!opened && !open(f, l.dataLeft, l.dataRight, 0)) break;
      opened = true;
      v.top = f.box[0];
      v.height = f.box[4];
      place(f, l, e);
      if (rows) progressRows(v, st, l.rows);
      else progress1D(v, st);
    }
    if (!opened) continue;
    progressFlush(v);
    v.backend.popClip();
  }
  const moving = progressEnd();
  f.animating = moving !== 0;
  f.animStep[0] = (moving & FADE) !== 0 ? 0 : 1000 / Math.max(1, fps);
}
/** Each layer's `overlay` member per entry in view (linear, wrapped), clipped to the data area and the entry's band. */
function overlays(f: ComposerFrame, l: SectionLayout, s: ViewportState) {
  const m = f.model as ViewModel;
  const rv = m.resolved;
  const v = f.v;
  for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) {
    const t = l.entryTrack[e];
    if (l.entryKind[e] !== TRACK || banded(f, l, e) || !boxOf(f, l, e)) continue;
    let opened = false;
    const one = SectionLayouts.splitRow(l, e) >= 0;
    for (
      let k = l.entryLayer[e];
      k < rv.layers.count && rv.layers.track[k] === t && !(one && k > l.entryLayer[e]);
      k++
    ) {
      const r = f.renderers[k];
      if (!r?.overlay || f.failed.has(k) || (l.entryRow[e] === -2 && k === rv.tracks.split[t])) continue;
      if (!opened && !open(f, l.dataLeft, l.dataRight, 0)) break;
      opened = true;
      v.top = f.box[0];
      v.height = f.box[4];
      place(f, l, e);
      prime(f, k, l, e);
      try {
        r.overlay(v, stateOf(f, k, l, e), rv.layers.params[k]);
      } catch (err) {
        fail(f, k, err);
      } finally {
        release(f);
      }
    }
    if (opened) v.backend.popClip();
  }
}
/** The keyboard cursor's ring at `cursorRect`, just outside its cell or header, in `--seq-focus-ring`. */
function cursorRing(f: ComposerFrame) {
  const c = f.ctx.ui.keyCursor.value;
  const g = f.cursor;
  const v = f.v;
  if (!c || c.section !== f.section || !cursorRect(f.ctx, c, g, true)) return;
  const w = Math.max(1, Math.round(1.5 * v.dpr)) | 0;
  const x0 = (Math.round(g[0] * v.dpr) - w) | 0;
  const y0 = (Math.round(g[1] * v.dpr) - w) | 0;
  const e = f.edges;
  e.reset();
  e.begin(TOKEN_FOCUS_RING, 255, w);
  e.push(x0, y0, (Math.round((g[0] + g[2]) * v.dpr) + w - x0) | 0, (Math.round((g[1] + g[3]) * v.dpr) + w - y0) | 0);
  v.backend.rects(e);
  e.reset();
}
/**
 * On rings, the window of the linear section linked by `follow` (ring-window.ts): a faint sector over the rings with
 * its radial sides, its arc along the rim (an outermost ruler's baseline) and a grip at each end, through the ring map;
 * across a cyclic origin in two parts.
 */
function ringWindow(f: ComposerFrame, l: SectionLayout) {
  const g = f.rw;
  if (!ringWindowFrame(f.ctx, index(f), g)) return;
  const v = f.v;
  const d = v.dpr;
  const lw = Math.max(1, Math.round(d)) | 0;
  ringBounds(l, f.rb);
  const rim = Math.round(f.rb[0] * d) | 0;
  const y1 = Math.round(f.rb[1] * d) | 0;
  const span = g[RW_HI] - g[RW_LO];
  const wraps = g[RW_CYCLIC] === 1 && g[RW_X1] > g[RW_HI];
  const xa = Math.round(g[RW_X0] * d) | 0;
  const xb = Math.round((wraps ? g[RW_HI] : g[RW_X1]) * d) | 0;
  const xc = Math.round(g[RW_LO] * d) | 0;
  const xe = Math.round((wraps ? g[RW_X1] - span : g[RW_X1]) * d) | 0;
  const fills = f.fills;
  fills.reset();
  // The whole ring (the linked section fits): its grips only, at the start.
  if (g[RW_X1] - g[RW_X0] < span - 0.5) {
    fills.begin(WINDOW_FILL);
    fills.push(xa, rim, Math.max(1, xb - xa), y1 - rim);
    if (wraps && xe > xc) fills.push(xc, rim, xe - xc, y1 - rim);
    const t = Math.max(2, Math.round(1.5 * d)) | 0;
    fills.begin(RW_RIM);
    fills.push(xa, rim - (t >> 1), Math.max(1, xb - xa), t);
    if (wraps && xe > xc) fills.push(xc, rim - (t >> 1), xe - xc, t);
    fills.begin(WINDOW_EDGE_INK);
    fills.push(xa, rim, lw, y1 - rim);
    fills.push(xe - lw, rim, lw, y1 - rim);
  }
  // The grips: pills `GRIP_ALONG` × `GRIP_ACROSS` CSS px straddling the rim, on a halo of the ground.
  const lit = gripsLit(litPart(f.ctx.ui.window.current, index(f), RW_KIND), false);
  const ga = Math.max(2, Math.round(3 * d)) | 0;
  const gc = Math.round(10 * d) | 0;
  fills.begin(RW_HALO);
  fills.push(xa - (ga >> 1) - lw, rim - (gc >> 1) - lw, ga + 2 * lw, gc + 2 * lw);
  fills.push(xe - (ga >> 1) - lw, rim - (gc >> 1) - lw, ga + 2 * lw, gc + 2 * lw);
  fills.begin(lit & 1 ? RW_GRIP_LIT : RW_GRIP);
  fills.push(xa - (ga >> 1), rim - (gc >> 1), ga, gc);
  fills.begin(lit & 2 ? RW_GRIP_LIT : RW_GRIP);
  fills.push(xe - (ga >> 1), rim - (gc >> 1), ga, gc);
  v.backend.rects(fills);
  fills.reset();
}
/** The navigator's view: its parts, the window as drawn and hit, the grips lit, hot while the pointer is over it. */
function navOf(
  f: ComposerFrame,
  be: RenderBackend,
  l: SectionLayout,
  av: SectionAxisView,
  s: ViewportState,
): NavigatorEntry | null {
  const i = index(f);
  const e = f.ctx.layout.navigatorOf(i);
  if (!e || !av.axis || l.geometry !== LINEAR) return null;
  const n = f.nview.from(av);
  const ctx = f.ctx;
  n.backend = be;
  n.text = ctx.text;
  n.parts = e.parts;
  windowFrame(av, l, NAVIGATOR, f.g);
  n.x0 = f.g[8];
  n.x1 = f.g[9];
  const p = ctx.pointer.current;
  const d = ctx.ui.drag.value;
  // Nowhere to go (the whole axis in view at the zoom limit): quiet, no hover, no lit grip; read per frame.
  n.idle = windowIdle(s, av.axis);
  n.hot = !n.idle && (railGrown(d, l, p, i) || (d === null && p.inside && p.section === i && p.y < l.navHeight));
  n.lit = n.idle ? 0 : gripsLit(litPart(ctx.ui.window.current, i, NAVIGATOR), false);
  n.marks = f.marks;
  if (f.namesAt !== av.axis) names(f, av);
  return e;
}
/** Edge-rate: per block its label, else its key. */
function names(f: ComposerFrame, av: SectionAxisView) {
  const m = f.model as ViewModel;
  const ax = av.axis as Axis;
  const b0 = m.resolved.sections.blockStart[index(f)];
  f.namesAt = ax;
  f.nview.names = ax.names.map((k, b) => name(f, m.resolved.blocks.label[b0 + b]) || k);
}
/** The navigator's row as a pass of its own. */
function navBand(f: ComposerFrame, pass: number, l: SectionLayout, av: SectionAxisView) {
  const n = f.nview;
  const W = Math.ceil(l.width * av.dpr) | 0;
  const H = Math.ceil(l.navHeight * av.dpr) | 0;
  n.backend.beginPass(pass, W, H, av.dpr);
  n.backend.clear(0, 0, W, H);
  n.top = n.gripTop = 0;
  n.height = n.gripBottom = l.navHeight;
}
