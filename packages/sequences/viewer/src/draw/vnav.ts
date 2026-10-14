/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vertical overview and the scrollbar (seqspec:viewer#the-render-backend): the composer's other two passes, the
// surface's `vnav` and `vnavOverlay`, drawn into any `RenderBackend`. The overview is the rows layer's `overview` of
// every row in miniature, as an image of the box (`KeyedImage`) built when the layout, the rows layer's regions, the
// fonts, the model or the skin change, behind the window; the window itself, the scrollbar's thumb, both pairs of grips
// and the rows the channels light go on the overlay every frame. Both passes lie in a box of their own, the rows' band
// at the section's right edge, in device px with its origin there (`box`; `renderSvg` offsets them by it). Per frame
// nothing is allocated: numbers in typed arrays and fields, no double passed to a helper (RN1's draw pattern).

import { Lanes } from '@molstar/seq-core/format';
import type { SectionAxisView } from '../axis.js';
import type { ChannelMarks } from '../channels.js';
import type { SeqViewerContext } from '../context.js';
import type { Skin, ViewModel } from '../graph.js';
import { SCROLLBAR } from '../hit.js';
import { hasScrollbar, type RenderProps, type SectionLayout, SectionLayouts } from '../layout.js';
import { MarkBit, type SectionMarks } from '../marks.js';
import type { LayerRenderer } from '../registry.js';
import { gripsLit, litPart, WINDOW_FRAME, WindowKind, windowFrame } from '../windows.js';
import { Pass, RectBatch, type RenderBackend, SkinToken } from './backend.js';
import { Dirty } from './composer.js';
import { MapView } from './miniature.js';
import { drawGrips, gripAt, KeyedImage, WINDOW_WASH, windowBox, windowOverMap } from './nav.js';

// Enum-like members read once: a module member read per frame allocates under the test runner's module wrapper
// (lib/seqstar/AGENTS.md, Tests; VL-15).
const P_VNAV = Pass.vnav;
const P_VNAV_OVERLAY = Pass.vnavOverlay;
const T_SURFACE_ALT = SkinToken.surfaceAlt;
const T_INK_QUIET = SkinToken.inkQuiet;
const T_INK_FAINT = SkinToken.inkFaint;
const T_SELECTION_BORDER = SkinToken.selectionBorder;
const T_HIGHLIGHT = SkinToken.highlight;
const W_OVERVIEW = WindowKind.overview;
const W_SCROLLBAR = WindowKind.scrollbar;
const M_ROW = MarkBit.row;
const M_SELECTION = MarkBit.selection;
const M_HIGHLIGHT = MarkBit.highlight;
const D_NAV = Dirty.nav;
const D_OVERLAY = Dirty.overlay;
const rects = new RectBatch(16);
const marks = new RectBatch(16);
/** The image's source and destination (device px). */
const img = new Float64Array(8);

export class VerticalNav {
  /** The box both passes' canvases cover in the section, device px: left, top, width, height; width 0: none to draw. */
  readonly box = new Int32Array(4);
  private readonly ctx: SeqViewerContext;
  private readonly section: string;
  private readonly channels: ChannelMarks;
  private readonly map = new MapView();
  private zeros = new Float64Array(0);
  private readonly key = new KeyedImage(8);
  /** The overview and the scrollbar as windows.ts frames them. */
  private readonly f = new Float64Array(WINDOW_FRAME);
  private readonly s = new Float64Array(WINDOW_FRAME);
  /** The pass as drawn: layout version, box width, height, DPR, image. */
  private readonly at = new Float64Array(5).fill(Number.NaN);
  /** The overview's origin in the box and its size (device px): x, width. */
  private readonly ovl = new Int32Array(2);
  private model: ViewModel | null = null;
  private skin: Skin | null = null;
  private failed = false;
  /** Bumped when the model or the skin changes (a key's number for both). */
  private gen = 0;
  private entry: { state: unknown; r: LayerRenderer; c: unknown } | null = null;
  /** The rows layer's data row → display row (-1: not shown), at layout version `displayAt`. */
  private display = new Int32Array(0);
  private displayAt = -1;

  constructor(ctx: SeqViewerContext, section: string, channels: ChannelMarks) {
    this.ctx = ctx;
    this.section = section;
    this.channels = channels;
  }

  /** The box for layout `l` at `dpr` into `box`: the scrollbar and the overview from the section's right edge, below the pinned rows. */
  place(l: SectionLayout, p: RenderProps, dpr: number): void {
    const b = this.box;
    const bar = hasScrollbar(p) && l.vscroll ? SCROLLBAR : 0;
    if (bar + l.vnavWidth <= 0 || l.height <= l.scrolledTop) {
      b[0] = b[1] = b[2] = b[3] = 0;
      return;
    }
    const W = Math.ceil(l.width * dpr) | 0;
    const H = Math.ceil(l.height * dpr) | 0;
    b[0] = Math.round((l.width - l.vnavWidth - bar) * dpr) | 0;
    b[1] = Math.round(l.scrolledTop * dpr) | 0;
    b[2] = Math.max(0, W - b[0]);
    b[3] = Math.max(0, H - b[1]);
    this.ovl[0] = (Math.round((l.width - l.vnavWidth) * dpr) | 0) - b[0];
    this.ovl[1] = l.vnavWidth > 0 ? b[2] - this.ovl[0] : 0;
  }

  /** The next overview pass draws whole (a canvas was resized). */
  redraw(): void {
    this.at[0] = Number.NaN;
  }

  /** Draws what `flags` (`Dirty`) ask: the overview on `nav`, its window and the scrollbar on `overlay`. */
  frame(be: RenderBackend, flags: number): void {
    const ctx = this.ctx;
    const live = ctx.layout.get(this.section);
    const av = ctx.layout.axis(this.section);
    const model = ctx.graph.assembled.value?.view ?? null;
    if (model !== this.model) this.setModel(model);
    if (ctx.skin.value !== this.skin) this.setSkin(ctx.skin.value);
    if (!live || !av?.axis || !model) return;
    const l = live.current;
    const p = ctx.layout.props(this.section);
    this.place(l, p, av.dpr);
    if (this.box[2] <= 0 || this.box[3] <= 0) return;
    if (flags & D_NAV) this.nav(be, l, av, live.version.value);
    if (flags & D_OVERLAY) this.overlay(be, l, av, p);
  }

  dispose(): void {
    this.key.dispose();
    this.drop();
    this.model = this.skin = null;
  }

  private setModel(m: ViewModel | null) {
    this.model = m;
    this.failed = false;
    this.gen++;
    this.drop();
    this.at[0] = Number.NaN;
  }
  private setSkin(skin: Skin) {
    this.skin = skin;
    this.gen++;
    this.drop();
    this.at[0] = Number.NaN;
  }
  private drop() {
    this.entry?.r.cache?.dispose?.(this.entry.c);
    this.entry = null;
    this.key.dispose();
  }

  // ---- The overview ------------------------------------------------------------------------------------------------

  private nav(be: RenderBackend, l: SectionLayout, av: SectionAxisView, version: number) {
    const W = this.box[2];
    const H = this.box[3];
    const id = l.vnavWidth > 0 ? this.image(be, l, av, version) : -1;
    const a = this.at;
    if (a[0] === version && a[1] === W && a[2] === H && a[3] === av.dpr && a[4] === id) return;
    a[0] = version;
    a[1] = W;
    a[2] = H;
    a[3] = av.dpr;
    a[4] = id;
    be.beginPass(P_VNAV, W, H, av.dpr);
    be.clear(0, 0, W, H);
    if (l.vnavWidth > 0) {
      const x = this.ovl[0];
      const w = this.ovl[1];
      rects.reset();
      rects.begin(T_SURFACE_ALT);
      rects.push(x, 0, w, H);
      be.rects(rects);
      if (id >= 0) {
        img[0] = img[1] = img[5] = 0;
        img[2] = img[6] = w;
        img[3] = img[7] = H;
        img[4] = x;
        be.image(id, img, false);
      }
    }
    be.endPass();
  }

  /** The overview's image under its key (the box's size, the layout, its layer's regions, the fonts, the model, the skin); -1: none. */
  private image(be: RenderBackend, l: SectionLayout, av: SectionAxisView, version: number): number {
    const ctx = this.ctx;
    const m = this.model as ViewModel;
    const r = l.rows;
    const w = this.ovl[1];
    const h = this.box[3];
    const rend = r.layer >= 0 ? ctx.registry.get('renderers', m.prepared[r.layer].renderer) : undefined;
    if (!rend?.overview || this.failed || w < 2 || h < 2 || !av.axis) return -1;
    const images = be.images;
    const k = this.key;
    k.next[0] = w;
    k.next[1] = h;
    k.next[2] = av.dpr;
    k.next[3] = version;
    // The rows layer's own regions (another layer's arrival leaves the image); a layer not streamed: -1.
    const stream = ctx.streams.get(m.set.strings.get(m.resolved.layers.ref[r.layer]));
    k.next[4] = stream === null ? -1 : stream.version;
    k.next[5] = ctx.text.version.value;
    k.next[6] = this.gen;
    if (!k.stale(images)) return k.id;
    const into = images.beginLayer(w, h);
    const ok = this.render(into, rend, l, av, w);
    const id = images.endLayer();
    if (ok) k.keep(images, id);
    else {
      images.dispose(id);
      k.keep(images, -1);
    }
    return k.id;
  }

  /** The rows layer's `overview` into a `w` device px wide image, as a map: one column of cells per alignment column. */
  private render(into: RenderBackend, rend: LayerRenderer, l: SectionLayout, av: SectionAxisView, w: number) {
    const ctx = this.ctx;
    const m = this.model as ViewModel;
    const ax = av.axis as NonNullable<SectionAxisView['axis']>;
    const k = l.rows.layer;
    const p = m.prepared[k];
    if (this.zeros.length !== ax.count) this.zeros = new Float64Array(ax.count);
    const v = this.map.of(av, this.zeros, 0, w / av.dpr / Math.max(1, ax.cells), w / av.dpr);
    v.backend = into;
    v.skin = ctx.skin.value;
    v.props = ctx.layout.props(this.section);
    v.text = ctx.text;
    v.top = 0;
    v.height = this.box[3] / av.dpr;
    v.lanes.fill(0);
    v.layer = k;
    v.setThemes(p.themes);
    v.marks = this.channels.marks(ctx.layout.sectionIndex(this.section)) ?? v.marks;
    v.stream = ctx.streams.get(m.set.strings.get(m.resolved.layers.ref[k]));
    v.cache = rend.cache ? this.cacheOf(rend, p.state) : null;
    try {
      return rend.overview?.(v, p.state, m.resolved.layers.params[k]) === true;
    } catch (err) {
      this.failed = true;
      ctx.report(
        'seqv.renderer.error',
        `layer ${k} threw while drawn in the vertical overview and is omitted: ${String(err)}`,
      );
      return false;
    } finally {
      v.cache = v.stream = null;
    }
  }
  private cacheOf(r: LayerRenderer, state: unknown): unknown {
    if (this.entry?.state !== state) {
      this.entry?.r.cache?.dispose?.(this.entry.c);
      this.entry = { state, r, c: r.cache?.create(state, this.ctx.skin.value) };
    }
    return (this.entry as { c: unknown }).c;
  }

  // ---- The window, the scrollbar and the marks ---------------------------------------------------------------------

  private overlay(be: RenderBackend, l: SectionLayout, av: SectionAxisView, p: RenderProps) {
    const ctx = this.ctx;
    const d = av.dpr;
    const ox = this.box[0];
    const oy = this.box[1];
    const W = this.box[2];
    const H = this.box[3];
    const lw = Math.max(1, Math.round(d)) | 0;
    const i = ctx.layout.sectionIndex(this.section);
    const zoom = p.verticalZoom;
    be.beginPass(P_VNAV_OVERLAY, W, H, d);
    be.clear(0, 0, W, H);
    if (l.vnavWidth > 0) {
      const f = this.f;
      const x = this.ovl[0];
      windowFrame(av, l, W_OVERVIEW, f);
      // The map washed outside the window, the window at full color in a two-tone outline (reads on any map).
      const w = windowBox;
      w[0] = (Math.round(f[8] * d) | 0) - ox;
      w[1] = (Math.round(f[10] * d) | 0) - oy;
      w[2] = (Math.round(f[9] * d) | 0) - ox;
      w[3] = (Math.round(f[11] * d) | 0) - oy;
      w[4] = x;
      w[5] = 0;
      w[6] = W;
      w[7] = H;
      windowOverMap(be, lw, T_SURFACE_ALT, WINDOW_WASH);
      this.rowMarks(be, l, x);
      // Grips where its ends resize: left and right (columns) always, top and bottom (rows) with `verticalZoom`.
      const lit = litPart(ctx.ui.window.current, i, W_OVERVIEW);
      const a = gripAt;
      a[0] = w[0];
      a[1] = w[2];
      a[2] = x;
      a[3] = W;
      a[4] = (w[1] + w[3]) >> 1;
      a[5] = w[1];
      a[6] = w[3];
      a[7] = d;
      drawGrips(be, false, gripsLit(lit, false));
      if (zoom) {
        a[0] = w[1];
        a[1] = w[3];
        a[2] = 0;
        a[3] = H;
        a[4] = (w[0] + w[2]) >> 1;
        a[5] = w[0];
        a[6] = w[2];
        drawGrips(be, true, gripsLit(lit, true));
      }
    }
    if (hasScrollbar(p) && l.vscroll) {
      const s = this.s;
      windowFrame(av, l, W_SCROLLBAR, s);
      // The thumb at its true share (inset 1 px where it is long enough), at least 2 device px so a sliver shows.
      const inset = s[11] - s[10] >= 4 ? 1 : 0;
      const y0 = (Math.round((s[10] + inset) * d) | 0) - oy;
      const th = Math.max(2, (Math.round((s[11] - inset) * d) | 0) - oy - y0);
      const bx = (Math.round((l.width - l.vnavWidth - SCROLLBAR) * d) | 0) - ox;
      const bw = Math.max(1, Math.round(SCROLLBAR * d) | 0);
      const inner = Math.max(1, Math.round((SCROLLBAR - 4) * d) | 0);
      rects.reset();
      rects.begin(ctx.ui.drag.value?.owner === 'scrollbar' ? T_INK_QUIET : T_INK_FAINT);
      rects.push(bx + ((bw - inner) >> 1), Math.min(y0, H - th), inner, th);
      be.rects(rects);
      // Its ends resize rows only with `verticalZoom`: grips there, none otherwise (a grip says what a drag does).
      if (zoom) {
        const g = gripAt;
        g[0] = y0;
        g[1] = y0 + th;
        g[2] = 0;
        g[3] = H;
        g[4] = bx + (bw >> 1);
        g[5] = bx;
        g[6] = bx + bw;
        g[7] = d;
        drawGrips(be, true, gripsLit(litPart(ctx.ui.window.current, i, W_SCROLLBAR), true));
      }
    }
    be.endPass();
  }
  /** The data row → display row map, rebuilt per layout version (edge-rate). */
  private displayRows(l: SectionLayout, version: number): Int32Array {
    if (this.displayAt === version) return this.display;
    this.displayAt = version;
    const r = l.rows;
    const lanes = (this.model as ViewModel).resolved.lanes;
    let n = 0;
    for (let q = 0; q < r.count; q++) {
      const lane = SectionLayouts.laneOf(r, q);
      if (lane >= 0) n = Math.max(n, Lanes.rowOf(lanes, lane) + 1);
    }
    if (this.display.length < n) this.display = new Int32Array(n);
    this.display.fill(-1);
    for (let q = 0; q < r.count; q++) {
      const lane = SectionLayouts.laneOf(r, q);
      if (lane >= 0) this.display[Lanes.rowOf(lanes, lane)] = q;
    }
    return this.display;
  }
  /** The rows the channels light (selection, highlights) as bars across the overview: at least 2 device px tall. */
  private rowMarks(be: RenderBackend, l: SectionLayout, x: number) {
    const m: SectionMarks | null = this.channels.marks(this.ctx.layout.sectionIndex(this.section));
    const r = l.rows;
    if (!m || r.layer < 0 || r.count === 0 || m.rows.count === 0) return;
    const rm = m.rows;
    const map = this.displayRows(l, this.ctx.layout.get(this.section)?.version.value ?? -1);
    const H = this.box[3];
    const w = this.ovl[1];
    marks.reset();
    for (let q = 0; q < rm.count; q++) {
      const bit = rm.bit[q] & ~M_ROW;
      if (bit !== M_SELECTION && bit !== M_HIGHLIGHT) continue;
      marks.begin(bit === M_SELECTION ? T_SELECTION_BORDER : T_HIGHLIGHT);
      for (let row = rm.rowStart[q]; row < rm.rowStop[q]; ) {
        const first = row < map.length ? map[row] : -1;
        let n = 1;
        while (first >= 0 && row + n < rm.rowStop[q] && row + n < map.length && map[row + n] === first + n) n++;
        row += n;
        if (first < 0) continue;
        const y0 = Math.floor((first * H) / r.count) | 0;
        const y1 = Math.max(y0 + 2, Math.floor(((first + n) * H) / r.count) | 0);
        marks.push(x, y0, w, Math.min(y1, H) - y0);
      }
    }
    if (marks.n > 0) be.rects(marks);
  }
}
