/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section composer (seqspec:viewer#the-render-backend): a section's passes drawn into any `RenderBackend`,
// DOM-free, so the surface drives it over Canvas2D and a test or an SVG figure over a recording. It knows no renderer
// by name: each layer goes through its registry entry (the prepared layer's renderer, so a fallback draws as itself):
// `draw`, `tiny` below `minHeight` (else its extents as bars), `rows` for dense or lazy rows, `gutter`, the header's
// `rowModel`, `bands` and `header.draw`. The main pass draws into the backend's back target and is presented whole; a
// pan at the same zoom copies the previous image (`Copy`) and draws the uncovered strips, then the gutters once. A
// layer's stream version redraws only its entries, within the columns of the regions that changed. The overlay draws
// what the channels light (`ChannelMarks`, `SectionMarks`): each block's marks over its lane through one `paintMark`,
// row marks across header, data columns and data under them, the hovered row at half strength over the data, renderers'
// `overlay` members and the keyboard cursor's ring; under them, streamed layers' skeletons, error bands and fade-ins
// (draw/progress.ts; `animating` asks for the next overlay frame). Every draw runs in a try/finally: a throw omits the
// layer, reported once, and the view keeps none of its inputs. Per frame nothing is allocated: numbers in typed arrays
// and fields, no closure in a per-frame function, no double passed to a helper (RN1's draw pattern).
//
// Its parts beside it: composer-frame.ts (`ComposerFrame`, the state, which this class extends, and the geometry
// every pass shares), composer-entries.ts, composer-headers.ts, composer-marks.ts and composer-overlay.ts. The ring
// labels are placed here (edge rate), not in ring-labels.ts, which the hit imports and which imports no layout.

import { AxisGeometry, type SectionAxisView, type ShownIndex } from '../axis.js';
import { ChannelMarks } from '../channels.js';
import type { Skin, ViewModel } from '../graph.js';
import { EntryKind, type SectionLayout } from '../layout.js';
import type { SectionMarks } from '../marks.js';
import { lanesInView } from '../streams.js';
import { Lod, type ViewportState } from '../viewport.js';
import { BackendCaps, Pass, type RenderBackend } from './backend.js';
import { bandCopies, cutRows, entry } from './composer-entries.js';
import {
  ALL,
  acrossOf,
  banded,
  bands,
  breaks,
  ComposerFrame,
  gutterWidth,
  index,
  next,
  Part,
  setClip,
  view,
} from './composer-frame.js';
import { heads } from './composer-headers.js';
import { nav, overlay } from './composer-overlay.js';
import { placedRingLabels, type RingLabelInput } from './ring-labels.js';
import { nextFrameId } from './view.js';

// Enum-like members read once (a read per frame allocated under the test runner's module wrapper).
const CIRCULAR = AxisGeometry.circular;
/** The imported `ALL`, read once: an imported binding read per frame is a module member read under the test runner. */
const PART_ALL = ALL;
const LINEAR = AxisGeometry.linear;
const WRAPPED_GEOMETRY = AxisGeometry.wrapped;
const CAPS_BACK = BackendCaps.Back;
const CAPS_COPY = BackendCaps.Copy;
const TRACK = EntryKind.track;
const LETTERS = Lod.letters;
const PART_CLEAR = Part.clear;
const PART_DATA = Part.data;
const PART_GUTTER = Part.gutter;
const PART_HEADER = Part.header;
const PASS_MAIN = Pass.main;

/** What a frame draws: the tracks (main pass), the overlay, the navigator. */
export const Dirty = { tracks: 1, overlay: 2, nav: 4, all: 7 } as const;

export type { ComposerStats } from './composer-frame.js';

export class SectionComposer extends ComposerFrame {
  /** Draws what `flags` (`Dirty`) ask into `be`, whole; the tracks only where the view, layout, data or skin changed. */
  frame(be: RenderBackend, flags: number): void {
    const ctx = this.ctx;
    const live = ctx.layout.get(this.section);
    const av = ctx.layout.axis(this.section);
    const vp = ctx.viewports.get(this.section);
    const model = ctx.graph.assembled.value?.view ?? null;
    if (model !== this.model) this.setModel(model);
    if (ctx.skin.value !== this.skin) this.setSkin(ctx.skin.value);
    this.stats.frames++;
    this.main.frameId = this.ov.frameId = nextFrameId();
    if (!live || !av?.axis || !vp || !model) {
      // Nothing drawn, nothing moving: a section that left the view asks for no more frames (VL-14).
      this.animating = false;
      this.again = false;
      return;
    }
    const l = live.current;
    this.lay = l;
    // Rings draw the axis unrolled through the polar map (ring.ts).
    if (l.geometry === CIRCULAR) be = this.rings(be, l, av, live.version.value);
    view(this, this.v, be, av);
    this.v.numberingFrom(model.resolved, ctx.layout.sectionIndex(this.section));
    this.v.lanes[0] = this.v.lanes[1] = 0;
    lanesInView(l, this.v, this.v.lanes);
    this.marks = this.marksOf(av);
    this.ext[2] = Math.round(l.dataLeft * av.dpr) | 0;
    this.ext[3] = Math.round(l.dataRight * av.dpr) | 0;
    if (flags & Dirty.tracks) {
      this.v.more = false;
      this.tracks(be, l, vp.live.current, live.version.value);
      // A draw left work for another frame (a dot plot's tiles past its budget): the surface draws the tracks again,
      // those layers' entries alone (`pend`).
      this.again = this.v.more;
    }
    if (flags & Dirty.nav) nav(this, be, l, av, vp.live.current, live.version.value);
    if (flags & Dirty.overlay) overlay(this, be, l, av, vp.live.current, live.version.value);
  }
  /** The ring map at a new layout or DPR only (edge rate: doubles kept off the frame), and the backend through it. */
  private rings(be: RenderBackend, l: SectionLayout, av: SectionAxisView, version: number): RenderBackend {
    if (this.ringAt[0] !== version || this.ringAt[1] !== av.dpr || this.ring.view !== av) {
      this.ringAt[0] = version;
      this.ringAt[1] = av.dpr;
      this.ringBe.map(this.ring.from(l, av, av.dpr));
    }
    this.ringBe.to = be;
    return this.ringBe;
  }
  /** The next tracks frame draws whole, and the navigator again (tests compare copies with it; benchmarks time it). */
  redraw(): void {
    this.drawn[4] = this.navAt[0] = Number.NaN;
  }
  dispose(): void {
    this.setModel(null);
    this.skin = null;
    if (this.own) this.channels?.dispose();
    this.mini.dispose();
    this.channels = this.marks = null;
  }
  /** The section's marks, its changed channels rebuilt; the hover's primary again on a new axis. */
  private marksOf(av: SectionAxisView): SectionMarks | null {
    this.channels ??= new ChannelMarks(this.ctx);
    const m = this.channels.marks(index(this));
    if (m && av.axis && this.primaryAt !== av.axis) m.updatePrimary(av, av.axis.count);
    this.primaryAt = av.axis;
    return m;
  }
  /** Edge-rate: the layers' entries and streams, failures and caches of states no longer prepared. */
  private setModel(model: ViewModel | null) {
    this.model = model;
    this.failed.clear();
    this.undrawn.clear();
    this.unringed.clear();
    this.redraw();
    const keep = new Set(model?.prepared.map((p) => p.state));
    for (const [state, e] of this.caches) {
      if (keep.has(state)) continue;
      e.r.cache?.dispose?.(e.c);
      this.caches.delete(state);
    }
    const p = model?.prepared ?? [];
    this.renderers = p.map((x) => this.ctx.registry.get('renderers', x.renderer));
    this.refs = p.map((_, k) => model?.set.strings.get(model.resolved.layers.ref[k]) ?? '');
    this.streams = p.map(() => null);
    this.streamAt = new Float64Array(p.length).fill(-1);
    this.fresh = new Uint8Array(p.length);
    this.pend = new Uint8Array(p.length);
    this.cols = new Int32Array(2 * p.length);
    acrossOf(this, model);
  }
  private setSkin(skin: Skin) {
    this.skin = skin;
    for (const e of this.caches.values()) e.r.cache?.dispose?.(e.c);
    this.caches.clear();
    this.redraw();
  }
  private tracks(be: RenderBackend, l: SectionLayout, s: ViewportState, version: number) {
    const d = this.drawn;
    const dpr = s.dpr;
    const W = Math.ceil(l.width * dpr) | 0;
    const H = Math.ceil(l.height * dpr) | 0;
    const dx = Math.round((d[0] - s.scrollX) * dpr) | 0;
    const dy = Math.round((d[1] - s.scrollY) * dpr) | 0;
    // A highlight of style `band` changed (rare): the tracks whole, as on a skin change.
    const m = this.marks;
    const hl = this.channels?.rebuilds.highlight ?? -1;
    if (hl !== this.bandAt && (this.bandsDrawn || (m !== null && m.band.count > 0))) d[4] = Number.NaN;
    this.bandAt = hl;
    this.bandsDrawn = m !== null && m.band.count > 0;
    const same = d[2] === s.cell && d[3] === s.rowScale && d[4] === version && d[5] === dpr && d[6] === W && d[7] === H;
    const arrived = this.arrivals(l, s);
    if (same && dx === 0 && dy === 0 && !arrived) return;
    if (this.version !== version) this.lineText.length = 0;
    this.version = version;
    const props = this.v.props;
    const banded = l.rows.layer >= 0 && (l.rows.lazy || l.rows.height < props.denseRowHeight);
    // A band copies unless it draws a tree; its regions pending or not (a strip equals a full draw), the changes since
    // the last frame redrawn after the copy.
    const copies = !banded || bandCopies(this, l);
    const fits = Math.abs(dx) < (l.dataRight - l.dataLeft) * dpr && Math.abs(dy) < (l.height - l.scrolledTop) * dpr;
    be.beginPass(PASS_MAIN, W, H, dpr);
    if (l.geometry === CIRCULAR) this.full(be, l, s, W, H);
    else if (same && dx === 0 && dy === 0) this.redraws(be, l, s, false);
    else if (same && fits && props.blit && copies && (be.caps & CAPS_COPY) !== 0) this.copy(be, l, s, dx, dy);
    else this.full(be, l, s, W, H);
    if (be.caps & CAPS_BACK) be.present(0, 0, W, H);
    be.endPass();
    d[0] = s.scrollX;
    d[1] = s.scrollY;
    d[2] = s.cell;
    d[3] = s.rowScale;
    d[4] = version;
    d[5] = dpr;
    d[6] = W;
    d[7] = H;
  }
  /**
   * Layers whose stream version moved (a region arrived, failed or left): `fresh`, with the columns the changes cover
   * (`cols`; -1: all); whether one is in view.
   */
  private arrivals(l: SectionLayout, s: ViewportState): boolean {
    const rv = (this.model as ViewModel).resolved;
    const c = this.cols;
    for (let k = 0; k < this.refs.length; k++) {
      const st = this.ctx.streams.get(this.refs[k]);
      this.streams[k] = st;
      const moved = st !== null && st.version !== this.streamAt[k];
      const pend = this.pend[k] === 1;
      this.pend[k] = 0;
      this.fresh[k] = moved || pend ? 1 : 0;
      if (pend) c[2 * k] = c[2 * k + 1] = -1;
      if (!st || !moved) continue;
      const n = st.changed ? st.changed(this.streamAt[k], this.changes) : -1;
      c[2 * k] = n < 0 || pend ? -1 : 0x7fffffff;
      c[2 * k + 1] = -1;
      for (let q = 0; q < n && !pend; q++) {
        c[2 * k] = Math.min(c[2 * k], this.changes[2 * q]);
        c[2 * k + 1] = Math.max(c[2 * k + 1], this.changes[2 * q + 1]);
      }
      this.streamAt[k] = st.version;
    }
    for (let e = next(l, s, -1); e >= 0; e = next(l, s, e))
      for (let k = l.entryLayer[e]; k >= 0 && k < rv.layers.count && rv.layers.track[k] === l.entryTrack[e]; k++)
        if (this.fresh[k] === 1) return true;
    return false;
  }
  private full(be: RenderBackend, l: SectionLayout, s: ViewportState, W: number, H: number) {
    be.clear(0, 0, W, H);
    setClip(this, 0, W, 0, H);
    bands(this, l);
    this.dense = false;
    for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) entry(this, l, e, PART_ALL);
    breaks(this, l);
    heads(this, l);
    if (l.geometry === CIRCULAR) this.labelRing(l);
    this.stats.fullDraws++;
  }
  /** Rings: the labeled items' labels outside the rings, placed again at a new layout, DPR or model (edge rate). */
  private labelRing(l: SectionLayout) {
    const dpr = this.v.dpr;
    const m = this.model as ViewModel;
    if (this.labelsAt[0] !== this.version || this.labelsAt[1] !== dpr || this.labelsOf !== m) {
      this.labelsAt[0] = this.version;
      this.labelsAt[1] = dpr;
      this.labelsOf = m;
      this.placeRingLabels(l, m);
      this.ringLabels.at = this.version;
      placedRingLabels.set(l, this.ringLabels);
    }
    this.ringLabels.draw(this.ringBe.to, this.ctx.text, dpr);
  }
  /**
   * The candidates from each layer whose renderer lists them (`ringLabels`), at their items' middle angles and from
   * their rings' outer edges; an item's runs (across a cyclic origin) once, at the longer. Their radius: past an
   * outermost ruler's numbers.
   */
  private placeRingLabels(l: SectionLayout, m: ViewModel) {
    const av = this.ctx.layout.axis(this.section);
    const ring = this.ring;
    const rv = m.resolved;
    const dpr = this.v.dpr;
    const items: (RingLabelInput & { extent: number })[] = [];
    const byItem = new Map<string, number>();
    for (let e = 0; av && e < l.entryCount; e++) {
      if (l.entryKind[e] !== TRACK) continue;
      const t = l.entryTrack[e];
      const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e];
      const r0 = ring.r - (top - ring.top) * ring.scale;
      for (let k = l.entryLayer[e]; k >= 0 && k < rv.layers.count && rv.layers.track[k] === t; k++) {
        this.renderers[k]?.ringLabels?.(m.prepared[k].state, rv.layers.params[k], {
          add: (b, start, stop, text, a, row) => {
            // The middle position's cell (an item's run may cross a rotated block's origin on screen), a boundary's edge.
            const q = av.shownOf(b, stop > start ? (start + stop - 1) >> 1 : start);
            if (q < 0) return;
            const x0 = av.x(b, q as ShownIndex);
            const x1 = stop > start ? av.x(b, (q + 1) as ShownIndex) : x0;
            const extent = stop - start;
            const id = `${k}:${a}:${m.set.annotations[a]?.item[row] ?? row}`;
            const was = byItem.get(id);
            if (was !== undefined && items[was].extent >= extent) return;
            const ux = (x0 + x1) / 2;
            const at = { layer: k, annotation: a, row, ux, uy: top + 0.5 };
            const x = { angle: ring.angleAt(ux * dpr), r0, text, priority: 0, extent, ...at };
            const priority = k * 2 ** 24 + 2 ** 24 - 1 - Math.min(extent, 2 ** 24 - 1);
            if (was === undefined) byItem.set(id, items.push({ ...x, priority }) - 1);
            else items[was] = { ...x, priority };
          },
        });
      }
    }
    const k0 = l.pinnedCount;
    const ruler = k0 < l.entryCount && l.ringTop > l.scrolledTop + l.entryTop[k0] + 0.5;
    const t = this.ctx.text;
    const numbers = ruler ? t.numberWidth(t.font('ui', 11), Math.max(1, av?.axis?.cells ?? 1)) + 10 : 0;
    this.ringLabels.place(
      items,
      {
        cx: l.ringCx,
        cy: l.ringCy,
        edge: l.ringR + 2,
        r: l.ringR + numbers + 10,
        left: 0,
        right: l.width,
        top: l.scrolledTop,
        bottom: l.height,
        line: 13,
      },
      t,
    );
    // Labels a side cannot hold are hidden, and said once per section (a wider section, or fewer labeled layers).
    if (this.ringLabels.hidden > 0 && !this.labelsHidden) {
      this.labelsHidden = true;
      this.ctx.report(
        'seqv.render.labels',
        `${this.ringLabels.hidden} of ${items.length} labels on the rings of section "${this.section}" do not fit beside them: hidden`,
      );
    }
  }
  /**
   * A pan or scroll at the same zoom: the data, pinned and header regions moved by one copy each (one snapshot), the
   * gutters cleared; then the uncovered strip and band, the entries a region's edge cuts and unstable layers whole, and
   * the gutters once.
   */
  private copy(be: RenderBackend, l: SectionLayout, s: ViewportState, dx: number, dy: number) {
    const v = this.v;
    const dpr = v.dpr;
    const W = this.drawn[6] | 0;
    const H = Math.ceil(l.height * dpr) | 0;
    const L = Math.round(l.dataLeft * dpr) | 0;
    const R = Math.round(l.dataRight * dpr) | 0;
    const P = Math.round(l.navHeight * dpr) | 0;
    const T = Math.round(l.scrolledTop * dpr) | 0;
    const G = Math.round((l.dataLeft - gutterWidth(this, l)) * dpr) | 0;
    be.copy(dx, 0, L, P, R - L, T - P);
    be.copy(dx, dy, L, T, R - L, H - T);
    be.copy(0, dy, 0, T, G, H - T);
    be.clear(G, P, L - G, H - P);
    be.clear(R, P, W - R, H - P);
    this.stats.copies++;
    if (dx !== 0) {
      // The uncovered strip, widened to whole cells at letter zoom (no letter cut where strip and copy meet).
      let s0 = dx > 0 ? L : R + dx;
      let s1 = dx > 0 ? L + dx : R;
      if (v.lod === LETTERS) {
        if (dx > 0) s1 = Math.min(R, this.cellEdge(s1, 1));
        else s0 = Math.max(L, this.cellEdge(s0, -1));
      }
      be.clear(s0, P, s1 - s0, H - P);
      setClip(this, s0, s1, P, H);
      bands(this, l);
      const lo = v.visLo;
      const hi = v.visHi;
      this.strip(s0, s1);
      this.dense = false;
      for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) entry(this, l, e, PART_DATA);
      breaks(this, l);
      v.visLo = lo;
      v.visHi = hi;
    }
    if (dy !== 0) {
      // The band uncovered by the scroll, widened to whole entries (no row or label cut where band and copy meet).
      let y0 = dy > 0 ? T : H + dy;
      let y1 = dy > 0 ? T + dy : H;
      for (let e = next(l, s, l.pinnedCount - 1); e >= 0; e = next(l, s, e)) {
        let top = (l.scrolledTop + l.entryTop[e] - s.scrollY) * dpr;
        let bottom = top + l.entryHeight[e] * dpr;
        if (bottom <= y0 || top >= y1) continue;
        if (banded(this, l, e)) {
          // A band widens to its whole rows, not the entry (a lazy band is one entry).
          const h = l.rows.height * dpr;
          bottom = Math.min(bottom, top + Math.ceil((y1 - top) / h) * h);
          top += Math.max(0, Math.floor((y0 - top) / h)) * h;
        }
        y0 = Math.max(T, Math.min(y0, Math.floor(top)));
        y1 = Math.min(H, Math.max(y1, Math.ceil(bottom)));
      }
      be.clear(0, y0, W, y1 - y0);
      setClip(this, 0, W, y0, y1);
      bands(this, l);
      this.dense = false;
      for (let e = next(l, s, l.pinnedCount - 1); e >= 0; e = next(l, s, e)) entry(this, l, e, PART_HEADER | PART_DATA);
      breaks(this, l);
    }
    this.redraws(be, l, s, true);
  }
  /**
   * After a move, the entries a scrolled region's edge cuts (Skia drops a clipped glyph's edge rows: drawn under the
   * clip a full draw has) and the unstable layers' entries whole, then the gutters once; and with or without a move,
   * the entries of layers whose stream changed.
   */
  private redraws(be: RenderBackend, l: SectionLayout, s: ViewportState, moved: boolean) {
    const dpr = this.v.dpr;
    const W = this.drawn[6] | 0;
    const H = Math.ceil(l.height * dpr) | 0;
    for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) {
      const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - s.scrollY;
      const bottom = top + l.entryHeight[e];
      const cut = (top < l.scrolledTop && bottom > l.scrolledTop) || (top < l.height && bottom > l.height);
      // Its visible band: what is redrawn, its break lines included (never drawn twice over).
      const y0 = Math.round(Math.max(top, e < l.pinnedCount ? l.navHeight : l.scrolledTop) * dpr) | 0;
      const y1 = Math.round(Math.min(bottom, e < l.pinnedCount ? l.scrolledTop : l.height) * dpr) | 0;
      this.dense = false;
      if (y1 <= y0) continue;
      // A band whose rows layer changed in every column (its row labels arrived, or the log overran): whole, its
      // labels and data columns with it.
      const k = l.rows.layer;
      if (k >= 0 && this.fresh[k] === 1 && this.cols[2 * k] < 0 && banded(this, l, e)) {
        be.clear(0, y0, W, y1 - y0);
        setClip(this, 0, W, y0, y1);
        bands(this, l);
        entry(this, l, e, PART_HEADER | PART_DATA);
        breaks(this, l);
        continue;
      }
      const rows = moved && e >= l.pinnedCount && cut && l.rows.lazy && banded(this, l, e);
      if (rows) cutRows(this, be, l, s, W, y0, y1);
      if (moved && e >= l.pinnedCount && cut && !rows) {
        // A chip after another on its row (seqspec:viewer#the-text-view): the row's first entry of the track cleared
        // it, and a clear of its own would wipe the chips drawn before it.
        if (!chipAfter(l, e)) be.clear(0, y0, W, y1 - y0);
        setClip(this, 0, W, y0, y1);
        bands(this, l);
        entry(this, l, e, PART_HEADER | PART_DATA);
      } else if (this.redrawn(l, e, moved)) {
        // Within the columns its streams' changes cover (`span`), else the whole data box.
        setClip(this, this.span[0], this.span[1], y0, y1);
        entry(this, l, e, PART_DATA | PART_CLEAR);
      } else if (!rows) continue;
      breaks(this, l);
    }
    setClip(this, 0, W, 0, H);
    for (let e = next(l, s, -1); moved && e >= 0; e = next(l, s, e)) entry(this, l, e, PART_GUTTER);
  }
  /**
   * An entry drawn again: a layer on it is unstable (after a move: whole), or its stream changed (the changes' columns
   * in view); its device x range into `span`.
   */
  private redrawn(l: SectionLayout, e: number, moved: boolean): boolean {
    const m = this.model as ViewModel;
    const rv = m.resolved;
    const v = this.v;
    const t = l.entryTrack[e];
    const sp = this.span;
    const L = Math.round(l.dataLeft * v.dpr) | 0;
    const R = Math.round(l.dataRight * v.dpr) | 0;
    if (l.entryKind[e] !== TRACK || t < 0) return false;
    sp[0] = R;
    sp[1] = L;
    for (let k = l.entryLayer[e]; k < rv.layers.count && rv.layers.track[k] === t; k++) {
      const r = this.renderers[k];
      const st = r?.stable;
      // A band's rows layer copies (`bandCopies`): only the other layers on it count.
      const band = k === l.rows.layer && banded(this, l, e);
      const whole =
        moved && r && !band && !(typeof st === 'function' ? st(m.prepared[k].state, rv.layers.params[k]) : st);
      if (whole || (this.fresh[k] === 1 && (this.cols[2 * k] < 0 || l.geometry !== LINEAR))) {
        sp[0] = L;
        sp[1] = R;
        return true;
      }
      if (this.fresh[k] === 0) continue;
      for (let b = 0; v.axis && b < v.axis.count; b++) {
        if (v.axis.stacked && b !== l.entryBlock[e]) continue;
        const n = v.runToShown(b, this.cols[2 * k], this.cols[2 * k + 1], this.shown);
        for (let q = 0; q < n; q++) {
          sp[0] = Math.max(L, Math.min(sp[0], v.xd(b, this.shown[2 * q])));
          sp[1] = Math.min(R, Math.max(sp[1], v.xd(b, this.shown[2 * q + 1])));
        }
      }
    }
    return sp[1] > sp[0];
  }
  /** The device px cell edge at or beyond `px` in direction `dir`, on the block under it. */
  private cellEdge(px: number, dir: number): number {
    const v = this.v;
    const ax = v.axis;
    for (let b = 0; ax && b < ax.count; b++) {
      const k = (px / v.dpr - v.blockX[b]) / (v.cell * v.span[b]);
      if (k < 0 || k > ax.shown[b]) continue;
      return v.xd(b, dir > 0 ? Math.ceil(k) : Math.floor(k));
    }
    return px;
  }
  /** The view's shown range per block cut to device px [x0, x1), a position either side kept (linear only). */
  private strip(x0: number, x1: number) {
    const v = this.v;
    const n = v.visLo.length;
    if (this.stripLo.length !== n) {
      this.stripLo = new Int32Array(n);
      this.stripHi = new Int32Array(n);
    }
    for (let b = 0; b < n; b++) {
      const w = v.cell * v.span[b];
      const lo = Math.floor((x0 / v.dpr - v.blockX[b]) / w) - 1;
      const hi = Math.ceil((x1 / v.dpr - v.blockX[b]) / w) + 1;
      this.stripLo[b] = Math.max(v.visLo[b], Math.min(v.visHi[b], lo));
      this.stripHi[b] = Math.max(this.stripLo[b], Math.min(v.visHi[b], hi));
    }
    v.visLo = this.stripLo;
    v.visHi = this.stripHi;
  }
}

/** Wrapped: entry e is a chip's after another on its row (its line shares the previous line's row). */
function chipAfter(l: SectionLayout, e: number): boolean {
  return l.geometry === WRAPPED_GEOMETRY && l.entryRow[e] >= 0 && l.lines.follows[l.entryRow[e]] === 1;
}
