/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section surface (seqspec:viewer#the-render-backend): a section's canvases, size, DPR and schedule, driving the
// DOM-free `SectionComposer` over the Canvas2D backend. It owns pixels, never state. Every subscription
// ORs `Dirty` flags and asks for one frame; the frame reads the layout (built if stale), sizes the canvases from it and
// hands the flags to the composer, which draws only what changed. One `ChannelMarks` per surface, synced by the
// composer once a frame. A canvas whose size changes loses its image, so the canvases are sized only inside a frame,
// from the layout the composer draws, and that frame draws every pass whole: the new size and its whole image land in
// one task, and the page never paints a blank or half-drawn canvas. A resize is drawn in the ResizeObserver's callback,
// before the paint of the frame whose layout changed the host. Per event and per frame nothing is allocated: pre-bound
// callbacks, numbers in typed arrays (the draw pattern).

import { ChannelMarks } from "../channels.js";
import { sectionLabel } from "../chrome.js";
import type { SeqViewerContext } from "../context.js";
import { type DrawResources, Pass } from "../draw/backend.js";
import { Dirty, SectionComposer } from "../draw/composer.js";
import { drawResources } from "../draw/resources.js";
import { VerticalNav } from "../draw/vnav.js";
import type { SectionLayout } from "../layout.js";
import { ringWindowTarget } from "../ring-window.js";
import type { Live, Unsubscribable } from "../state.js";
import type { SectionViewport } from "../viewport.js";
import { railGrown, WindowPart } from "../windows.js";
import { Canvas2DBackend, Canvas2DShared, type Paint2D } from "./backend.js";
import { DomInput } from "./dom-input.js";
import { joinShared } from "./shared.js";

export interface SurfaceOptions {
  /** Device pixels per CSS px, fixed; default the window's, followed as it changes (`watchDpr`). */
  dpr?: number;
  /** False: keep the context's text measurer (seqspec:viewer#text-measuring). */
  measurer?: boolean;
  /** Canvases read back often (pixel tests): rasterized on the CPU from the start, as the composer's scratch is. */
  readback?: boolean;
  /**
   * False: every letter through `fillText`, not the glyph atlas (default true). The atlas is the context's, shared by
   * its surfaces: the first surface mounted on a context decides it.
   */
  atlas?: boolean;
}
export interface Surface {
  readonly element: HTMLElement;
  /** Draws what is pending now, synchronously (tests, benchmarks, hidden pages). */
  renderNow(): void;
  dispose(): void;
}
export interface SurfaceStats {
  /** Frames drawn, and those that sized a canvas (each drawn whole). */
  frames: number;
  resizes: number;
  /** Main-thread ms of the last frame drawn (sizing and the composer's draws; the raster after it not included). */
  lastMs: number;
}

export function mountSection(
  ctx: SeqViewerContext,
  section: string,
  host: HTMLElement,
  opts: SurfaceOptions = {},
): SectionSurface {
  return new SectionSurface(ctx, section, host, opts);
}

/** What `watchDpr` reads: a window. */
export type DprHost = Pick<Window, "matchMedia" | "devicePixelRatio">;
/**
 * Calls `changed` with the new device pixel ratio whenever it changes, with or without a size change (a window moved to
 * another display, a browser zoom): a `resolution` query of the current ratio, re-armed at each change.
 */
export function watchDpr(win: DprHost, changed: (dpr: number) => void): () => void {
  let query: MediaQueryList | null = null;
  const arm = () => {
    query?.removeEventListener("change", fire);
    query = win.matchMedia(`(resolution: ${win.devicePixelRatio}dppx)`);
    query.addEventListener("change", fire);
  };
  const fire = () => {
    arm();
    changed(win.devicePixelRatio);
  };
  arm();
  return () => {
    query?.removeEventListener("change", fire);
    query = null;
  };
}

/** Passes with a canvas of their own, in stacking order (`Pass`): the section's four, then the vertical overview's two. */
const PASSES = 6;
/** `frameId` while a frame reads the layout: what it triggers is ORed in and drawn by this frame. */
const DRAWING = -2;

/** One per context, shared by its surfaces (`Canvas2DShared`: images, the scratch canvas, path glyphs). */
const shared = new WeakMap<SeqViewerContext, { s: Canvas2DShared; n: number }>();
function joinImages(ctx: SeqViewerContext, atlas: boolean | undefined): Canvas2DShared {
  const held = shared.get(ctx);
  if (held) {
    held.n++;
    return held.s;
  }
  // A font's px by id, parsed once from its string (ids are never reused).
  let px = new Float64Array(16).fill(Number.NaN);
  const base = drawResources(ctx);
  const res: DrawResources = {
    ...base,
    fontPx: (id) => {
      if (id >= px.length) {
        const grown = new Float64Array(2 * id + 2).fill(Number.NaN);
        grown.set(px);
        px = grown;
      }
      if (Number.isNaN(px[id])) px[id] = base.fontPx(id);
      return px[id];
    },
  };
  // A font load bumps the text version: the glyph atlas is dropped then.
  const s = new Canvas2DShared(res, { atlas, fonts: () => ctx.text.version.value });
  shared.set(ctx, { s, n: 1 });
  return s;
}
function leaveImages(ctx: SeqViewerContext) {
  const held = shared.get(ctx);
  if (held && --held.n === 0) shared.delete(ctx);
}

export class SectionSurface implements Surface {
  readonly element: HTMLDivElement;
  readonly section: string;
  readonly stats: SurfaceStats = { frames: 0, resizes: 0, lastMs: 0 };
  /** The surface's channels, one per surface, synced by the composer once a frame. */
  readonly channels: ChannelMarks;
  readonly composer: SectionComposer;
  private readonly vnav: VerticalNav;
  private readonly ctx: SeqViewerContext;
  private readonly opts: SurfaceOptions;
  /** By `Pass`: main, overlay, nav, nav overlay, vnav, vnav overlay. */
  private readonly canvases: HTMLCanvasElement[] = [];
  private readonly paints: Paint2D[] = [];
  private readonly backend: Canvas2DBackend;
  /** Per canvas as sized: width, height (device px). */
  private readonly sized = new Int32Array(2 * PASSES);
  /** Where the vertical overview's canvases lie in the section, device px: left, top. */
  private readonly placed = new Int32Array(2);
  /** The host as last measured and given to the viewport: width, height (CSS px), DPR. */
  private readonly size = new Float64Array(3);
  /** The rail is grown or the pointer is in the navigator's row (1), and the window part lit here (-1: none). */
  private readonly lit = new Int32Array([0, -1]);
  private dirty: number = Dirty.all;
  private frameId = -1;
  private disposed = false;
  private vp: SectionViewport | null = null;
  private live: Live<SectionLayout> | null = null;
  private readonly bound: Unsubscribable[] = [];
  /** On rings, the viewport the ring's window shows (ring-window.ts) and its subscription: its moves redraw the overlay. */
  private windowOf: unknown = null;
  private windowSub: Unsubscribable | null = null;
  private readonly subs: Unsubscribable[] = [];
  private readonly cleanup: (() => void)[] = [];

  constructor(ctx: SeqViewerContext, section: string, host: HTMLElement, opts: SurfaceOptions = {}) {
    if (ctx.disposed) throw new Error("SectionSurface: the context is disposed");
    this.ctx = ctx;
    this.section = section;
    this.opts = opts;
    const el = document.createElement("div");
    // One tab stop, a named group with the role description, wherever it is mounted (seqspec:interaction#accessibility).
    el.tabIndex = 0;
    el.setAttribute("role", "group");
    el.setAttribute("aria-roledescription", "sequence view");
    el.style.cssText =
      "position:relative;width:100%;height:100%;overflow:hidden;outline:none;touch-action:none;user-select:none";
    host.appendChild(el);
    this.element = el;
    for (let p = 0; p < PASSES; p++) {
      const c = document.createElement("canvas");
      // 0 × 0 until sized (a canvas is 300 × 150 by default): the vertical overview's stay empty while the rows fit.
      c.width = c.height = 0;
      c.style.cssText = "position:absolute;left:0;top:0;width:0;height:0;pointer-events:none";
      el.appendChild(c);
      this.canvases.push(c);
      this.paints.push(c.getContext("2d", { willReadFrequently: !!opts.readback }) as CanvasRenderingContext2D);
    }
    this.backend = new Canvas2DBackend(this.paints, joinImages(ctx, opts.atlas), 1 << Pass.main);
    this.channels = new ChannelMarks(ctx);
    this.composer = new SectionComposer(ctx, section, this.channels);
    this.composer.live = true;
    this.vnav = new VerticalNav(ctx, section, this.channels);
    const s = this.subs;
    // Set only when it changes: a host watching the DOM sees no write for an equal label.
    s.push(
      ctx.document.shown.subscribe(() => {
        const label = sectionLabel(ctx, section);
        if (el.getAttribute("aria-label") !== label) el.setAttribute("aria-label", label);
      }),
    );
    s.push(ctx.viewports.sections.subscribe(this.bind));
    // A context that completes takes its surfaces with it.
    s.push(
      ctx.graph.assembled.subscribe(this.onAll),
      ctx.graph.assembled.subscribe(this.bindWindow),
      ctx.skin.subscribe({ next: this.onAll, complete: this.onComplete }),
    );
    s.push(ctx.hover.live.version.subscribe(this.onOverlay), ctx.selection.current.subscribe(this.onOverlay));
    s.push(ctx.selection.pending.version.subscribe(this.onOverlay), ctx.ui.keyCursor.subscribe(this.onOverlay));
    s.push(ctx.systems.version.subscribe(this.onOverlay), ctx.highlight.version.subscribe(this.onHighlight));
    s.push(ctx.ui.rowHover.subscribe(this.onOverlay));
    s.push(ctx.streams.version.subscribe(this.onStreams));
    s.push(ctx.pointer.version.subscribe(this.onPointer), ctx.ui.drag.subscribe(this.onPointer));
    s.push(ctx.ui.window.version.subscribe(this.onWindow), ctx.ui.active.subscribe(this.onActive));
    s.push(ctx.ui.cursor(section).subscribe(this.onCursor));
    const input = new DomInput(el, ctx.input, this.index);
    this.cleanup.push(joinShared(ctx, el, opts), () => input.stop(!this.ctx.disposed));
    const ro = new ResizeObserver(this.onResize);
    ro.observe(el);
    this.cleanup.push(() => ro.disconnect());
    if (opts.dpr === undefined) this.cleanup.push(watchDpr(window, this.onDpr));
    this.measure(this.dprNow());
  }

  /** ORs `flags` (`Dirty`) in and asks for one frame; none while paused (`ui.active` false) or disposed. */
  invalidate(flags: number): void {
    this.dirty |= flags;
    if (this.frameId === -1 && this.dirty !== 0 && !this.disposed && this.ctx.ui.active.value)
      this.frameId = this.ctx.scheduler.requestFrame(this.onFrame);
  }
  renderNow(): void {
    if (this.frameId >= 0) this.ctx.scheduler.cancelFrame(this.frameId);
    this.frameId = -1;
    if (!this.disposed) this.frame();
  }
  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    if (this.frameId >= 0) this.ctx.scheduler.cancelFrame(this.frameId);
    if (this.animId >= 0) this.ctx.scheduler.cancelFrame(this.animId);
    this.frameId = this.animId = -1;
    for (const s of this.subs) s.unsubscribe();
    for (const s of this.bound) s.unsubscribe();
    this.windowSub?.unsubscribe();
    for (const c of this.cleanup) c();
    this.composer.dispose();
    this.vnav.dispose();
    this.channels.dispose();
    leaveImages(this.ctx);
    for (const c of this.canvases) c.width = c.height = 0;
    this.element.remove();
    this.vp = this.live = null;
  }

  // ---- Bindings and flags -------------------------------------------------------------------------------------------

  /** The section's viewport and layout, rebound when either object is replaced (the section left and came back). */
  private readonly bind = () => {
    const vp = this.ctx.viewports.get(this.section);
    const live = this.ctx.layout.get(this.section);
    if (vp === this.vp && live === this.live) return;
    for (const s of this.bound) s.unsubscribe();
    this.bound.length = 0;
    this.vp = vp;
    this.live = live;
    if (vp && live) {
      const z = this.size;
      if (z[0] > 0 || z[1] > 0) vp.setSize(z[0], z[1], z[2]);
      this.bound.push(vp.live.version.subscribe(this.onAll), live.version.subscribe(this.onAll));
    }
    this.composer.redraw();
    this.vnav.redraw();
    this.bindWindow();
    this.invalidate(Dirty.all);
  };
  /** A ring's window follows the linked section's viewport: rebound when the view or the sections change (edge rate). */
  private readonly bindWindow = () => {
    const vp = this.disposed ? null : ringWindowTarget(this.ctx, this.ctx.layout.sectionIndex(this.section));
    if (vp === this.windowOf) return;
    this.windowSub?.unsubscribe();
    this.windowOf = vp;
    this.windowSub = vp ? vp.live.version.subscribe(this.onOverlay) : null;
  };
  private readonly index = () => this.ctx.layout.sectionIndex(this.section);
  private readonly onComplete = () => this.dispose();
  private readonly onCursor = (c: string) => {
    this.element.style.cursor = c;
  };
  private readonly onAll = () => this.invalidate(Dirty.all);
  private readonly onOverlay = () => this.invalidate(Dirty.overlay);
  /** A highlight of style `band` is drawn in the main pass: the composer redraws the tracks only if one changed. */
  private readonly onHighlight = () => this.invalidate(Dirty.overlay | Dirty.tracks);
  /**
   * Regions arrived, failed or left: the composer redraws the arrived layers' entries, the navigator if its map changed,
   * and the overlay (skeletons, error bands, fade-ins).
   */
  private readonly onStreams = () => this.invalidate(Dirty.all);
  private readonly onActive = (on: boolean) => {
    if (on) this.invalidate(0);
  };
  /** The overlay only when the rail's or the navigator row's hot state changes (the composer's `navOf` test). */
  private readonly onPointer = () => {
    const i = this.ctx.layout.sectionIndex(this.section);
    const l = this.live?.current ?? null;
    const p = this.ctx.pointer.current;
    const d = this.ctx.ui.drag.value;
    const hot = railGrown(d, l, p, i) || (d === null && l !== null && p.inside && p.section === i && p.y < l.navHeight);
    if ((hot ? 1 : 0) === this.lit[0]) return;
    this.lit[0] = hot ? 1 : 0;
    this.invalidate(Dirty.overlay);
  };
  /** The lit window part (`ui.window`), when it changes in this section. */
  private readonly onWindow = () => {
    const w = this.ctx.ui.window.current;
    const at =
      w.section === this.ctx.layout.sectionIndex(this.section) && w.part !== WindowPart.none ? w.kind * 8 + w.part : -1;
    if (at === this.lit[1]) return;
    this.lit[1] = at;
    this.invalidate(Dirty.overlay);
  };

  // ---- Size and DPR -------------------------------------------------------------------------------------------------

  private dprNow(): number {
    return this.opts.dpr ?? window.devicePixelRatio ?? 1;
  }
  private readonly onResize = () => this.measure(this.dprNow());
  private readonly onDpr = (dpr: number) => this.measure(this.opts.dpr ?? dpr);
  /**
   * The host's size and DPR to the viewport, then drawn now: the ResizeObserver calls back after layout and before
   * paint, so the canvases take the new size with their whole image in the frame that resized the host.
   */
  private measure(dpr: number) {
    const w = this.element.clientWidth;
    const h = this.element.clientHeight;
    const z = this.size;
    if (w === z[0] && h === z[1] && dpr === z[2]) return;
    z[0] = w;
    z[1] = h;
    z[2] = dpr;
    this.vp?.setSize(w, h, dpr);
    if (this.ctx.ui.active.value) this.renderNow();
  }
  /**
   * Each canvas at its pass's size, from the layout the composer draws (its `W`, `H`: the layout's CSS px × DPR,
   * rounded up, shown at that size over CSS px); the vertical overview's two lie in the box `vnav` places. True when one
   * changed: the frame then draws every pass whole.
   */
  private sizeCanvases(l: SectionLayout, dpr: number): boolean {
    const W = Math.ceil(l.width * dpr) | 0;
    const H = Math.ceil(l.height * dpr) | 0;
    const N = Math.ceil(l.navHeight * dpr) | 0;
    const box = this.vnav.box;
    this.vnav.place(l, this.ctx.layout.props(this.section), dpr);
    const moved = this.placed[0] !== box[0] || this.placed[1] !== box[1];
    let changed = moved;
    for (let p = 0; p < PASSES; p++) {
      const v = p >= Pass.vnav;
      const w = v ? box[2] : W;
      const h = v ? box[3] : p === Pass.nav || p === Pass.navOverlay ? N : H;
      const c = this.canvases[p];
      if (v && moved) {
        c.style.left = `${box[0] / dpr}px`;
        c.style.top = `${box[1] / dpr}px`;
      }
      if (this.sized[2 * p] === w && this.sized[2 * p + 1] === h) continue;
      c.width = this.sized[2 * p] = w;
      c.height = this.sized[2 * p + 1] = h;
      c.style.width = `${w / dpr}px`;
      c.style.height = `${h / dpr}px`;
      changed = true;
    }
    this.placed[0] = box[0];
    this.placed[1] = box[1];
    return changed;
  }

  // ---- Frames -------------------------------------------------------------------------------------------------------

  /** The next animated overlay's time (ms), and its polling frame (-1: none). */
  private readonly due = new Float64Array(1);
  private animId = -1;
  private readonly onAnim = () => {
    this.animId = -1;
    if (this.disposed || !this.ctx.ui.active.value) return;
    if (this.ctx.scheduler.now() < this.due[0] - 1) {
      this.animId = this.ctx.scheduler.requestFrame(this.onAnim);
      return;
    }
    // Due: drawn in this frame (a frame already asked for takes the flag with its own).
    this.dirty |= Dirty.overlay;
    if (this.frameId === -1) this.frame();
  };
  private readonly onFrame = () => {
    this.frameId = -1;
    if (!this.disposed && this.ctx.ui.active.value) this.frame();
  };
  /** Reads the layout (built if stale), sizes the canvases from it, then the composer draws the flags gathered. */
  private frame() {
    const ctx = this.ctx;
    this.frameId = DRAWING;
    // A moved view's re-hit first: the hover it changes is drawn in this frame, with the view
    // (seqspec:interaction#input-events).
    ctx.input.flush?.();
    const live = this.live === null ? null : ctx.layout.get(this.section);
    const vp = this.vp;
    const model = ctx.graph.assembled.value?.view ?? null;
    this.frameId = -1;
    if (!live || !vp || !model) {
      // Nothing to show: no document, or the section left the view; the flags wait for one.
      if (!model)
        for (let p = 0; p < PASSES; p++) this.paints[p].clearRect(0, 0, this.sized[2 * p], this.sized[2 * p + 1]);
      return;
    }
    // A host of no size (hidden, collapsed) shows nothing, and a canvas of no size cannot be drawn from: the flags wait.
    const dpr = vp.live.current.dpr;
    if (live.current.width * dpr < 1 || live.current.height * dpr < 1) return;
    const t0 = performance.now();
    if (this.sizeCanvases(live.current, dpr)) {
      this.composer.redraw();
      this.vnav.redraw();
      this.dirty = Dirty.all;
      this.stats.resizes++;
    }
    const flags = this.dirty;
    this.dirty = 0;
    this.composer.frame(this.backend, flags);
    this.vnav.frame(this.backend, flags);
    // A draw left work (a dot plot's tiles): the tracks again next frame, that layer's entries alone.
    if (this.composer.again) {
      this.composer.again = false;
      this.invalidate(Dirty.tracks);
    }
    // A skeleton shimmering or a region fading in: the overlay again, next frame for a fade, else at the shimmer's rate
    // (`Stream.ShimmerFps`): frames are polled until it is due, nothing drawn (pre-bound, nothing allocated).
    if (this.composer.animating) {
      this.due[0] = ctx.scheduler.now() + this.composer.animStep[0];
      if (this.animId === -1) this.animId = ctx.scheduler.requestFrame(this.onAnim);
    }
    this.stats.frames++;
    this.stats.lastMs = performance.now() - t0;
  }
}
