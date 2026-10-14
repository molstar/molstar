/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// (its `mini`, `miniList`, `miniDraw` and `mapTracks`)
// The section's tracks in miniature (seqspec:viewer#extension-points `seqv:map`, the navigators' `miniature`), the
// composer's: each track's layers (the named tracks, else every layer whose renderer is not `miniature: false`) in a
// lane of an image, over the whole axis laid out as a map (`MapView`: the navigator's `navX0s` and `navCell`, never the
// zoom), through the renderer's `overview`, else its `draw`. Images are `KeyedImage`s under the map's size, the layout,
// the regions, the fonts and a generation bumped with the model or the skin: a pan or a zoom builds none, and nothing
// to draw is kept as none. Edge-rate but `stamp`, so it may allocate.

import type { Axis, SectionAxisView } from "../axis.js";
import type { SeqViewerContext } from "../context.js";
import type { Skin, ViewModel } from "../graph.js";
import type { LayerRenderer } from "../registry.js";
import { Lod } from "../viewport.js";
import type { RenderBackend } from "./backend.js";
import { KeyedImage, mapSize, type NavView } from "./nav.js";
import { DrawView, type TracksHost } from "./view.js";

/**
 * A draw view laid out as a map of the whole axis: every shown position of each block, block b from `x0s[b] − left` at
 * `cell` px a cell, over [0, `width`) CSS px; nothing scrolled, no lanes in view. The miniatures draw through it, and
 * so does the vertical overview's map (`rows.overview` takes its tile levels from the section's overview column,
 * `overviewLevels`, whatever box it is drawn into).
 */
export class MapView {
  readonly v = new DrawView();
  private x0s = new Float64Array(0);
  private lo = new Int32Array(0);
  private hi = new Int32Array(0);
  of(av: SectionAxisView, x0s: Float64Array, left: number, cell: number, width: number): DrawView {
    const ax = av.axis as Axis;
    if (this.x0s.length !== ax.count) {
      this.x0s = new Float64Array(ax.count);
      this.lo = new Int32Array(ax.count);
      this.hi = new Int32Array(ax.count);
    }
    for (let b = 0; b < ax.count; b++) {
      this.x0s[b] = x0s[b] - left;
      this.hi[b] = ax.shown[b];
    }
    const v = this.v.from(av);
    v.blockX = v.blockX0 = this.x0s;
    v.visLo = v.viewLo = this.lo;
    v.visHi = v.viewHi = this.hi;
    v.scrollDev = 0;
    v.cell = cell;
    v.lod = cell >= 3 ? Lod.cells : cell >= 1 ? Lod.thin : Lod.aggregate;
    v.left = 0;
    v.right = width;
    v.miniature = true;
    return v;
  }
}

const NONE: readonly string[] = [];
const size = new Int32Array(2);

export class Miniature implements TracksHost {
  private readonly ctx: SeqViewerContext;
  private readonly section: string;
  private readonly map = new MapView();
  private readonly layers: number[] = [];
  /** Per prepared state its renderer's cache, and the layers that threw: dropped with the model or the skin. */
  private readonly caches = new Map<unknown, { r: LayerRenderer; c: unknown }>();
  private readonly failed = new Set<number>();
  private readonly key = new KeyedImage(8);
  private model: ViewModel | null = null;
  private skin: Skin | null = null;
  /** Bumped when the model or the skin changes (a key's number for both). */
  private gen = 0;

  constructor(ctx: SeqViewerContext, section: string) {
    this.ctx = ctx;
    this.section = section;
  }

  /** Into `out[at..at + 3]`: the layout, regions, fonts and generation versions (a new model or skin drops the caches). */
  stamp(out: Float64Array, at: number): void {
    const ctx = this.ctx;
    const m = ctx.graph.assembled.value?.view ?? null;
    if (m !== this.model || ctx.skin.value !== this.skin) this.reset(m, ctx.skin.value);
    out[at] = ctx.layout.get(this.section)?.version.value ?? -1;
    out[at + 1] = ctx.streams.version.value;
    out[at + 2] = ctx.text.version.value;
    out[at + 3] = this.gen;
  }

  /** `seqv:map`'s layer image (its renderer keeps it): the tracks over the map `v` lays out. */
  draw(into: RenderBackend, refs: readonly string[], v: DrawView, w: number, h: number): boolean {
    if (!this.pick(refs)) return false;
    this.render(into, v.navX0s, v.left, v.navCell, w, h);
    return true;
  }

  /** The navigator's map image for `n` (its band, its `tracks`), built when its key changed; -1: none to show. */
  image(n: NavView): number {
    const p = n.parts;
    if (p.onRuler || !n.axis || !(p.axis.includes("overview") || p.axis.includes("miniature"))) return -1;
    const images = n.backend.images;
    const k = this.key;
    mapSize(n, size);
    k.next[0] = size[0];
    k.next[1] = size[1];
    k.next[2] = n.dpr;
    k.next[3] = n.cell;
    this.stamp(k.next, 4);
    if (!k.stale(images)) return k.id;
    if (!this.pick(p.tracks ?? NONE)) k.keep(images, -1);
    else {
      this.render(images.beginLayer(size[0], size[1]), n.x0s, n.left, n.cell, size[0], size[1]);
      k.keep(images, images.endLayer());
    }
    return k.id;
  }

  dispose(): void {
    this.key.dispose();
    this.reset(null, null);
  }

  /** The layers to draw into `this.layers`; false: none. */
  private pick(refs: readonly string[]): boolean {
    const ctx = this.ctx;
    const m = this.model;
    this.layers.length = 0;
    if (!m || !ctx.layout.axis(this.section)?.axis) return false;
    const rv = m.resolved;
    const section = ctx.layout.sectionIndex(this.section);
    for (let k = 0; k < rv.layers.count; k++) {
      const t = rv.layers.track[k];
      if (rv.tracks.section[t] !== section || rv.tracks.hidden[t] || this.failed.has(k)) continue;
      const r = ctx.registry.get("renderers", m.prepared[k].renderer);
      if (!r || !(r.draw || r.overview)) continue;
      if (refs.length > 0 ? refs.includes(m.set.strings.get(rv.tracks.ref[t])) : r.miniature !== false)
        this.layers.push(k);
    }
    return this.layers.length > 0;
  }
  private reset(m: ViewModel | null, skin: Skin | null): void {
    for (const e of this.caches.values()) e.r.cache?.dispose?.(e.c);
    this.caches.clear();
    this.failed.clear();
    this.model = m;
    this.skin = skin;
    this.gen++;
  }

  /** The picked layers in lanes over a `w` × `h` device px layer, the blocks from `x0s` (surface px, the map from `left`). */
  private render(into: RenderBackend, x0s: Float64Array, left: number, cell: number, w: number, h: number): void {
    const ctx = this.ctx;
    const m = this.model as ViewModel;
    const av = ctx.layout.axis(this.section) as SectionAxisView;
    const v = this.map.of(av, x0s, left, cell, w / av.dpr);
    v.backend = into;
    v.skin = ctx.skin.value;
    v.props = ctx.layout.props(this.section);
    v.text = ctx.text;
    const lane = h / v.dpr / this.layers.length;
    for (let i = 0; i < this.layers.length; i++) {
      const k = this.layers[i];
      const p = m.prepared[k];
      const r = ctx.registry.get("renderers", p.renderer) as LayerRenderer;
      const params = m.resolved.layers.params[k];
      v.layer = k;
      v.setThemes(p.themes);
      v.stream = ctx.streams.get(m.set.strings.get(m.resolved.layers.ref[k]));
      v.cache = r.cache ? this.cacheOf(r, p.state) : null;
      v.top = i * lane;
      v.height = lane;
      into.pushClip(0, Math.round(v.top * v.dpr), w, Math.round(lane * v.dpr));
      try {
        if (!r.overview?.(v, p.state, params)) r.draw?.(v, p.state, params);
      } catch (err) {
        this.failed.add(k);
        ctx.report("seqv.renderer.error", `layer ${k} threw while drawn in miniature and is omitted: ${String(err)}`);
      } finally {
        v.cache = v.stream = null;
        into.popClip();
      }
    }
  }
  private cacheOf(r: LayerRenderer, s: unknown): unknown {
    let e = this.caches.get(s);
    if (!e) {
      e = { r, c: r.cache?.create(s, this.ctx.skin.value) };
      this.caches.set(s, e);
    }
    return e.c;
  }
}
