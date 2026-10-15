/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The circular geometry's backend (seqspec:viewer#the-render-backend): a section on rings is drawn as its axis unrolled
// (ring.ts), and this wrapper maps each batch onto the circle before the real backend draws it, so every built-in draws
// on a ring as it draws linearly, its reductions per device column becoming reductions per angular column. Rects become
// annular sectors (a stroked one its four edges), path points and glyph and text anchors move to their polar places (a
// line long in
// x is cut so it bends with the ring), images are left out. Clips and partial clears are the unrolled frame's, so a
// ring pass clips to the whole pass and draws whole. Allocation-free: its batches are reused and grow only.

import type { Ring } from '../ring.js';
import {
  type BackendImages,
  GlyphBatch,
  PathBatch,
  PathVerb,
  PathVerbCoords,
  type RectBatch,
  type RenderBackend,
  SectorBatch,
  TextAlign,
  TextBaseline,
  TextBatch,
} from './backend.js';

/** A label group's baselines: beside the ring, below its anchor (at the top), above it (at the bottom). */
const BASES = [TextBaseline.middle, TextBaseline.bottom, TextBaseline.top] as const;
/** Device px of unrolled x a path line is cut into, so it follows the ring. */
const STEP = 6;

export class RingBackend implements RenderBackend {
  /** The backend drawn into. */
  to: RenderBackend = null as unknown as RenderBackend;
  private readonly s = new SectorBatch();
  private readonly p = new PathBatch();
  private readonly g = new GlyphBatch();
  private readonly t = new TextBatch();
  private readonly pt = new Float64Array(2);
  // The map, device px: angle = a0 + ka·(x − x0); radius = r0 − (y − y0)·k.
  private a0 = 0;
  private ka = 0;
  private x0 = 0;
  private r0 = 0;
  private y0 = 0;
  private k = 1;
  private d = 1;
  /** The outermost ring's inner radius (`Ring.edge`), and the outer radius of the ink drawn in it this pass. */
  private edge = 0;
  private ink = 0;
  private cx = 0;
  private cy = 0;
  private w = 0;
  private h = 0;

  get caps(): number {
    return this.to.caps;
  }
  get images(): BackendImages {
    return this.to.images;
  }
  /** Draws through `ring`'s map (onto `to`, set per frame); at a new map only. */
  map(ring: Ring): this {
    const d = ring.dpr;
    this.a0 = ring.theta0;
    this.ka = (2 * Math.PI) / (ring.w * d);
    this.x0 = ring.left * d;
    this.r0 = ring.r * d;
    this.y0 = ring.top * d;
    this.k = ring.scale;
    this.d = d;
    this.edge = ring.edge;
    this.cx = ring.cx;
    this.cy = ring.cy;
    return this;
  }

  beginPass(pass: number, width: number, height: number, dpr: number): void {
    this.ink = 0;
    this.w = width;
    this.h = height;
    this.to.beginPass(pass, width, height, dpr);
  }
  endPass(): void {
    this.to.endPass();
  }
  clear(x: number, y: number, w: number, h: number): void {
    if (x <= 0 && y <= 0 && x + w >= this.w && y + h >= this.h) this.to.clear(x, y, w, h);
  }
  pushClip(): void {
    this.to.pushClip(0, 0, this.w, this.h);
  }
  popClip(): void {
    this.to.popClip();
  }
  copy(): void {}
  present(x: number, y: number, w: number, h: number): void {
    this.to.present(x, y, w, h);
  }
  stale(sx: number, sy: number, sw: number, sh: number, dx: number, dy: number, dw: number, dh: number): void {
    this.to.stale(sx, sy, sw, sh, dx, dy, dw, dh);
  }

  rects(b: RectBatch): void {
    const s = this.s;
    s.reset();
    for (let i = 0; i < b.n; i++) {
      s.begin(b.color[i], b.alpha[i]);
      const a0 = this.a0 + this.ka * (b.x[i] - this.x0);
      const a1 = this.a0 + this.ka * (b.x[i] + b.w[i] - this.x0);
      const ro = this.r0 - (b.y[i] - this.y0) * this.k;
      const ri = Math.max(0, this.r0 - (b.y[i] + b.h[i] - this.y0) * this.k);
      const st = b.stroke[i];
      // Ink in the outermost ring (the ruler's ticks): its labels go just past it.
      if (ri >= this.edge && ro > this.ink) this.ink = ro;
      // Filled: one sector; stroked: its outer and inner arcs and its two radial sides, each `st` device px.
      const da = st / Math.max(1, (ri + ro) / 2);
      for (let edge = st === 0 ? 4 : 0; edge < 5; edge++) {
        const j = s.add();
        s.cx[j] = this.cx;
        s.cy[j] = this.cy;
        s.r0[j] = edge === 0 ? Math.max(ri, ro - st) : ri;
        s.r1[j] = edge === 1 ? Math.min(ro, ri + st) : ro;
        s.a0[j] = edge === 3 ? Math.max(a0, a1 - da) : a0;
        s.a1[j] = edge === 2 ? Math.min(a1, a0 + da) : a1;
        if (st !== 0 && edge === 3) break;
      }
    }
    if (s.n > 0) this.to.sectors(s);
  }
  sectors(b: SectorBatch): void {
    this.to.sectors(b);
  }

  paths(b: PathBatch): void {
    const p = this.p;
    p.reset();
    let px = 0;
    let py = 0;
    for (let q = 0; q < b.n; q++) {
      p.begin(b.color[q], b.alpha[q], b.mode[q], b.width[q]);
      if (b.dash[q] > 0) p.dashed(b.dash[q], b.gap[q], b.phase[q]);
      const v1 = q + 1 < b.n ? b.start[q + 1] : b.nv;
      let c = b.cstart[q];
      for (let v = b.start[q]; v < v1; v++) {
        const verb = b.verbs[v];
        const x = b.coords[c];
        const y = b.coords[c + 1];
        if (verb === PathVerb.line) {
          // Cut along x so it bends with the ring; each piece's end mapped.
          const n = Math.max(1, Math.ceil(Math.abs(x - px) / STEP));
          for (let j = 1; j <= n; j++) {
            this.pt[0] = px + ((x - px) * j) / n;
            this.pt[1] = py + ((y - py) * j) / n;
            this.point(p, PathVerb.line);
          }
        } else if (verb === PathVerb.move) {
          this.pt[0] = x;
          this.pt[1] = y;
          this.point(p, PathVerb.move);
        } else if (verb === PathVerb.close) p.verb(PathVerb.close);
        else this.mapped(p, verb, b.coords, c);
        if (verb !== PathVerb.close && verb !== PathVerb.arc && verb !== PathVerb.roundRect) {
          px = b.coords[c + PathVerbCoords[verb] - 2];
          py = b.coords[c + PathVerbCoords[verb] - 1];
        }
        c += PathVerbCoords[verb];
      }
    }
    if (p.n > 0) this.to.paths(p);
  }
  /** The unrolled point in `pt` (no double crosses the call), mapped, as `verb`'s coordinates. */
  private point(p: PathBatch, verb: number) {
    const o = p.verb(verb);
    const a = this.a0 + this.ka * (this.pt[0] - this.x0);
    const r = Math.max(0, this.r0 - (this.pt[1] - this.y0) * this.k);
    p.coords[o] = this.cx + r * Math.cos(a);
    p.coords[o + 1] = this.cy + r * Math.sin(a);
  }
  /** A curve's points mapped one by one; an arc's or a round rect's anchor moved, its size kept. */
  private mapped(p: PathBatch, verb: number, from: Float64Array, c: number) {
    const o = p.verb(verb);
    const n = PathVerbCoords[verb];
    for (let j = 0; j < n; j++) p.coords[o + j] = from[c + j];
    const pairs = verb === PathVerb.quad || verb === PathVerb.cubic ? n / 2 : 1;
    for (let j = 0; j < pairs; j++) {
      const a = this.a0 + this.ka * (from[c + 2 * j] - this.x0);
      const r = Math.max(0, this.r0 - (from[c + 2 * j + 1] - this.y0) * this.k);
      p.coords[o + 2 * j] = this.cx + r * Math.cos(a);
      p.coords[o + 2 * j + 1] = this.cy + r * Math.sin(a);
    }
  }

  glyphs(b: GlyphBatch): void {
    const g = this.g;
    g.reset();
    g.append(b, 0, b.n);
    for (let i = 0; i < g.n; i++) {
      const a = this.a0 + this.ka * (g.x[i] - this.x0);
      const r = Math.max(0, this.r0 - (g.y[i] - this.y0) * this.k);
      g.x[i] = this.cx + r * Math.cos(a);
      g.y[i] = this.cy + r * Math.sin(a);
    }
    if (g.n > 0) this.to.glyphs(g);
  }
  /**
   * Text anchors at their polar places, upright. Text in the outermost ring (the ruler's numbers) is the ring's label:
   * set just outside the ring, aligned away from the center by its angle, so it never covers the ring.
   */
  text(b: TextBatch): void {
    let inside = 0;
    for (let i = 0; i < b.n; i++) if (this.r0 - (b.y[i] - this.y0) * this.k < this.edge) inside++;
    if (inside > 0) this.place(b, -1);
    // Outside: a batch per (alignment, baseline) the angles ask.
    for (let g = 0; inside < b.n && g < 9; g++) this.place(b, g);
  }
  /** Items of `b` mapped into a batch: those inside the outer ring's edge (g < 0), or outside in label group g. */
  private place(b: TextBatch, g: number) {
    const t = this.t;
    t.reset();
    const out = g >= 0;
    t.style(b.font, out ? g % 3 : b.align, out ? BASES[(g / 3) | 0] : b.baseline, b.angle);
    for (let i = 0; i < b.n; i++) {
      const a = this.a0 + this.ka * (b.x[i] - this.x0);
      const c = Math.cos(a);
      const s = Math.sin(a);
      let r = Math.max(0, this.r0 - (b.y[i] - this.y0) * this.k);
      if (r >= this.edge !== out) continue;
      if (out) {
        // Right of the anchor on the right half, left on the left; above it at the top, below at the bottom.
        const ga = c > 0.35 ? TextAlign.left : c < -0.35 ? TextAlign.right : TextAlign.center;
        const gb = s < -0.35 ? 1 : s > 0.35 ? 2 : 0;
        if (ga !== g % 3 || gb !== ((g / 3) | 0)) continue;
        r = (this.ink > 0 ? this.ink : this.r0) + 5 * this.d;
      }
      const j = t.add(b.text[i]);
      t.color[j] = b.color[i];
      t.alpha[j] = b.alpha[i];
      t.sx[j] = b.sx[i];
      t.sy[j] = b.sy[i];
      t.x[j] = this.cx + r * c;
      t.y[j] = this.cy + r * s;
    }
    if (t.n > 0) this.to.text(t);
  }
  /** Images (tiles, miniatures) have no polar form here: left out. */
  image(): void {}
}
