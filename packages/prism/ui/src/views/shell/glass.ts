/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { WordmarkGeometry } from '../../components/icon/wordmark.js';

// The splash's wordmark as glass (`prism-field.ts` draws it): the word cut into triangular shards, which fly in and
// snap together when the splash starts, and fly apart when it is pressed, spinning and slowing, then come back
// together. The shards are drawn as the wordmark is (the same path data, line and the logo's gradient), each its piece
// of it, its edges glinting while it moves; whole, the page's own wordmark shows instead (`onWhole`). In drawing units
// (the wordmark's viewBox: 1 unit a pixel at 1×); everything per frame is preallocated.

/** Columns and rows of the grid the word is cut on (each cell two triangles, the inner corners moved). */
const Cols = 11;
const Rows = 4;
const Shards = Cols * Rows * 2;
/** Seconds a shard takes to come home, and how far apart their starts are. */
const Assemble = 1;
const Stagger = 0.35;
/** Seconds the pieces stay apart before they come back. */
const Hold = 0.8;
/** How fast the pieces slow when apart (per second). */
const Drag = 3;

export type GlassPhase = 'assembling' | 'whole' | 'broken';

/** A small seeded generator (mulberry32). */
function generator(seed: number): () => number {
  let a = seed >>> 0 || 1;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** easeOutBack: past home a little, then settled: a piece snapping into place. */
function snap(p: number): number {
  const c = 1.4;
  const q = p - 1;
  return 1 + (c + 1) * q * q * q + c * q * q;
}

export class Glass {
  readonly #random: () => number;
  readonly #word: Path2D[];
  // The shards: their outline (around their centroid) and centroid, where they are (offset from home, rotation), how
  // they move apart, where an assembly started from, and its stagger.
  readonly #shard: Path2D[] = [];
  readonly #cx = new Float32Array(Shards);
  readonly #cy = new Float32Array(Shards);
  readonly #ox = new Float32Array(Shards);
  readonly #oy = new Float32Array(Shards);
  readonly #rot = new Float32Array(Shards);
  readonly #vx = new Float32Array(Shards);
  readonly #vy = new Float32Array(Shards);
  readonly #vr = new Float32Array(Shards);
  readonly #fx = new Float32Array(Shards);
  readonly #fy = new Float32Array(Shards);
  readonly #fr = new Float32Array(Shards);
  readonly #delay = new Float32Array(Shards);
  #phase: GlassPhase = 'whole';
  #since = 0;
  #time = 0;
  /** Told when the glass becomes whole (the page's wordmark shows) or stops being whole. */
  onWhole: (whole: boolean) => void = () => {};
  /** Told when it snaps back together after a break (the light pops). */
  onSnap: () => void = () => {};

  constructor(seed = 11) {
    this.#random = generator(seed);
    this.#word = WordmarkGeometry.paths.map((d) => new Path2D(d));
    // The grid over the word (with its line), its inner corners moved: broken glass, not a lattice.
    const x0 = -1;
    const y0 = -1;
    const w = WordmarkGeometry.width;
    const h = WordmarkGeometry.height;
    const px: number[] = [];
    const py: number[] = [];
    for (let r = 0; r <= Rows; r++) {
      for (let c = 0; c <= Cols; c++) {
        const inner = r > 0 && r < Rows && c > 0 && c < Cols;
        px.push(x0 + (c / Cols) * w + (inner ? (this.#random() - 0.5) * 0.7 * (w / Cols) : 0));
        py.push(y0 + (r / Rows) * h + (inner ? (this.#random() - 0.5) * 0.7 * (h / Rows) : 0));
      }
    }
    const at = (c: number, r: number) => r * (Cols + 1) + c;
    let k = 0;
    const add = (a: number, b: number, c: number) => {
      const cx = (px[a] + px[b] + px[c]) / 3;
      const cy = (py[a] + py[b] + py[c]) / 3;
      const p = new Path2D();
      p.moveTo(px[a] - cx, py[a] - cy);
      p.lineTo(px[b] - cx, py[b] - cy);
      p.lineTo(px[c] - cx, py[c] - cy);
      p.closePath();
      this.#shard.push(p);
      this.#cx[k] = cx;
      this.#cy[k] = cy;
      k++;
    };
    for (let r = 0; r < Rows; r++) {
      for (let c = 0; c < Cols; c++) {
        const [a, b, d, e] = [at(c, r), at(c + 1, r), at(c, r + 1), at(c + 1, r + 1)];
        // Each cell split along one diagonal or the other.
        if (this.#random() < 0.5) {
          add(a, b, e);
          add(a, e, d);
        } else {
          add(a, b, d);
          add(b, e, d);
        }
      }
    }
  }

  get phase(): GlassPhase {
    return this.#phase;
  }

  /** Whole: the star splits the light. */
  get whole(): boolean {
    return this.#phase === 'whole';
  }

  /** Scattered around the word, then flying in: the splash's start. */
  scatter(): void {
    for (let i = 0; i < Shards; i++) {
      const a = this.#random() * Math.PI * 2;
      const d = 14 + this.#random() * 30;
      this.#fx[i] = Math.cos(a) * d;
      this.#fy[i] = Math.sin(a) * d * 0.7;
      this.#fr[i] = (this.#random() * 2 - 1) * 2.4;
      this.#delay[i] = this.#random() * Stagger;
    }
    this.#begin('assembling');
  }

  /** Pressed at (x, y) (drawing units): the pieces fly apart from there, faster the nearer, spinning. */
  shatter(x: number, y: number): void {
    for (let i = 0; i < Shards; i++) {
      const dx = this.#cx[i] + this.#ox[i] - x;
      const dy = this.#cy[i] + this.#oy[i] - y;
      const d = Math.hypot(dx, dy) || 1;
      const v = (26 + this.#random() * 34) / (1 + d * 0.04);
      this.#vx[i] = (dx / d) * v + (this.#random() - 0.5) * 6;
      this.#vy[i] = (dy / d) * v + (this.#random() - 0.5) * 6;
      this.#vr[i] = (this.#random() * 2 - 1) * 6;
    }
    this.#begin('broken');
  }

  /** Whether (x, y) (drawing units) is on the word (its box, a little around). */
  hits(x: number, y: number): boolean {
    return x >= -2 && x <= WordmarkGeometry.width + 1 && y >= -2 && y <= WordmarkGeometry.height + 1;
  }

  #begin(phase: GlassPhase): void {
    const wasWhole = this.#phase === 'whole';
    this.#phase = phase;
    this.#since = this.#time;
    if (wasWhole && phase !== 'whole') this.onWhole(false);
  }

  /** Back together from where the pieces are now. */
  #gather(): void {
    for (let i = 0; i < Shards; i++) {
      this.#fx[i] = this.#ox[i];
      this.#fy[i] = this.#oy[i];
      this.#fr[i] = this.#rot[i];
      this.#delay[i] = this.#random() * Stagger;
    }
    this.#begin('assembling');
  }

  step(dt: number): void {
    this.#time += dt;
    const since = this.#time - this.#since;
    if (this.#phase === 'broken') {
      const drag = Math.exp(-Drag * dt);
      for (let i = 0; i < Shards; i++) {
        this.#ox[i] += this.#vx[i] * dt;
        this.#oy[i] += this.#vy[i] * dt;
        this.#rot[i] += this.#vr[i] * dt;
        this.#vx[i] *= drag;
        this.#vy[i] *= drag;
        this.#vr[i] *= drag;
      }
      if (since >= Hold) this.#gather();
    } else if (this.#phase === 'assembling') {
      let done = true;
      for (let i = 0; i < Shards; i++) {
        const p = Math.min(1, Math.max(0, (since - this.#delay[i]) / Assemble));
        if (p < 1) done = false;
        const k = 1 - snap(p);
        this.#ox[i] = this.#fx[i] * k;
        this.#oy[i] = this.#fy[i] * k;
        this.#rot[i] = this.#fr[i] * k;
      }
      if (done) {
        this.#ox.fill(0);
        this.#oy.fill(0);
        this.#rot.fill(0);
        this.#phase = 'whole';
        this.onWhole(true);
        this.onSnap();
      }
    }
  }

  /**
   * Draws the pieces (nothing while whole: the page's wordmark shows): the word's line in `stroke` (the logo's gradient,
   * made in drawing units), each piece its part of it, and the piece's edges in `edge`, glinting while it moves. The
   * context is in canvas CSS pixels; the word's drawing origin is at (x, y), `scale` px a unit.
   */
  draw(
    ctx: CanvasRenderingContext2D,
    x: number,
    y: number,
    scale: number,
    stroke: CanvasGradient | string,
    edge: string,
    alpha: number,
  ): void {
    if (this.#phase === 'whole') return;
    const since = this.#time - this.#since;
    const line = WordmarkGeometry.stroke(scale);
    for (let i = 0; i < Shards; i++) {
      const off = Math.hypot(this.#ox[i], this.#oy[i]) + Math.abs(this.#rot[i]) * 6;
      // A shard flying in fades in as it starts.
      const shown =
        this.#phase === 'assembling' ? Math.min(1, Math.max(0, (since - this.#delay[i]) / (Assemble * 0.3))) : 1;
      if (shown <= 0) continue;
      ctx.save();
      ctx.translate(x + (this.#cx[i] + this.#ox[i]) * scale, y + (this.#cy[i] + this.#oy[i]) * scale);
      ctx.rotate(this.#rot[i]);
      ctx.scale(scale, scale);
      const shard = this.#shard[i];
      // Its edges: a glint while it is away from home.
      const glint = Math.min(1, off / 8) * 0.22 * alpha * shown;
      if (glint > 0.01) {
        ctx.globalAlpha = glint;
        ctx.strokeStyle = edge;
        ctx.lineWidth = 0.35;
        ctx.stroke(shard);
      }
      ctx.clip(shard);
      ctx.translate(-this.#cx[i], -this.#cy[i]);
      ctx.globalAlpha = alpha * shown;
      ctx.strokeStyle = stroke;
      ctx.lineWidth = line;
      ctx.lineJoin = 'miter';
      for (let w = 0; w < this.#word.length; w++) ctx.stroke(this.#word[w]);
      ctx.restore();
    }
    ctx.globalAlpha = 1;
  }
}
