/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { WordmarkGeometry } from './wordmark.js';

// The Mol* wordmark as glass: the word cut into triangular shards that drift in from around it and settle together, fly
// apart when pressed (spinning and slowing) and snap back; under the pointer the pieces near it lean away and tilt a
// little, wobbling back on a spring. The shards are drawn as the wordmark is (the same paths, line and the logo's
// gradient), each its piece of it, its edges glinting while it moves. In drawing units (the wordmark's viewBox: 1 unit a
// pixel at 1×). Everything a frame touches is preallocated: a step and a draw allocate nothing.

/** Columns and rows of the grid the word is cut on (each cell two triangles, the inner corners moved). */
const Cols = 11;
const Rows = 4;
const Shards = Cols * Rows * 2;

/** How a gathering goes: seconds a shard takes, and how far apart their starts are. */
export interface GlassTiming {
  readonly duration: number;
  readonly stagger: number;
}

/** The first assembly: pieces drifting in and settling (eased out, no overshoot), their starts spread. */
export const GlassGather: GlassTiming = { duration: 2.4, stagger: 0.7 };
/** After a press: quick, a little past home. */
const SnapBack: GlassTiming = { duration: 0.45, stagger: 0.12 };
/** Seconds the pieces stay apart before they come back. */
const Hold = 0.8;
/** How fast the pieces slow when apart (per second). */
const Drag = 3;
/** The pointer: its reach (units), the most a piece leans (units) and tilts (radians), the spring (a little bouncy). */
const HoverReach = 14;
const HoverLean = 1.4;
const HoverTilt = 0.1;
const Stiffness = 140;
const Damping = 11;
/** Sprites a row in the pieces' atlas. */
const AtlasCols = 11;

export type GlassPhase = 'assembling' | 'whole' | 'broken';

/** A small seeded generator (mulberry32): the same cut for the same seed. */
function generator(seed: number): () => number {
  let a = seed >>> 0 || 1;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** easeOutBack: past home a little, then settled. */
function snap(p: number): number {
  const c = 1.4;
  const q = p - 1;
  return 1 + (c + 1) * q * q * q + c * q * q;
}

/** The word's paths, made once for every glass (they never change). */
let wordPaths: Path2D[] | undefined;
function word(): Path2D[] {
  wordPaths ??= WordmarkGeometry.paths.map((d) => new Path2D(d));
  return wordPaths;
}

/**
 * The cut (the shards' outlines, middles and corners) for a seed: made once and shared by every glass of that seed, as
 * are its atlases (below). Why: a heading's glass is made each time its view shows; cutting and rendering the pieces
 * anew each time was most of what showing it cost.
 */
interface Cut {
  readonly shard: readonly Path2D[];
  readonly cx: Float32Array;
  readonly cy: Float32Array;
  /** The largest reach of a corner from its piece's middle (units): a sprite's size. */
  readonly reach: number;
  /** Atlases by px a unit and the stroke's colors (the last few kept). */
  readonly atlases: Map<string, { canvas: HTMLCanvasElement; side: number }>;
}

const cuts = new Map<number, Cut>();

function cutFor(seed: number): Cut {
  let cut = cuts.get(seed);
  if (cut) return cut;
  const random = generator(seed);
  const x0 = -1;
  const y0 = -1;
  const w = WordmarkGeometry.width;
  const h = WordmarkGeometry.height;
  const px: number[] = [];
  const py: number[] = [];
  for (let r = 0; r <= Rows; r++) {
    for (let c = 0; c <= Cols; c++) {
      const inner = r > 0 && r < Rows && c > 0 && c < Cols;
      px.push(x0 + (c / Cols) * w + (inner ? (random() - 0.5) * 0.7 * (w / Cols) : 0));
      py.push(y0 + (r / Rows) * h + (inner ? (random() - 0.5) * 0.7 * (h / Rows) : 0));
    }
  }
  const shard: Path2D[] = [];
  const cx = new Float32Array(Shards);
  const cy = new Float32Array(Shards);
  let reach = 0;
  let k = 0;
  const at = (c: number, r: number) => r * (Cols + 1) + c;
  const add = (a: number, b: number, c: number) => {
    const mx = (px[a] + px[b] + px[c]) / 3;
    const my = (py[a] + py[b] + py[c]) / 3;
    const p = new Path2D();
    p.moveTo(px[a] - mx, py[a] - my);
    p.lineTo(px[b] - mx, py[b] - my);
    p.lineTo(px[c] - mx, py[c] - my);
    p.closePath();
    shard.push(p);
    for (const i of [a, b, c]) reach = Math.max(reach, Math.abs(px[i] - mx), Math.abs(py[i] - my));
    cx[k] = mx;
    cy[k] = my;
    k++;
  };
  for (let r = 0; r < Rows; r++) {
    for (let c = 0; c < Cols; c++) {
      const [a, b, d, e] = [at(c, r), at(c + 1, r), at(c, r + 1), at(c + 1, r + 1)];
      // Each cell split along one diagonal or the other.
      if (random() < 0.5) {
        add(a, b, e);
        add(a, e, d);
      } else {
        add(a, b, d);
        add(b, e, d);
      }
    }
  }
  cut = { shard, cx, cy, reach, atlases: new Map() };
  cuts.set(seed, cut);
  return cut;
}

/** The pieces' look: the logo's four colors along the word, and the edges' color. */
export interface GlassColors {
  readonly logo: readonly [string, string, string, string];
  readonly edge: string;
}

export class Glass {
  /** The generator's state (mulberry32), in a typed array: `#random` is a method a step's loop inlines (no boxing). */
  readonly #seed = new Int32Array(1);
  readonly #cut: Cut;
  readonly #travel: number;
  // Where each piece is (offset from home, rotation), how it moves apart, where an assembly started from, its stagger.
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
  // The pointer's lean on each piece (units, radians) and its speed, on a spring.
  readonly #hx = new Float32Array(Shards);
  readonly #hy = new Float32Array(Shards);
  readonly #hr = new Float32Array(Shards);
  readonly #hvx = new Float32Array(Shards);
  readonly #hvy = new Float32Array(Shards);
  readonly #hvr = new Float32Array(Shards);
  #assembly: GlassTiming = GlassGather;
  #phase: GlassPhase = 'whole';
  // The clock and the pointer as numbers in one typed array (time, since, pointer x, y): a step writes them without
  // boxing (fields holding fractions can allocate on each write).
  readonly #t = new Float64Array(4);
  /** An assembly that fades the pieces in as they start (one from nothing; not a snap back). */
  #fadeIn = false;
  #pointer = false;
  /** Whole and nothing moving: nothing needs a frame (the word is drawn whole unless a piece leans, `#leaned`). */
  #still = true;
  // The stroke made for the last colors (a gradient in drawing units), and those colors' key.
  #stroke: CanvasGradient | string = '#000';
  #strokeKey = '';
  #strokeColors: GlassColors | undefined;
  // The atlas in use and what it was made for: looked up again only when the scale or the colors change (no key string
  // a frame).
  #atlas: { canvas: HTMLCanvasElement; side: number } | undefined;
  #atlasU = 0;
  #atlasKey = '';

  /** The pieces there are. */
  readonly count = Shards;

  /** `travel`: how far the pieces fly (a factor: 1 for a splash's, less for a heading's). */
  constructor(seed = 11, travel = 1) {
    this.#seed[0] = (seed * 7919 + 1) | 0 || 1;
    this.#cut = cutFor(seed);
    this.#travel = travel;
  }

  #random(): number {
    const a = (this.#seed[0] + 0x6d2b79f5) | 0;
    this.#seed[0] = a;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  /** Whole and still: nothing to draw again until something happens (no frames while nothing moves). */
  get resting(): boolean {
    return this.#still;
  }

  get phase(): GlassPhase {
    return this.#phase;
  }

  /** Scattered around the word, then gathering (over `timing`: the splash's, by default). */
  scatter(timing: GlassTiming = GlassGather): void {
    for (let i = 0; i < Shards; i++) {
      const a = this.#random() * Math.PI * 2;
      const d = (18 + this.#random() * 36) * this.#travel;
      this.#fx[i] = Math.cos(a) * d;
      this.#fy[i] = Math.sin(a) * d * 0.7;
      this.#fr[i] = (this.#random() * 2 - 1) * 1.6;
      this.#delay[i] = this.#random() * timing.stagger;
      this.#ox[i] = this.#fx[i];
      this.#oy[i] = this.#fy[i];
      this.#rot[i] = this.#fr[i];
    }
    this.#assembly = timing;
    this.#fadeIn = true;
    this.#begin('assembling');
  }

  /** The pointer over the glass at (x, y) (drawing units). */
  setPointer(x: number, y: number): void {
    this.#pointer = true;
    this.#t[2] = x;
    this.#t[3] = y;
  }

  clearPointer(): void {
    this.#pointer = false;
  }

  /** Pressed at (x, y): the pieces fly apart from there, faster the nearer, spinning; again while apart, further. */
  shatter(x: number, y: number): void {
    const keep = this.#phase === 'broken' ? 1 : 0;
    const { cx, cy } = this.#cut;
    for (let i = 0; i < Shards; i++) {
      const dx = cx[i] + this.#ox[i] - x;
      const dy = cy[i] + this.#oy[i] - y;
      const d = Math.sqrt(dx * dx + dy * dy) || 1;
      const v = ((26 + this.#random() * 34) / (1 + d * 0.04)) * this.#travel;
      this.#vx[i] = this.#vx[i] * keep + (dx / d) * v + (this.#random() - 0.5) * 6;
      this.#vy[i] = this.#vy[i] * keep + (dy / d) * v + (this.#random() - 0.5) * 6;
      this.#vr[i] = this.#vr[i] * keep + (this.#random() * 2 - 1) * 6;
    }
    this.#begin('broken');
  }

  #begin(phase: GlassPhase): void {
    this.#phase = phase;
    this.#t[1] = this.#t[0];
    this.#still = false;
  }

  /** Back together from where the pieces are (they are all there to see: none fades in). */
  #gather(): void {
    this.#fx.set(this.#ox);
    this.#fy.set(this.#oy);
    this.#fr.set(this.#rot);
    for (let i = 0; i < Shards; i++) this.#delay[i] = this.#random() * SnapBack.stagger;
    this.#assembly = SnapBack;
    this.#fadeIn = false;
    this.#begin('assembling');
  }

  /** Advances `dt` seconds. */
  step(dt: number): void {
    this.#t[0] += dt;
    const since = this.#t[0] - this.#t[1];
    const leaning = this.#lean(dt);
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
      const smooth = this.#assembly !== SnapBack;
      for (let i = 0; i < Shards; i++) {
        const p = Math.min(1, Math.max(0, (since - this.#delay[i]) / this.#assembly.duration));
        if (p < 1) done = false;
        const q = 1 - p;
        const k = smooth ? q * q * q * q * q : 1 - snap(p);
        this.#ox[i] = this.#fx[i] * k;
        this.#oy[i] = this.#fy[i] * k;
        this.#rot[i] = this.#fr[i] * k;
      }
      if (done) {
        this.#ox.fill(0);
        this.#oy.fill(0);
        this.#rot.fill(0);
        this.#phase = 'whole';
      }
    }
    this.#still = this.#phase === 'whole' && !leaning;
  }

  /** The pointer's lean, on a spring: whether any piece still leans or moves. */
  #lean(dt: number): boolean {
    const h = Math.min(dt, 1 / 30);
    const { cx, cy } = this.#cut;
    let moving = false;
    let leaned = false;
    for (let i = 0; i < Shards; i++) {
      let tx = 0;
      let ty = 0;
      let tr = 0;
      if (this.#pointer) {
        const dx = cx[i] + this.#ox[i] - this.#t[2];
        const dy = cy[i] + this.#oy[i] - this.#t[3];
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d < HoverReach) {
          const g = 1 - d / HoverReach;
          const f = g * g;
          const inv = 1 / (d || 1);
          tx = dx * inv * HoverLean * f;
          ty = dy * inv * HoverLean * f;
          tr = (dx >= 0 ? 1 : -1) * HoverTilt * f;
        }
      }
      this.#hvx[i] += (Stiffness * (tx - this.#hx[i]) - Damping * this.#hvx[i]) * h;
      this.#hvy[i] += (Stiffness * (ty - this.#hy[i]) - Damping * this.#hvy[i]) * h;
      this.#hvr[i] += (Stiffness * (tr - this.#hr[i]) - Damping * this.#hvr[i]) * h;
      this.#hx[i] += this.#hvx[i] * h;
      this.#hy[i] += this.#hvy[i] * h;
      this.#hr[i] += this.#hvr[i] * h;
      // A frame is needed while a piece moves or is not yet where the pointer holds it; a piece held leaned by a resting
      // pointer needs none (asking for frames then would be a loop while nothing changes).
      if (
        Math.abs(this.#hvx[i]) + Math.abs(this.#hvy[i]) > 0.02 ||
        Math.abs(this.#hvr[i]) > 0.002 ||
        Math.abs(tx - this.#hx[i]) + Math.abs(ty - this.#hy[i]) > 0.01 ||
        Math.abs(tr - this.#hr[i]) > 0.001
      ) {
        moving = true;
      }
      if (Math.abs(this.#hx[i]) + Math.abs(this.#hy[i]) > 0.01 || Math.abs(this.#hr[i]) > 0.001) leaned = true;
    }
    this.#leaned = leaned;
    return moving;
  }

  /** Some piece leans (under the pointer, or settling back): drawn in pieces even when nothing moves. */
  #leaned = false;

  /** The stroke for these colors: the logo's gradient along the word, made again only when they change. */
  #strokeFor(ctx: CanvasRenderingContext2D, colors: GlassColors): CanvasGradient | string {
    // The same colors object: nothing to compare (a frame allocates nothing).
    if (colors === this.#strokeColors) return this.#stroke;
    this.#strokeColors = colors;
    const key = colors.logo.join(',');
    if (key === this.#strokeKey) return this.#stroke;
    const g = ctx.createLinearGradient(0, 0, WordmarkGeometry.width - 2, 0);
    g.addColorStop(0, colors.logo[0]);
    g.addColorStop(0.35, colors.logo[1]);
    g.addColorStop(0.6, colors.logo[2]);
    g.addColorStop(1, colors.logo[3]);
    this.#stroke = g;
    this.#strokeKey = key;
    return g;
  }

  /**
   * Draws the glass: the word whole in one go while at rest, else piece by piece (each piece's part of the word from the
   * shared atlas, one `drawImage` a piece; its edges glinting while it moves). The context is in CSS pixels times `dpr`,
   * the word's origin at (x, y), `scale` px a unit.
   */
  draw(ctx: CanvasRenderingContext2D, x: number, y: number, scale: number, colors: GlassColors, dpr = 1): void {
    const stroke = this.#strokeFor(ctx, colors);
    const paths = word();
    if (this.#still && !this.#leaned) {
      ctx.save();
      ctx.translate(x, y);
      ctx.scale(scale, scale);
      ctx.strokeStyle = stroke;
      ctx.lineWidth = WordmarkGeometry.stroke(scale);
      ctx.lineJoin = 'miter';
      for (let w = 0; w < paths.length; w++) ctx.stroke(paths[w]);
      ctx.restore();
      return;
    }
    const { cx, cy, shard } = this.#cut;
    const u = scale * dpr;
    // Looked up when the scale or the colors change, not when there is none (no canvas to make one on: not retried).
    if (this.#atlasU !== u || this.#atlasKey !== this.#strokeKey) {
      this.#atlas = atlasFor(this.#cut, scale, dpr, stroke, this.#strokeKey);
      this.#atlasU = u;
      this.#atlasKey = this.#strokeKey;
    }
    const atlas = this.#atlas;
    const since = this.#t[0] - this.#t[1];
    ctx.save();
    ctx.lineWidth = 0.35;
    ctx.strokeStyle = colors.edge;
    for (let i = 0; i < Shards; i++) {
      const ox = this.#ox[i] + this.#hx[i];
      const oy = this.#oy[i] + this.#hy[i];
      const rot = this.#rot[i] + this.#hr[i];
      const shown =
        this.#phase === 'assembling' && this.#fadeIn
          ? Math.min(1, Math.max(0, (since - this.#delay[i]) / (this.#assembly.duration * 0.35)))
          : 1;
      if (shown <= 0) continue;
      const c = Math.cos(rot) * u;
      const sn = Math.sin(rot) * u;
      ctx.setTransform(c, sn, -sn, c, (x + (cx[i] + ox) * scale) * dpr, (y + (cy[i] + oy) * scale) * dpr);
      const off = Math.abs(this.#ox[i]) + Math.abs(this.#oy[i]) + Math.abs(this.#rot[i]) * 6;
      const lean = Math.abs(this.#hx[i]) + Math.abs(this.#hy[i]) + Math.abs(this.#hr[i]) * 6;
      const glint = Math.min(1, off / 8 + lean / 3) * 0.22 * shown;
      if (glint > 0.01) {
        ctx.globalAlpha = glint;
        ctx.stroke(shard[i]);
      }
      ctx.globalAlpha = shown;
      if (atlas) {
        const side = atlas.side;
        ctx.drawImage(
          atlas.canvas,
          (i % AtlasCols) * side,
          ((i / AtlasCols) | 0) * side,
          side,
          side,
          -side / (2 * u),
          -side / (2 * u),
          side / u,
          side / u,
        );
      } else {
        ctx.save();
        ctx.clip(shard[i]);
        ctx.translate(-cx[i], -cy[i]);
        ctx.strokeStyle = stroke;
        ctx.lineWidth = WordmarkGeometry.stroke(scale);
        ctx.lineJoin = 'miter';
        for (let w = 0; w < paths.length; w++) ctx.stroke(paths[w]);
        ctx.restore();
      }
    }
    ctx.restore();
    ctx.globalAlpha = 1;
  }
}

/**
 * The pieces' atlas for a cut at this scale, pixel ratio and colors: each piece's part of the word, its middle in the
 * middle of its sprite, rendered once and shared (a heading shown again, another of the same size, reuses it); the last
 * four kept. None where there is no canvas.
 */
function atlasFor(
  cut: Cut,
  scale: number,
  dpr: number,
  stroke: CanvasGradient | string,
  colors: string,
): { canvas: HTMLCanvasElement; side: number } | undefined {
  const u = scale * dpr;
  const key = `${u}|${colors}`;
  const kept = cut.atlases.get(key);
  if (kept) return kept;
  if (typeof document === 'undefined') return undefined;
  const canvas = document.createElement('canvas');
  const g = canvas.getContext('2d');
  if (!g) return undefined;
  const pad = WordmarkGeometry.stroke(scale) / 2 + 0.5;
  const side = Math.ceil(2 * (cut.reach + pad) * u);
  canvas.width = side * AtlasCols;
  canvas.height = side * Math.ceil(Shards / AtlasCols);
  g.strokeStyle = stroke;
  g.lineWidth = WordmarkGeometry.stroke(scale);
  g.lineJoin = 'miter';
  const paths = word();
  for (let i = 0; i < Shards; i++) {
    g.setTransform(u, 0, 0, u, (i % AtlasCols) * side + side / 2, ((i / AtlasCols) | 0) * side + side / 2);
    g.save();
    g.clip(cut.shard[i]);
    g.translate(-cut.cx[i], -cut.cy[i]);
    for (let w = 0; w < paths.length; w++) g.stroke(paths[w]);
    g.restore();
  }
  const atlas = { canvas, side };
  if (cut.atlases.size >= 4) cut.atlases.delete(cut.atlases.keys().next().value!);
  cut.atlases.set(key, atlas);
  return atlas;
}
