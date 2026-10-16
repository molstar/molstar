/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { WordmarkGeometry } from '../../components/icon/wordmark.js';

// The splash's wordmark as glass (`prism-field.ts` draws it): the word cut into triangular shards, which drift in from
// around it and settle together when the splash starts, and fly apart when it is pressed, spinning and slowing, then snap
// back; the pointer over it makes the pieces near it lean away and tilt a little, wobbling back on a spring. The
// shards are drawn as the wordmark is (the same path data, line and the logo's gradient), each its piece of it, its
// edges glinting while it moves; the canvas draws the word throughout (the page's wordmark only places it). In drawing units
// (the wordmark's viewBox: 1 unit a pixel at 1×); everything per frame is preallocated.

/** Columns and rows of the grid the word is cut on (each cell two triangles, the inner corners moved). */
const Cols = 11;
const Rows = 4;
const Shards = Cols * Rows * 2;
/**
 * Two ways to come together. The first assembly is smooth: pieces drifting in from around the word and settling (eased
 * out, no overshoot), over a while, their starts spread. After a press they snap back: quick, a little past home.
 * Seconds a shard takes, and how far apart their starts are.
 */
export interface GlassTiming {
  readonly duration: number;
  readonly stagger: number;
}
const Gather: GlassTiming = { duration: 2.4, stagger: 0.7 };
const SnapBack: GlassTiming = { duration: 0.45, stagger: 0.12 };
/** Seconds the pieces stay apart before they come back. */
const Hold = 0.8;
/** How fast the pieces slow when apart (per second). */
const Drag = 3;
/**
 * The pointer over the glass: how far it reaches (drawing units), how far a piece leans away from it at most (units) and
 * how much it tilts (radians), and the spring the pieces wobble back on (stiffness, damping: a little bouncy).
 */
const HoverReach = 14;
const HoverLean = 1.4;
const HoverTilt = 0.1;
const Stiffness = 140;
const Damping = 11;

/** Sprites a row in the pieces' atlas. */
const AtlasCols = 11;

/** The mean distance of the pieces from home (drawing units) at which `assembled` is 1/e. */
const AssembledScale = 3;

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
  // The pointer's lean on each piece (units, radians) and its speed, on a spring.
  readonly #hx = new Float32Array(Shards);
  readonly #hy = new Float32Array(Shards);
  readonly #hr = new Float32Array(Shards);
  readonly #hvx = new Float32Array(Shards);
  readonly #hvy = new Float32Array(Shards);
  readonly #hvr = new Float32Array(Shards);
  /** The assembly under way: smooth (the first) or a snap back. */
  #assembly: GlassTiming = Gather;
  #phase: GlassPhase = 'whole';
  #since = 0;
  #time = 0;
  /** An assembly that fades the pieces in as they start (the first one, from nothing; not a reassembly). */
  #fadeIn = false;
  /** The pointer over the glass (drawing units), if it is. */
  #pointer = false;
  #pu = 0;
  #pv = 0;
  /** Whether any piece is away from home or leaning (otherwise the word is drawn whole, in one go). */
  #still = true;
  /** Told when the glass becomes whole or stops being whole (the star splits the light only while whole). */
  onWhole: (whole: boolean) => void = () => {};
  /** Told when it snaps back together after a break (the light pops). */
  onSnap: () => void = () => {};

  /** How far the pieces travel (a factor: 1 for the empty canvas's, less for a small wordmark). */
  readonly #travel: number;

  /** The pieces there are. */
  readonly count = Shards;
  // The pieces' atlas (`#atlasFor`): its canvas, the px a unit and the stroke it was made for, a sprite's size.
  #atlas: HTMLCanvasElement | undefined;
  #atlasU = 0;
  #atlasStroke: CanvasGradient | string | undefined;
  #spriteW = 0;
  #spriteH = 0;
  /** Each piece's corners about its middle (drawing units; six a piece: x and y of each corner). */
  readonly corners = new Float32Array(Shards * 6);

  constructor(seed = 11, travel = 1) {
    this.#random = generator(seed);
    this.#travel = travel;
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
      this.corners.set([px[a] - cx, py[a] - cy, px[b] - cx, py[b] - cy, px[c] - cx, py[c] - cy], k * 6);
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

  /**
   * Writes piece `i`'s place now into `out` (five values): its middle (drawing units, the lean under the pointer
   * included), its rotation (radians), and its middle at home; for a host that follows the pieces (the quantum
   * splash's potential).
   */
  place(i: number, out: Float32Array): void {
    out[0] = this.#cx[i] + this.#ox[i] + this.#hx[i];
    out[1] = this.#cy[i] + this.#oy[i] + this.#hy[i];
    out[2] = this.#rot[i] + this.#hr[i];
    out[3] = this.#cx[i];
    out[4] = this.#cy[i];
  }

  /**
   * How assembled it is (0..1): 1 with every piece home, falling as they are away (their mean distance in drawing
   * units, a turn counted as distance), toward 0 when scattered; for a host whose light follows the glass.
   */
  get assembled(): number {
    let away = 0;
    for (let i = 0; i < Shards; i++) {
      const ox = this.#ox[i];
      const oy = this.#oy[i];
      away += Math.sqrt(ox * ox + oy * oy) + Math.abs(this.#rot[i]) * 6;
    }
    return Math.exp(-away / Shards / AssembledScale);
  }

  /**
   * How hard the glass pulls the light in as it snaps back (0..1): strongest as the pieces start home, gone as they land
   * (the splash's gravity; 0 otherwise).
   */
  get pull(): number {
    if (this.#phase !== 'assembling' || this.#assembly !== SnapBack) return 0;
    const span = SnapBack.duration + SnapBack.stagger;
    const p = Math.min(1, (this.#time - this.#since) / span);
    return (1 - p) * (1 - p);
  }

  /** Whole and still: nothing to draw again until something happens (no frames while nothing moves). */
  get resting(): boolean {
    return this.#still;
  }

  get phase(): GlassPhase {
    return this.#phase;
  }

  /** Whole: the star splits the light. */
  get whole(): boolean {
    return this.#phase === 'whole';
  }

  /** Scattered around the word, then gathering smoothly (over `timing`: the splash's, by default). */
  scatter(timing: GlassTiming = Gather): void {
    for (let i = 0; i < Shards; i++) {
      const a = this.#random() * Math.PI * 2;
      const d = (18 + this.#random() * 36) * this.#travel;
      this.#fx[i] = Math.cos(a) * d;
      this.#fy[i] = Math.sin(a) * d * 0.7;
      this.#fr[i] = (this.#random() * 2 - 1) * 1.6;
      this.#delay[i] = this.#random() * timing.stagger;
    }
    this.#assembly = timing;
    this.#fadeIn = true;
    this.#begin('assembling');
  }

  /** The pointer over the glass at (x, y) (drawing units); none: it left. */
  setPointer(x: number, y: number): void {
    this.#pointer = true;
    this.#pu = x;
    this.#pv = y;
  }

  clearPointer(): void {
    this.#pointer = false;
  }

  /** Pressed at (x, y) (drawing units): the pieces fly apart from there, faster the nearer, spinning. */
  shatter(x: number, y: number): void {
    // Again while in pieces (or snapping back): from where they are, the kick added to how they move, the wait anew.
    const moving = this.#phase === 'broken';
    for (let i = 0; i < Shards; i++) {
      const dx = this.#cx[i] + this.#ox[i] - x;
      const dy = this.#cy[i] + this.#oy[i] - y;
      const d = Math.hypot(dx, dy) || 1;
      const v = ((26 + this.#random() * 34) / (1 + d * 0.04)) * this.#travel;
      const keep = moving ? 1 : 0;
      this.#vx[i] = this.#vx[i] * keep + (dx / d) * v + (this.#random() - 0.5) * 6;
      this.#vy[i] = this.#vy[i] * keep + (dy / d) * v + (this.#random() - 0.5) * 6;
      this.#vr[i] = this.#vr[i] * keep + (this.#random() * 2 - 1) * 6;
    }
    this.#begin('broken');
  }

  /**
   * Whether (x, y) (drawing units) is on the glass: whole, the word's box (a little around); in pieces, near any of them
   * (a press there shatters them further).
   */
  hits(x: number, y: number): boolean {
    if (this.#phase === 'whole') {
      return x >= -2 && x <= WordmarkGeometry.width + 1 && y >= -2 && y <= WordmarkGeometry.height + 1;
    }
    for (let i = 0; i < Shards; i++) {
      const dx = this.#cx[i] + this.#ox[i] - x;
      const dy = this.#cy[i] + this.#oy[i] - y;
      if (dx * dx + dy * dy < 36) return true;
    }
    return false;
  }

  #begin(phase: GlassPhase): void {
    const wasWhole = this.#phase === 'whole';
    this.#phase = phase;
    this.#since = this.#time;
    this.#still = false;
    if (wasWhole && phase !== 'whole') this.onWhole(false);
  }

  /** Back together from where the pieces are now. */
  #gather(): void {
    for (let i = 0; i < Shards; i++) {
      this.#fx[i] = this.#ox[i];
      this.#fy[i] = this.#oy[i];
      this.#fr[i] = this.#rot[i];
      this.#delay[i] = this.#random() * SnapBack.stagger;
    }
    this.#assembly = SnapBack;
    // From where the pieces are: they are all there to see, none fades in.
    this.#fadeIn = false;
    this.#begin('assembling');
  }

  step(dt: number): void {
    this.#time += dt;
    const since = this.#time - this.#since;
    this.#lean(dt);
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
        const p = Math.min(1, Math.max(0, (since - this.#delay[i]) / this.#assembly.duration));
        if (p < 1) done = false;
        // Smooth: eased out (quintic), settling; a snap back: past home a little.
        const q = 1 - p;
        const k = this.#assembly === SnapBack ? 1 - snap(p) : q * q * q * q * q;
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
   * The pointer's lean: each piece near it leans away and tilts a little (more the nearer), on a spring, so it wobbles
   * back when the pointer moves on.
   */
  #lean(dt: number): void {
    const h = Math.min(dt, 1 / 30);
    let moving = false;
    for (let i = 0; i < Shards; i++) {
      let tx = 0;
      let ty = 0;
      let tr = 0;
      if (this.#pointer) {
        const dx = this.#cx[i] + this.#ox[i] - this.#pu;
        const dy = this.#cy[i] + this.#oy[i] - this.#pv;
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d < HoverReach) {
          const g = 1 - d / HoverReach;
          const f = g * g;
          const inv = 1 / (d || 1);
          tx = dx * inv * HoverLean * f;
          ty = dy * inv * HoverLean * f;
          // Turned away: a piece left of the pointer one way, right of it the other.
          tr = (dx >= 0 ? 1 : -1) * HoverTilt * f;
        }
      }
      this.#hvx[i] += (Stiffness * (tx - this.#hx[i]) - Damping * this.#hvx[i]) * h;
      this.#hvy[i] += (Stiffness * (ty - this.#hy[i]) - Damping * this.#hvy[i]) * h;
      this.#hvr[i] += (Stiffness * (tr - this.#hr[i]) - Damping * this.#hvr[i]) * h;
      this.#hx[i] += this.#hvx[i] * h;
      this.#hy[i] += this.#hvy[i] * h;
      this.#hr[i] += this.#hvr[i] * h;
      if (
        Math.abs(this.#hx[i]) + Math.abs(this.#hy[i]) > 0.01 ||
        Math.abs(this.#hr[i]) > 0.001 ||
        Math.abs(this.#hvx[i]) + Math.abs(this.#hvy[i]) > 0.02
      ) {
        moving = true;
      }
    }
    this.#still = this.#phase === 'whole' && !moving;
  }

  /**
   * Draws the glass: the word whole in one go while at rest, else piece by piece: the word's line in `stroke` (the
   * logo's gradient, made in drawing units), each piece its part of it, and the piece's edges in `edge`, glinting while
   * it moves. The context is in canvas CSS pixels (times `dpr`); the word's drawing origin is at (x, y), `scale` px a
   * unit. The pieces come from an atlas, each drawn once at this scale, pixel ratio and stroke (one `drawImage` a
   * piece a frame, not a clip and the word's strokes); without a canvas to make it on, each is clipped and stroked.
   */
  draw(
    ctx: CanvasRenderingContext2D,
    x: number,
    y: number,
    scale: number,
    stroke: CanvasGradient | string,
    edge: string,
    alpha: number,
    dpr = 1,
  ): void {
    const since = this.#time - this.#since;
    if (this.#still) {
      // Whole and at rest: the word in one go, as the wordmark draws it.
      ctx.save();
      ctx.translate(x, y);
      ctx.scale(scale, scale);
      ctx.globalAlpha = alpha;
      ctx.strokeStyle = stroke;
      ctx.lineWidth = WordmarkGeometry.stroke(scale);
      ctx.lineJoin = 'miter';
      for (let w = 0; w < this.#word.length; w++) ctx.stroke(this.#word[w]);
      ctx.restore();
      ctx.globalAlpha = 1;
      return;
    }
    const atlas = this.#atlasFor(scale, dpr, stroke);
    const u = scale * dpr;
    ctx.save();
    ctx.lineWidth = 0.35;
    ctx.strokeStyle = edge;
    for (let i = 0; i < Shards; i++) {
      const ox = this.#ox[i] + this.#hx[i];
      const oy = this.#oy[i] + this.#hy[i];
      const rot = this.#rot[i] + this.#hr[i];
      // A shard flying in (the first assembly) fades in as it starts.
      const shown =
        this.#phase === 'assembling' && this.#fadeIn
          ? Math.min(1, Math.max(0, (since - this.#delay[i]) / (this.#assembly.duration * 0.35)))
          : 1;
      if (shown <= 0) continue;
      // The piece's frame in device pixels: its middle where it is now, turned, a drawing unit `u` px.
      const c = Math.cos(rot) * u;
      const sn = Math.sin(rot) * u;
      ctx.setTransform(c, sn, -sn, c, (x + (this.#cx[i] + ox) * scale) * dpr, (y + (this.#cy[i] + oy) * scale) * dpr);
      // Its edges: a glint while it is away from home.
      const off = Math.abs(this.#ox[i]) + Math.abs(this.#oy[i]) + Math.abs(this.#rot[i]) * 6;
      const lean = Math.abs(this.#hx[i]) + Math.abs(this.#hy[i]) + Math.abs(this.#hr[i]) * 6;
      const glint = Math.min(1, off / 8 + lean / 3) * 0.22 * alpha * shown;
      if (glint > 0.01) {
        ctx.globalAlpha = glint;
        ctx.stroke(this.#shard[i]);
      }
      ctx.globalAlpha = alpha * shown;
      if (atlas) {
        const w = this.#spriteW;
        const h = this.#spriteH;
        ctx.drawImage(
          atlas,
          (i % AtlasCols) * w,
          ((i / AtlasCols) | 0) * h,
          w,
          h,
          -w / (2 * u),
          -h / (2 * u),
          w / u,
          h / u,
        );
      } else {
        ctx.save();
        ctx.clip(this.#shard[i]);
        ctx.translate(-this.#cx[i], -this.#cy[i]);
        ctx.strokeStyle = stroke;
        ctx.lineWidth = WordmarkGeometry.stroke(scale);
        ctx.lineJoin = 'miter';
        for (let w = 0; w < this.#word.length; w++) ctx.stroke(this.#word[w]);
        ctx.restore();
      }
    }
    ctx.restore();
    ctx.globalAlpha = 1;
  }

  /**
   * The pieces' atlas for this scale, pixel ratio and stroke: each piece's part of the word, its middle in the middle of
   * its sprite (`AtlasCols` a row), made again only when one of them changes; none where there is no canvas.
   */
  #atlasFor(scale: number, dpr: number, stroke: CanvasGradient | string): HTMLCanvasElement | undefined {
    const u = scale * dpr;
    if (this.#atlas && this.#atlasU === u && this.#atlasStroke === stroke) return this.#atlas;
    if (typeof document === 'undefined') return undefined;
    const canvas = this.#atlas ?? document.createElement('canvas');
    const g = canvas.getContext('2d');
    if (!g) return undefined;
    // A sprite fits the largest piece about its middle, with room for the line.
    let reach = 0;
    for (let k = 0; k < this.corners.length; k++) reach = Math.max(reach, Math.abs(this.corners[k]));
    const pad = WordmarkGeometry.stroke(scale) / 2 + 0.5;
    const side = Math.ceil(2 * (reach + pad) * u);
    this.#spriteW = side;
    this.#spriteH = side;
    canvas.width = side * AtlasCols;
    canvas.height = side * Math.ceil(Shards / AtlasCols);
    g.clearRect(0, 0, canvas.width, canvas.height);
    g.strokeStyle = stroke;
    g.lineWidth = WordmarkGeometry.stroke(scale);
    g.lineJoin = 'miter';
    for (let i = 0; i < Shards; i++) {
      g.setTransform(u, 0, 0, u, (i % AtlasCols) * side + side / 2, ((i / AtlasCols) | 0) * side + side / 2);
      g.save();
      g.clip(this.#shard[i]);
      g.translate(-this.#cx[i], -this.#cy[i]);
      for (let w = 0; w < this.#word.length; w++) g.stroke(this.#word[w]);
      g.restore();
    }
    g.setTransform(1, 0, 0, 1, 0, 0);
    this.#atlas = canvas;
    this.#atlasU = u;
    this.#atlasStroke = stroke;
    return canvas;
  }
}
