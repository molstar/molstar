/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The empty canvas's splash as a quantum system (beside prism-field.ts, the same surface; design/shell.md): the light is
// a wave, ψ on a 128 × 128 grid round the star, moved on by the Schrödinger equation (ħ = m = 1, H = −∇²/2 + V) with the
// split-step Fourier method: half a step of the potential, the kinetic step in momentum space (an FFT there and back),
// the other half. Wave packets come in from the left, each of its own wavenumber (its color); the Mol* star is a
// potential well, glass: a packet's refraction through it depends on its energy, so the star splits the light by
// wavenumber the way a prism does. The specks ride the wave (Bohmian trajectories: each moves with the probability
// current, ∇ of ψ's phase, a midpoint Runge–Kutta step) and take their color from the wavenumber they carry once they
// leave the straight way; under them the wave itself, |ψ|², faintly, its wavefronts in it. The pointer is a barrier the
// wave goes round, a press kicks it outward, a shatter takes the prism away until the glass is back (and its snap back
// is a well that draws the light in for a moment). The grid's edges absorb (no wrapping round).
//
// Allocation-free per frame (memory: render paths); it draws nothing while hidden or stopped, one still frame under
// reduced motion.

import { WordmarkGeometry } from '../../components/icon/wordmark.js';
import { Glass } from './glass.js';

export interface QuantumFieldOptions {
  /** Wave packets that come in per second. */
  rate: number;
  /** Their mean wavenumber (radians a cell): the light's mean color. */
  k: number;
  /** How far their wavenumbers spread round it (0..1 of `k`): the spectrum the prism splits. */
  spectrum: number;
  /** Their speed at the mean wavenumber (px/s): how fast time runs. */
  speed: number;
  /** The prism: the star's well's depth, in a mean packet's energies. */
  prism: number;
  /** The wave itself shows (|ψ|², 0: only the specks). */
  glow: number;
  /** The wavefronts show in it (0: the density alone). */
  ripples: number;
  /** Specks riding the wave. */
  specks: number;
  /** The most a speck shows (0..1). */
  alpha: number;
  /** A speck's radius (px). */
  size: number;
  /** The rays in the logo's colors; in the stream's color when off. */
  colors: boolean;
  /** It answers the pointer. */
  interactive: boolean;
  /** The pointer is a barrier: its radius (px) and height (in a mean packet's energies). */
  radius: number;
  barrier: number;
  /** How hard a press kicks the wave outward (radians a cell). */
  kick: number;
  seed: number;
}

export const QuantumFieldDefaults: QuantumFieldOptions = {
  rate: 2.4,
  k: 1.1,
  spectrum: 0.4,
  speed: 64,
  prism: 1.2,
  glow: 0.35,
  ripples: 0.45,
  specks: 650,
  alpha: 0.6,
  size: 1.3,
  colors: true,
  interactive: true,
  radius: 36,
  barrier: 2,
  kick: 1.4,
  seed: 3,
};

/** The grid: cells a side (a power of two, for the FFT), and px a cell. */
const N = 128;
const Bits = 7;
const Cell = 5;
const Cells = N * N;
/** The time step of the evolution (in the grid's units), and the most steps a frame. */
const Dt = 0.1;
const MaxSteps = 4;
/** The absorbing edge: cells wide, and how much it takes a step at the very edge. */
const Edge = 14;
const Absorb = 0.12;
/** The prism's radius (px), its arms' share of it, and how soft its edge is (cells). */
const PrismRadius = 34;
const PrismArms = 0.45;
const PrismSoft = 1.5;
/** A packet's half size (cells, along and across its way), and where the light comes in (cells from the left edge). */
const PacketAlong = 8;
const PacketAcross = 6;
const Source = Edge + 6;
/** The snap back's well: how far it reaches (px) and how deep it is (in a mean packet's energies, at full pull). */
const GravityRadius = 320;
const GravityDepth = 1.6;
/** Specks at most; how long one lives (s); below this share of the densest, a speck lets go. */
const Capacity = 1500;
const LifeMin = 2.2;
const LifeMax = 5;
const Thin = 0.004;
/** A speck going within this of straight on (radians) is light not yet split: the stream's color. */
const Straight = 0.22;
/** Alpha levels a frame batches specks into (per color). */
const Levels = 10;
const Colors = 5;
/** Seconds of the start the light comes in faster. */
const Intro = 2;

const ColorVars = [
  '--msp-color-text-subtle',
  '--msp-color-brand-1',
  '--msp-color-brand-2',
  '--msp-color-brand-3',
  '--msp-color-brand-4',
];

function generator(seed: number): () => number {
  let a = seed >>> 0 || 1;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// The FFT's tables: the bit reversal and the twiddles (e^{−2πij/N}).
const reversed = new Uint16Array(N);
for (let i = 0; i < N; i++) {
  let r = 0;
  for (let b = 0; b < Bits; b++) r |= ((i >> b) & 1) << (Bits - 1 - b);
  reversed[i] = r;
}
const twiddleRe = new Float32Array(N / 2);
const twiddleIm = new Float32Array(N / 2);
for (let j = 0; j < N / 2; j++) {
  twiddleRe[j] = Math.cos((-2 * Math.PI * j) / N);
  twiddleIm[j] = Math.sin((-2 * Math.PI * j) / N);
}

/** An in-place radix-2 FFT of N values at `offset`, `stride` apart (inverse: unscaled). */
function fft(re: Float32Array, im: Float32Array, offset: number, stride: number, inverse: boolean): void {
  for (let i = 0; i < N; i++) {
    const j = reversed[i];
    if (j <= i) continue;
    const a = offset + i * stride;
    const b = offset + j * stride;
    let t = re[a];
    re[a] = re[b];
    re[b] = t;
    t = im[a];
    im[a] = im[b];
    im[b] = t;
  }
  const sign = inverse ? -1 : 1;
  for (let size = 2; size <= N; size <<= 1) {
    const half = size >> 1;
    const step = N / size;
    for (let start = 0; start < N; start += size) {
      for (let k = 0; k < half; k++) {
        const wr = twiddleRe[k * step];
        const wi = sign * twiddleIm[k * step];
        const a = offset + (start + k) * stride;
        const b = a + half * stride;
        const xr = re[b] * wr - im[b] * wi;
        const xi = re[b] * wi + im[b] * wr;
        re[b] = re[a] - xr;
        im[b] = im[a] - xi;
        re[a] += xr;
        im[a] += xi;
      }
    }
  }
}

/** The 2D FFT of the N × N grid (rows, then columns; inverse: unscaled). Exported for its tests. */
export function fft2(re: Float32Array, im: Float32Array, inverse: boolean): void {
  for (let y = 0; y < N; y++) fft(re, im, y * N, 1, inverse);
  for (let x = 0; x < N; x++) fft(re, im, x, N, inverse);
}

// The kinetic step in momentum space, e^{−i k² Dt / 2}, with the inverse FFT's 1 / N² folded in; and the absorbing
// edge's share a step (1 inside).
const kineticRe = new Float32Array(Cells);
const kineticIm = new Float32Array(Cells);
const edgeKeep = new Float32Array(Cells);
/** How much of the wave shows (1 inside, easing to 0 across the absorbing edge): the grid's square never shows. */
const edgeShow = new Float32Array(Cells);
for (let y = 0; y < N; y++) {
  const ky = (2 * Math.PI * (y < N / 2 ? y : y - N)) / N;
  const ey = Math.max(0, Edge - Math.min(y, N - 1 - y)) / Edge;
  for (let x = 0; x < N; x++) {
    const kx = (2 * Math.PI * (x < N / 2 ? x : x - N)) / N;
    const phase = (-(kx * kx + ky * ky) * Dt) / 2;
    kineticRe[y * N + x] = Math.cos(phase) / Cells;
    kineticIm[y * N + x] = Math.sin(phase) / Cells;
    const ex = Math.max(0, Edge - Math.min(x, N - 1 - x)) / Edge;
    const e = Math.max(ex, ey);
    edgeKeep[y * N + x] = 1 - Absorb * e * e;
    const inside = 1 - e;
    edgeShow[y * N + x] = inside * inside * (3 - 2 * inside);
  }
}

export class QuantumField {
  #options: QuantumFieldOptions;
  readonly #canvas: HTMLCanvasElement;
  readonly #ctx: CanvasRenderingContext2D | null;
  #random: () => number;
  readonly #still: boolean;

  // The wave, the potential and its half step (with the edge), the current's velocity (cells a time unit) and the density.
  readonly #re = new Float32Array(Cells);
  readonly #im = new Float32Array(Cells);
  readonly #v = new Float32Array(Cells);
  readonly #halfRe = new Float32Array(Cells);
  readonly #halfIm = new Float32Array(Cells);
  readonly #vx = new Float32Array(Cells);
  readonly #vy = new Float32Array(Cells);
  readonly #rho = new Float32Array(Cells);
  /** The density's running sum, to draw where a speck starts. */
  readonly #cdf = new Float32Array(Cells);
  #rhoMax = 1e-6;
  #potential = true;
  #clock = 0;
  #emit = 0;

  // The specks: position (canvas px), age and life (s), color (0: the stream's; 1–4: the logo's), a phase for their size.
  readonly #x = new Float32Array(Capacity);
  readonly #y = new Float32Array(Capacity);
  readonly #age = new Float32Array(Capacity);
  readonly #life = new Float32Array(Capacity);
  readonly #shimmer = new Float32Array(Capacity);
  readonly #color = new Uint8Array(Capacity);
  #alive = 0;
  readonly #bucket = new Uint16Array(Capacity);
  readonly #bucketOf = new Uint16Array(Capacity);
  readonly #counts = new Uint16Array(Colors * Levels);
  readonly #starts = new Uint16Array(Colors * Levels);

  // The wave's picture: a pixel a cell, drawn scaled up (smoothed).
  readonly #image: ImageData | undefined;
  readonly #picture: HTMLCanvasElement | undefined;
  readonly #pictureCtx: CanvasRenderingContext2D | null = null;

  readonly #colors: string[] = ['#888', '#7a4f2f', '#e0ae2c', '#6f9a80', '#3f86ad'];
  readonly #rgb = new Uint8Array(Colors * 3);
  #gradient: CanvasGradient | undefined;
  #width = 0;
  #height = 0;
  #dpr = 1;
  /** The star's middle (canvas px), and the grid's corner. */
  #ox = 0;
  #oy = 0;
  #gx = 0;
  #gy = 0;
  #time = 0;
  #last = 0;
  #frame = 0;
  #running = false;
  #visible = true;
  #colorsRead = Number.NEGATIVE_INFINITY;

  #px = 0;
  #py = 0;
  #pointer = false;

  #glass: Glass | undefined;
  #wordX = 0;
  #wordY = 0;
  #wordScale = 1;
  #wasWhole = true;
  #wasPull = 0;
  #bursts = 0;

  constructor(canvas: HTMLCanvasElement, options: Partial<QuantumFieldOptions> = {}, init: { still?: boolean } = {}) {
    this.#canvas = canvas;
    this.#ctx = canvas.getContext('2d');
    this.#options = { ...QuantumFieldDefaults, ...options };
    this.#random = generator(this.#options.seed);
    this.#still = !!init.still;
    if (typeof document !== 'undefined' && this.#ctx) {
      const picture = document.createElement('canvas');
      picture.width = N;
      picture.height = N;
      const pctx = picture.getContext('2d');
      if (pctx) {
        this.#picture = picture;
        this.#pictureCtx = pctx;
        this.#image = pctx.createImageData(N, N);
      }
    }
  }

  setOptions(options: Partial<QuantumFieldOptions>): void {
    const seed = this.#options.seed;
    this.#options = { ...this.#options, ...options };
    if (this.#options.seed !== seed) {
      this.#random = generator(this.#options.seed);
      this.#re.fill(0);
      this.#im.fill(0);
      this.#alive = 0;
    }
    this.#potential = true;
    if (this.#still) this.#drawStill();
  }

  setOrigin(x: number, y: number): void {
    this.#ox = x;
    this.#oy = y;
    // The star a little left of the grid's middle: more room for the light it sends on.
    this.#gx = x - 0.47 * N * Cell;
    this.#gy = y - 0.5 * N * Cell;
    this.#potential = true;
    if (this.#still) this.#drawStill();
  }

  resize(): void {
    const w = this.#canvas.clientWidth;
    const h = this.#canvas.clientHeight;
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    if (w === this.#width && h === this.#height && dpr === this.#dpr) return;
    this.#width = w;
    this.#height = h;
    this.#dpr = dpr;
    this.#canvas.width = Math.round(w * dpr);
    this.#canvas.height = Math.round(h * dpr);
    this.readColors();
    if (this.#still) this.#drawStill();
  }

  readColors(): void {
    if (typeof getComputedStyle !== 'function') return;
    const el = this.#canvas;
    const before = el.style.color;
    for (let i = 0; i < Colors; i++) {
      el.style.color = `var(${ColorVars[i]})`;
      const v = getComputedStyle(el).color;
      if (v) this.#colors[i] = v;
      const m = /rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)/.exec(v);
      if (m) {
        this.#rgb[i * 3] = +m[1];
        this.#rgb[i * 3 + 1] = +m[2];
        this.#rgb[i * 3 + 2] = +m[3];
      }
    }
    el.style.color = before;
    if (this.#ctx) {
      const g = this.#ctx.createLinearGradient(0, 0, WordmarkGeometry.width - 2, 0);
      g.addColorStop(0, this.#colors[1]);
      g.addColorStop(0.35, this.#colors[2]);
      g.addColorStop(0.6, this.#colors[3]);
      g.addColorStop(1, this.#colors[4]);
      this.#gradient = g;
    }
  }

  setWord(x: number, y: number, scale: number): boolean {
    this.#wordX = x;
    this.#wordY = y;
    this.#wordScale = scale;
    if (this.#glass) return true;
    if (this.#still || typeof Path2D !== 'function') return false;
    const glass = new Glass(this.#options.seed + 8);
    // Back together: the prism is there again, and lets out a little of the light it held.
    glass.onSnap = () => {
      this.#bursts += 4;
    };
    this.#glass = glass;
    glass.scatter();
    return true;
  }

  shatterAt(x: number, y: number): boolean {
    const glass = this.#glass;
    if (!glass || !this.#options.interactive || this.#still) return false;
    const u = (x - this.#wordX) / this.#wordScale;
    const v = (y - this.#wordY) / this.#wordScale;
    if (!glass.hits(u, v)) return false;
    glass.shatter(u, v);
    this.#kickAt(x, y, this.#options.kick * 1.6, 90);
    this.#potential = true;
    return true;
  }

  pointerMove(x: number, y: number, _now: number): void {
    if (!this.#options.interactive || this.#still) return;
    this.#px = x;
    this.#py = y;
    this.#pointer = true;
    this.#potential = true;
    this.#glass?.setPointer((x - this.#wordX) / this.#wordScale, (y - this.#wordY) / this.#wordScale);
  }

  pointerLeave(): void {
    this.#pointer = false;
    this.#potential = true;
    this.#glass?.clearPointer();
  }

  pointerDown(x: number, y: number, now: number): void {
    if (!this.#options.interactive || this.#still) return;
    this.pointerMove(x, y, now);
    this.#kickAt(x, y, this.#options.kick, 70);
  }

  start(): void {
    this.#running = true;
    if (this.#still) this.#drawStill();
    else this.#schedule();
  }

  stop(): void {
    this.#running = false;
    if (this.#frame) cancelAnimationFrame(this.#frame);
    this.#frame = 0;
  }

  setVisible(visible: boolean): void {
    this.#visible = visible;
    if (visible) this.#schedule();
  }

  dispose(): void {
    this.stop();
  }

  #schedule(): void {
    if (this.#still || !this.#running || !this.#visible || this.#frame) return;
    this.#last = 0;
    this.#frame = requestAnimationFrame(this.#tick);
  }

  readonly #tick = (now: number) => {
    this.#frame = 0;
    if (!this.#running || !this.#visible) return;
    const dt = this.#last ? Math.min(0.1, (now - this.#last) / 1000) : 1 / 60;
    this.#last = now;
    if (now - this.#colorsRead > 1500) {
      this.#colorsRead = now;
      this.readColors();
    }
    this.#step(dt);
    this.#draw();
    this.#frame = requestAnimationFrame(this.#tick);
  };

  #drawStill(): void {
    if (!this.#width) return;
    this.#re.fill(0);
    this.#im.fill(0);
    this.#alive = 0;
    this.#time = 0;
    this.#random = generator(this.#options.seed);
    for (let t = 0; t < 6; t += 1 / 20) this.#step(1 / 20);
    this.#draw();
  }

  /** Time units a second: a packet of the mean wavenumber goes `speed` px a second (its group velocity is k). */
  #rate(): number {
    const o = this.#options;
    return o.speed / (Math.max(0.2, o.k) * Cell);
  }

  /** Adds a wave packet at grid (cx, cy) going `angle` with wavenumber `k`. */
  #packet(cx: number, cy: number, angle: number, k: number, amplitude: number): void {
    const re = this.#re;
    const im = this.#im;
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    const kx = k * c;
    const ky = k * s;
    const reach = 3 * Math.max(PacketAlong, PacketAcross);
    const x0 = Math.max(1, Math.floor(cx - reach));
    const x1 = Math.min(N - 2, Math.ceil(cx + reach));
    const y0 = Math.max(1, Math.floor(cy - reach));
    const y1 = Math.min(N - 2, Math.ceil(cy + reach));
    const a2 = 2 * PacketAlong * PacketAlong;
    const c2 = 2 * PacketAcross * PacketAcross;
    for (let y = y0; y <= y1; y++) {
      const dy = y - cy;
      for (let x = x0; x <= x1; x++) {
        const dx = x - cx;
        const along = dx * c + dy * s;
        const across = -dx * s + dy * c;
        const g = amplitude * Math.exp(-(along * along) / a2 - (across * across) / c2);
        if (g < 1e-4) continue;
        const phase = kx * dx + ky * dy;
        const i = y * N + x;
        re[i] += g * Math.cos(phase);
        im[i] += g * Math.sin(phase);
      }
    }
  }

  /** Multiplies the wave round (x, y) (canvas px) by an outward phase: momentum away from the press. */
  #kickAt(x: number, y: number, kick: number, radius: number): void {
    const cx = (x - this.#gx) / Cell;
    const cy = (y - this.#gy) / Cell;
    const r = radius / Cell;
    const x0 = Math.max(0, Math.floor(cx - r));
    const x1 = Math.min(N - 1, Math.ceil(cx + r));
    const y0 = Math.max(0, Math.floor(cy - r));
    const y1 = Math.min(N - 1, Math.ceil(cy + r));
    for (let gy = y0; gy <= y1; gy++) {
      for (let gx = x0; gx <= x1; gx++) {
        const dx = gx - cx;
        const dy = gy - cy;
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d >= r) continue;
        const f = 1 - d / r;
        const phase = kick * d * f * f;
        const i = gy * N + gx;
        const c = Math.cos(phase);
        const s = Math.sin(phase);
        const a = this.#re[i];
        const b = this.#im[i];
        this.#re[i] = a * c - b * s;
        this.#im[i] = a * s + b * c;
      }
    }
  }

  /** The potential (the prism while the glass is whole, the pointer's barrier, the snap back's well) and its half step. */
  #buildPotential(pull: number): void {
    const o = this.#options;
    const v = this.#v;
    v.fill(0);
    const e0 = (o.k * o.k) / 2;
    const sx = (this.#ox - this.#gx) / Cell;
    const sy = (this.#oy - this.#gy) / Cell;
    if (!this.#glass || this.#glass.whole) {
      // The star: a well with eight arms.
      const r = PrismRadius / Cell;
      const depth = o.prism * e0;
      const reach = Math.ceil(r + PrismSoft + 1);
      for (let y = Math.max(0, Math.floor(sy - reach)); y <= Math.min(N - 1, sy + reach); y++) {
        for (let x = Math.max(0, Math.floor(sx - reach)); x <= Math.min(N - 1, sx + reach); x++) {
          const dx = x - sx;
          const dy = y - sy;
          const d = Math.sqrt(dx * dx + dy * dy);
          const c4 = Math.cos(4 * Math.atan2(dy, dx));
          const edge = r * (1 - PrismArms + PrismArms * c4 * c4);
          const t = Math.min(1, Math.max(0, (edge - d) / PrismSoft + 0.5));
          if (t > 0) v[y * N + x] -= depth * t * t * (3 - 2 * t);
        }
      }
    }
    if (this.#pointer && o.interactive && o.barrier > 0) {
      const cx = (this.#px - this.#gx) / Cell;
      const cy = (this.#py - this.#gy) / Cell;
      const sigma = Math.max(1, o.radius / Cell / 1.5);
      const s2 = 2 * sigma * sigma;
      const reach = 3 * sigma;
      const height = o.barrier * e0;
      for (let y = Math.max(0, Math.floor(cy - reach)); y <= Math.min(N - 1, cy + reach); y++) {
        for (let x = Math.max(0, Math.floor(cx - reach)); x <= Math.min(N - 1, cx + reach); x++) {
          const dx = x - cx;
          const dy = y - cy;
          v[y * N + x] += height * Math.exp(-(dx * dx + dy * dy) / s2);
        }
      }
    }
    if (pull > 0) {
      const cx = (this.#wordX + (WordmarkGeometry.width / 2 - 1) * this.#wordScale - this.#gx) / Cell;
      const cy = (this.#wordY + (WordmarkGeometry.height / 2 - 1) * this.#wordScale - this.#gy) / Cell;
      const r = GravityRadius / Cell;
      const depth = GravityDepth * e0 * pull;
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const dx = x - cx;
          const dy = y - cy;
          const d = Math.sqrt(dx * dx + dy * dy);
          if (d < r) v[y * N + x] -= depth * (1 - d / r);
        }
      }
    }
    for (let i = 0; i < Cells; i++) {
      const phase = (-v[i] * Dt) / 2;
      const keep = Math.sqrt(edgeKeep[i]);
      this.#halfRe[i] = Math.cos(phase) * keep;
      this.#halfIm[i] = Math.sin(phase) * keep;
    }
  }

  /** ψ times a half step of the potential. */
  #half(): void {
    const re = this.#re;
    const im = this.#im;
    const hr = this.#halfRe;
    const hi = this.#halfIm;
    for (let i = 0; i < Cells; i++) {
      const a = re[i];
      const b = im[i];
      re[i] = a * hr[i] - b * hi[i];
      im[i] = a * hi[i] + b * hr[i];
    }
  }

  /** One step of the evolution (split-step Fourier). */
  #evolve(): void {
    const re = this.#re;
    const im = this.#im;
    this.#half();
    fft2(re, im, false);
    for (let i = 0; i < Cells; i++) {
      const a = re[i];
      const b = im[i];
      re[i] = a * kineticRe[i] - b * kineticIm[i];
      im[i] = a * kineticIm[i] + b * kineticRe[i];
    }
    fft2(re, im, true);
    this.#half();
  }

  /** The density and the current's velocity (∇ of the phase, clamped where the wave is thin), and their running max. */
  #current(): void {
    const re = this.#re;
    const im = this.#im;
    const rho = this.#rho;
    const vx = this.#vx;
    const vy = this.#vy;
    const cap = 2.5 * Math.max(0.2, this.#options.k);
    const cap2 = cap * cap;
    let most = 0;
    let sum = 0;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const i = y * N + x;
        const a = re[i];
        const b = im[i];
        const p = a * a + b * b;
        rho[i] = p;
        if (p > most) most = p;
        sum += p;
        this.#cdf[i] = sum;
        if (x === 0 || y === 0 || x === N - 1 || y === N - 1 || p < 1e-9) {
          vx[i] = 0;
          vy[i] = 0;
          continue;
        }
        const dax = (re[i + 1] - re[i - 1]) / 2;
        const dbx = (im[i + 1] - im[i - 1]) / 2;
        const day = (re[i + N] - re[i - N]) / 2;
        const dby = (im[i + N] - im[i - N]) / 2;
        let u = (a * dbx - b * dax) / p;
        let w = (a * dby - b * day) / p;
        const s2 = u * u + w * w;
        if (s2 > cap2) {
          const k = cap / Math.sqrt(s2);
          u *= k;
          w *= k;
        }
        vx[i] = u;
        vy[i] = w;
      }
    }
    // The picture's scale follows the densest part, falling slowly.
    this.#rhoMax = Math.max(most, this.#rhoMax * 0.985, 1e-6);
  }

  /** Bilinear lookup of a grid field at canvas px (x, y). */
  #at(f: Float32Array, x: number, y: number): number {
    let gx = (x - this.#gx) / Cell;
    let gy = (y - this.#gy) / Cell;
    if (gx < 0) gx = 0;
    else if (gx > N - 1.001) gx = N - 1.001;
    if (gy < 0) gy = 0;
    else if (gy > N - 1.001) gy = N - 1.001;
    const i = gx | 0;
    const j = gy | 0;
    const fx = gx - i;
    const fy = gy - j;
    const k = j * N + i;
    return (f[k] * (1 - fx) + f[k + 1] * fx) * (1 - fy) + (f[k + N] * (1 - fx) + f[k + N + 1] * fx) * fy;
  }

  /** A color for light going (u, w) (cells a time unit): the stream's while straight on, else its wavenumber's band. */
  #colorFor(u: number, w: number): number {
    const o = this.#options;
    const s = Math.sqrt(u * u + w * w);
    if (s < 1e-3 || Math.abs(Math.atan2(w, u)) < Straight) return 0;
    const lo = o.k * (1 - o.spectrum);
    const band = (s - lo) / Math.max(1e-3, 2 * o.k * o.spectrum);
    return 1 + Math.min(3, Math.max(0, Math.floor(band * 4)));
  }

  /** A speck where the wave is (drawn from its density); false if there is next to none. */
  #spawn(i: number): boolean {
    const total = this.#cdf[Cells - 1];
    if (total < 1e-3) return false;
    const pick = this.#random() * total;
    let lo = 0;
    let hi = Cells - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.#cdf[mid] < pick) lo = mid + 1;
      else hi = mid;
    }
    const gx = lo % N;
    const gy = (lo / N) | 0;
    this.#x[i] = this.#gx + (gx + this.#random()) * Cell;
    this.#y[i] = this.#gy + (gy + this.#random()) * Cell;
    this.#age[i] = 0;
    this.#life[i] = LifeMin + this.#random() * (LifeMax - LifeMin);
    this.#shimmer[i] = this.#random() * Math.PI * 2;
    this.#color[i] = this.#colorFor(this.#vx[lo], this.#vy[lo]);
    return true;
  }

  #step(dt: number): void {
    const o = this.#options;
    this.#time += dt;
    this.#glass?.step(dt);
    const glass = this.#glass;
    const pull = glass ? glass.pull : 0;
    const whole = !glass || glass.whole;
    if (whole !== this.#wasWhole || pull > 0 || this.#wasPull > 0) this.#potential = true;
    this.#wasWhole = whole;
    this.#wasPull = pull;
    if (this.#potential) {
      this.#buildPotential(pull);
      this.#potential = false;
    }

    // The light comes in: packets at the left, toward the star, each its own wavenumber (faster at the start).
    const sx = (this.#ox - this.#gx) / Cell;
    const sy = (this.#oy - this.#gy) / Cell;
    const intro = this.#time < Intro ? 1 + 2 * (1 - this.#time / Intro) : 1;
    this.#emit += dt * o.rate * intro;
    while (this.#emit >= 1) {
      this.#emit -= 1;
      const y = sy + (this.#random() * 2 - 1) * 2;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      this.#packet(Source, y, Math.atan2(sy - y, sx - Source), k, 1);
    }
    // The glass back together: a little of the light it held, out along its arms.
    while (this.#bursts >= 1) {
      this.#bursts -= 1;
      const arm = (Math.floor(this.#random() * 8) * Math.PI) / 4;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      this.#packet(sx + Math.cos(arm) * 4, sy + Math.sin(arm) * 4, arm, k, 0.8);
    }

    // The evolution: as many steps as the time that passed takes.
    this.#clock += dt * this.#rate();
    let steps = 0;
    while (this.#clock >= Dt && steps < MaxSteps) {
      this.#clock -= Dt;
      this.#evolve();
      steps++;
    }
    if (steps === MaxSteps) this.#clock = 0;
    this.#current();

    // The specks ride the current (a midpoint step), and let go where the wave thins.
    const toPx = this.#rate() * Cell;
    const target = Math.min(Capacity, Math.max(0, Math.round(o.specks)));
    const thin = this.#rhoMax * Thin;
    const x0 = this.#gx;
    const y0 = this.#gy;
    const x1 = this.#gx + N * Cell;
    const y1 = this.#gy + N * Cell;
    for (let i = 0; i < this.#alive; i++) {
      const x = this.#x[i];
      const y = this.#y[i];
      const u1 = this.#at(this.#vx, x, y) * toPx;
      const w1 = this.#at(this.#vy, x, y) * toPx;
      const mx = x + (u1 * dt) / 2;
      const my = y + (w1 * dt) / 2;
      const u2 = this.#at(this.#vx, mx, my);
      const w2 = this.#at(this.#vy, mx, my);
      const nx = x + u2 * toPx * dt;
      const ny = y + w2 * toPx * dt;
      this.#x[i] = nx;
      this.#y[i] = ny;
      this.#age[i] += dt;
      if (this.#color[i] === 0) this.#color[i] = this.#colorFor(u2, w2);
      if (
        this.#age[i] > this.#life[i] ||
        nx < x0 ||
        ny < y0 ||
        nx >= x1 ||
        ny >= y1 ||
        this.#at(this.#rho, nx, ny) < thin
      ) {
        if (this.#alive > target || !this.#spawn(i)) {
          this.#kill(i);
          i--;
        }
      }
    }
    // More where there are fewer than asked for (a few a frame, so they come in gradually).
    for (let n = 0; n < 24 && this.#alive < target; n++) {
      if (!this.#spawn(this.#alive)) break;
      this.#alive++;
    }
  }

  #kill(i: number): void {
    const last = --this.#alive;
    this.#x[i] = this.#x[last];
    this.#y[i] = this.#y[last];
    this.#age[i] = this.#age[last];
    this.#life[i] = this.#life[last];
    this.#shimmer[i] = this.#shimmer[last];
    this.#color[i] = this.#color[last];
  }

  /** The wave's picture: |ψ|² in its light's colors, the wavefronts (ψ's real part) in it. */
  #paint(): void {
    const image = this.#image;
    if (!image || !this.#pictureCtx) return;
    const o = this.#options;
    const data = image.data;
    const rho = this.#rho;
    const inv = 1 / this.#rhoMax;
    const glow = Math.min(1, o.glow) * 255;
    for (let i = 0; i < Cells; i++) {
      const p = rho[i] * inv;
      const j = i * 4;
      if (p < 1e-4) {
        data[j + 3] = 0;
        continue;
      }
      const c = o.colors ? this.#colorFor(this.#vx[i], this.#vy[i]) * 3 : 0;
      data[j] = this.#rgb[c];
      data[j + 1] = this.#rgb[c + 1];
      data[j + 2] = this.#rgb[c + 2];
      const front = 1 - o.ripples * 0.5 * (1 - this.#re[i] / Math.sqrt(rho[i]));
      data[j + 3] = glow * Math.min(1, p ** 0.7) * front * edgeShow[i];
    }
    this.#pictureCtx.putImageData(image, 0, 0);
  }

  #draw(): void {
    const ctx = this.#ctx;
    if (!ctx) return;
    const o = this.#options;
    const dpr = this.#dpr;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, this.#canvas.width, this.#canvas.height);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    if (o.glow > 0 && this.#picture) {
      this.#paint();
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(this.#picture, this.#gx, this.#gy, N * Cell, N * Cell);
    }

    // The specks, bucketed by color and alpha level (counting sort): one path each.
    const counts = this.#counts;
    counts.fill(0);
    const inv = 1 / this.#rhoMax;
    for (let i = 0; i < this.#alive; i++) {
      const p = this.#at(this.#rho, this.#x[i], this.#y[i]) * inv;
      const age = this.#age[i];
      const fade = Math.min(1, age / 0.4, Math.max(0, (this.#life[i] - age) / 0.6));
      const a = Math.min(1, Math.max(0, p) ** 0.6 * 1.6) * fade * this.#at(edgeShow, this.#x[i], this.#y[i]);
      const level = Math.max(0, Math.min(Levels - 1, Math.round(a * (Levels - 1))));
      const c = o.colors ? this.#color[i] : 0;
      const b = c * Levels + level;
      this.#bucketOf[i] = b;
      counts[b]++;
    }
    let sum = 0;
    for (let b = 0; b < counts.length; b++) {
      this.#starts[b] = sum;
      sum += counts[b];
    }
    counts.fill(0);
    for (let i = 0; i < this.#alive; i++) {
      const b = this.#bucketOf[i];
      this.#bucket[this.#starts[b] + counts[b]++] = i;
    }
    const t = this.#time;
    for (let c = 0; c < Colors; c++) {
      ctx.fillStyle = this.#colors[c];
      for (let l = 1; l < Levels; l++) {
        const b = c * Levels + l;
        const n = counts[b];
        if (!n) continue;
        ctx.globalAlpha = ((c === 0 ? 0.5 : 1) * o.alpha * l) / (Levels - 1);
        ctx.beginPath();
        const s = this.#starts[b];
        for (let k = 0; k < n; k++) {
          const i = this.#bucket[s + k];
          const r = o.size * (c === 0 ? 0.7 : 1) * (1 + 0.25 * Math.sin(t * 2.2 + this.#shimmer[i]));
          const x = this.#x[i];
          const y = this.#y[i];
          ctx.moveTo(x + r, y);
          ctx.arc(x, y, r, 0, Math.PI * 2);
        }
        ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
    this.#glass?.draw(
      ctx,
      this.#wordX,
      this.#wordY,
      this.#wordScale,
      this.#gradient ?? this.#colors[1],
      this.#colors[0],
      1,
    );
  }
}
