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
// the other half. Emitters at the left radiate wave packets, each of its own wavenumber (its color), the flow as strong
// as the word is assembled (none while its pieces fly in; dimmer while it is broken). The wordmark is
// glass, a potential well carried by its pieces (`Glass`: the filled word, piece by piece where each is now): a packet's
// refraction through it depends on its energy, so the star splits the light by wavenumber as a prism does; a shatter
// scatters the well with the pieces and the snap back puts it together (drawing the light in for a moment). The specks
// are the wave's particles: |ψ|² is the chance of one being born there (as many a second as the wave's amount says);
// each then rides the wave (a Bohmian trajectory: it moves with the probability current, ∇ of ψ's phase, a midpoint
// Runge–Kutta step), fading with a half-life, colored by the wavenumber it carries once it leaves the straight way.
// Under them the wave itself, |ψ|², faintly, its wavefronts in it (added light on a dark canvas). The pointer is a
// barrier the wave goes round, a press kicks it outward. The grid's edges absorb (no wrapping round).
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
  /**
   * The glass: the word's well's depth, in a mean packet's energies. Its pieces carry it (a shatter scatters it, the
   * snap back puts it together); the star is the prism.
   */
  prism: number;
  /** The wave itself shows (|ψ|², 0: only the specks). */
  glow: number;
  /** The wavefronts show in it (0: the density alone). */
  ripples: number;
  /**
   * Specks born a second, at the wave's usual amount (more while there is more of it): each where the wave is, as
   * likely as it is dense there (|ψ|² is the chance of a speck), then riding it.
   */
  spawn: number;
  /** A speck's half-life (s): it fades as it rides, whatever the wave under it. */
  decay: number;
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
  spawn: 160,
  decay: 1.6,
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
/**
 * The word as glass, sampled once: samples a drawing unit, and how much it is softened (units): the well follows the
 * letters and the star, filled, a little rounded.
 */
const MaskRes = 4;
const MaskSoft = 0.8;
/**
 * What shows under the glass: the wave pools in its well (that is what refracts it), but in the picture the word stays
 * clear: no glow, no births, specks fading, where the glass covers a cell (its cover times this, at most 1).
 */
const CoverHide = 1.6;
/** Without the glass (reduced motion, no Path2D): the star alone, a well with eight arms of this radius (px). */
const PrismRadius = 34;
const PrismArms = 0.45;
const PrismSoft = 1.5;
/**
 * The light comes in from emitters, not all at once: each radiates one packet over `EmitTime` (the grid's time units),
 * its strength easing in and out, at the packet's own frequency (ω = k²/2) from a source `EmitAlong` × `EmitAcross`
 * cells (half sizes) whose phase runs along its way, so it radiates forward only. `EmitStrength` makes a packet about
 * as dense as a unit one. At most `Emitters` at once; the light comes in `Source` cells from the left edge, and shows
 * gradually over the next `BeamFade`.
 */
const Emitters = 12;
const EmitTime = 20;
const EmitAlong = 6;
const EmitAcross = 4;
const EmitStrength = 0.07;
const Source = Edge + 6;
const BeamFade = 26;
/** How fast a cell's color follows its light's (1/s). */
const TintEase = 5;
/** How much stronger the wave shows on a dark canvas. */
const DarkGlow = 3;
/** How fast the picture's scale follows the wave's densest part (1/s): no jump when a packet comes in. */
const ScaleEase = 1.2;
/** The snap back's well: how far it reaches (px) and how deep it is (in a mean packet's energies, at full pull). */
const GravityRadius = 320;
const GravityDepth = 1.6;
/** Specks at most; half-lives a speck lasts (it is gone at 1/32); seconds it takes to come in. */
const Capacity = 1500;
const Lifetimes = 5;
const Birth = 0.3;
/** How fast the wave's usual amount follows its amount now (1/s): the births' measure. */
const AmountEase = 0.25;
/** A speck going within this of straight on (radians) is light not yet split: the stream's color. */
const Straight = 0.22;
/** Alpha levels a frame batches specks into (per color). */
const Levels = 10;
const Colors = 5;
/**
 * The light follows the glass: the flow from the left is as strong as the word is assembled (`Glass.assembled`; none
 * while its pieces fly in, growing as they come together, dimming when it breaks), eased at this rate (1/s).
 */
const FlowEase = 3;
/**
 * The least the picture's scale is (|ψ|², about a full packet's): a weaker flow shows weaker instead of being scaled
 * back up.
 */
const RhoFloor = 0.5;

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
    // The beam shows gradually from where it comes in, as the light's stream does.
    const beam = Math.min(1, Math.max(0, (x - Source + 4) / BeamFade));
    edgeShow[y * N + x] = inside * inside * (3 - 2 * inside) * beam * beam * (3 - 2 * beam);
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
  #rhoMax = RhoFloor;
  /** How strong the flow from the left is now (0..1, the glass's assembly, eased). */
  #flow = 0;
  #potential = true;
  #clock = 0;
  #emit = 0;
  /** The evolution's own time (the grid's units). */
  #sim = 0;
  // The emitters: where (cells), which way, the wavenumber, the strength, when each started (`#sim`).
  readonly #emX = new Float32Array(Emitters);
  readonly #emY = new Float32Array(Emitters);
  readonly #emCos = new Float32Array(Emitters);
  readonly #emSin = new Float32Array(Emitters);
  readonly #emK = new Float32Array(Emitters);
  readonly #emAmp = new Float32Array(Emitters);
  readonly #emT0 = new Float32Array(Emitters);
  /** 1 for the flow's emitters (they follow the glass), 0 for the snap's bursts. */
  readonly #emFlow = new Uint8Array(Emitters);
  #emitters = 0;

  // The specks: position (canvas px), age (s), color (0: the stream's; 1–4: the logo's), a phase for their size; the
  // births owed, and the wave's usual amount (Σ|ψ|²).
  readonly #x = new Float32Array(Capacity);
  readonly #y = new Float32Array(Capacity);
  readonly #age = new Float32Array(Capacity);
  #births = 0;
  #amount = 0;
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
  /** The picture's color per cell, eased toward its light's (no flicker where the phase is noise); set once colors are read. */
  readonly #tintR = new Float32Array(Cells);
  readonly #tintG = new Float32Array(Cells);
  readonly #tintB = new Float32Array(Cells);
  #tinted = false;
  #frameDt = 1 / 60;
  /** On a dark canvas (its stream color light): the wave is light added to it, stronger. */
  #dark = false;
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
  /** The filled word (0..1, `MaskRes` samples a unit from (−1, −1)), and a piece's place (`Glass.place`). */
  #mask: Float32Array | undefined;
  #maskWidth = 0;
  #maskHeight = 0;
  readonly #pose = new Float32Array(5);
  /** How much of each cell the glass covers now (0..1, from its pieces; `CoverHide`). */
  readonly #cover = new Float32Array(Cells);
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
      this.#emitters = 0;
      this.#births = 0;
      this.#amount = 0;
      this.#flow = 0;
      this.#rhoMax = RhoFloor;
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
    this.#dark = 0.299 * this.#rgb[0] + 0.587 * this.#rgb[1] + 0.114 * this.#rgb[2] > 128;
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
    this.#makeMask();
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
    this.#frameDt = dt;
    this.#step(dt);
    this.#draw();
    this.#frame = requestAnimationFrame(this.#tick);
  };

  #drawStill(): void {
    if (!this.#width) return;
    this.#re.fill(0);
    this.#im.fill(0);
    this.#alive = 0;
    this.#emitters = 0;
    this.#births = 0;
    this.#amount = 0;
    this.#flow = 0;
    this.#rhoMax = RhoFloor;
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

  /** An emitter at grid (cx, cy) radiating a packet going `angle` with wavenumber `k` (none if all are busy). */
  #launch(cx: number, cy: number, angle: number, k: number, amplitude: number, flow: boolean): void {
    const e = this.#emitters;
    if (e >= Emitters) return;
    this.#emX[e] = cx;
    this.#emY[e] = cy;
    this.#emCos[e] = Math.cos(angle);
    this.#emSin[e] = Math.sin(angle);
    this.#emK[e] = k;
    this.#emAmp[e] = amplitude;
    this.#emT0[e] = this.#sim;
    this.#emFlow[e] = flow ? 1 : 0;
    this.#emitters = e + 1;
  }

  /** A step of the emitters' radiation: each adds its source, easing in and out over its time, then is done. */
  #radiate(): void {
    const re = this.#re;
    const im = this.#im;
    const a2 = 2 * EmitAlong * EmitAlong;
    const c2 = 2 * EmitAcross * EmitAcross;
    const reach = 3 * Math.max(EmitAlong, EmitAcross);
    for (let e = 0; e < this.#emitters; e++) {
      const tau = (this.#sim - this.#emT0[e]) / EmitTime;
      if (tau >= 1) {
        const last = --this.#emitters;
        this.#emX[e] = this.#emX[last];
        this.#emY[e] = this.#emY[last];
        this.#emCos[e] = this.#emCos[last];
        this.#emSin[e] = this.#emSin[last];
        this.#emK[e] = this.#emK[last];
        this.#emAmp[e] = this.#emAmp[last];
        this.#emT0[e] = this.#emT0[last];
        this.#emFlow[e] = this.#emFlow[last];
        e--;
        continue;
      }
      const ease = Math.sin(Math.PI * tau);
      const strength = this.#emAmp[e] * EmitStrength * ease * ease * Dt * (this.#emFlow[e] ? this.#flow : 1);
      const k = this.#emK[e];
      const c = this.#emCos[e];
      const sn = this.#emSin[e];
      const cx = this.#emX[e];
      const cy = this.#emY[e];
      // Its phase now: e^{−iωt}, ω = k²/2.
      const t = this.#sim - this.#emT0[e];
      const wt = (-k * k * t) / 2;
      const x0 = Math.max(1, Math.floor(cx - reach));
      const x1 = Math.min(N - 2, Math.ceil(cx + reach));
      const y0 = Math.max(1, Math.floor(cy - reach));
      const y1 = Math.min(N - 2, Math.ceil(cy + reach));
      for (let y = y0; y <= y1; y++) {
        const dy = y - cy;
        for (let x = x0; x <= x1; x++) {
          const dx = x - cx;
          const along = dx * c + dy * sn;
          const across = -dx * sn + dy * c;
          const g = strength * Math.exp(-(along * along) / a2 - (across * across) / c2);
          if (g < 1e-6) continue;
          const phase = k * along + wt;
          const i = y * N + x;
          re[i] += g * Math.cos(phase);
          im[i] += g * Math.sin(phase);
        }
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

  /** The word filled (its letters and star as solid glass), softened a little, as a mask in drawing units. */
  #makeMask(): void {
    if (typeof document === 'undefined') return;
    const w = (WordmarkGeometry.width + 2) * MaskRes;
    const h = (WordmarkGeometry.height + 2) * MaskRes;
    const c = document.createElement('canvas');
    c.width = w;
    c.height = h;
    const g = c.getContext('2d', { willReadFrequently: true });
    if (!g) return;
    g.filter = `blur(${MaskSoft * MaskRes}px)`;
    g.setTransform(MaskRes, 0, 0, MaskRes, MaskRes, MaskRes);
    g.fillStyle = '#000';
    g.strokeStyle = '#000';
    g.lineWidth = WordmarkGeometry.stroke(3);
    g.lineJoin = 'miter';
    for (const d of WordmarkGeometry.paths) {
      const path = new Path2D(d);
      g.fill(path);
      g.stroke(path);
    }
    const data = g.getImageData(0, 0, w, h).data;
    const mask = new Float32Array(w * h);
    for (let i = 0; i < mask.length; i++) mask[i] = data[i * 4 + 3] / 255;
    this.#mask = mask;
    this.#maskWidth = w;
    this.#maskHeight = h;
  }

  /** The glass's well: each grid cell (2 × 2 samples) inside a piece takes the word where that piece came from. */
  #glassPotential(depth: number): void {
    const glass = this.#glass;
    const mask = this.#mask;
    if (!glass || !mask) return;
    const v = this.#v;
    const pose = this.#pose;
    const corners = glass.corners;
    const scale = this.#wordScale;
    const mw = this.#maskWidth;
    const mh = this.#maskHeight;
    for (let i = 0; i < glass.count; i++) {
      glass.place(i, pose);
      const mx = pose[0];
      const my = pose[1];
      const c = Math.cos(pose[2]);
      const sn = Math.sin(pose[2]);
      const hx = pose[3];
      const hy = pose[4];
      const ax = corners[i * 6];
      const ay = corners[i * 6 + 1];
      const bx = corners[i * 6 + 2];
      const by = corners[i * 6 + 3];
      const qx = corners[i * 6 + 4];
      const qy = corners[i * 6 + 5];
      // The piece's box on the grid, from its corners where they are now.
      let x0 = N;
      let x1 = -1;
      let y0 = N;
      let y1 = -1;
      for (let k = 0; k < 3; k++) {
        const lx = corners[i * 6 + 2 * k];
        const ly = corners[i * 6 + 2 * k + 1];
        const gx = (this.#wordX + (mx + lx * c - ly * sn) * scale - this.#gx) / Cell;
        const gy = (this.#wordY + (my + lx * sn + ly * c) * scale - this.#gy) / Cell;
        if (gx < x0) x0 = gx;
        if (gx > x1) x1 = gx;
        if (gy < y0) y0 = gy;
        if (gy > y1) y1 = gy;
      }
      const cx0 = Math.max(0, Math.floor(x0));
      const cx1 = Math.min(N - 1, Math.ceil(x1));
      const cy0 = Math.max(0, Math.floor(y0));
      const cy1 = Math.min(N - 1, Math.ceil(y1));
      for (let y = cy0; y <= cy1; y++) {
        for (let x = cx0; x <= cx1; x++) {
          let sum = 0;
          for (let sy = 0; sy < 2; sy++) {
            for (let sx = 0; sx < 2; sx++) {
              // The sample in the piece's own frame (drawing units about its middle).
              const du = (this.#gx + (x + 0.25 + 0.5 * sx) * Cell - this.#wordX) / scale - mx;
              const dv = (this.#gy + (y + 0.25 + 0.5 * sy) * Cell - this.#wordY) / scale - my;
              const lx = du * c + dv * sn;
              const ly = -du * sn + dv * c;
              const e1 = (bx - ax) * (ly - ay) - (by - ay) * (lx - ax);
              const e2 = (qx - bx) * (ly - by) - (qy - by) * (lx - bx);
              const e3 = (ax - qx) * (ly - qy) - (ay - qy) * (lx - qx);
              if (!((e1 >= 0 && e2 >= 0 && e3 >= 0) || (e1 <= 0 && e2 <= 0 && e3 <= 0))) continue;
              // Where in the word it came from.
              const mi = Math.round((hx + lx + 1) * MaskRes);
              const mj = Math.round((hy + ly + 1) * MaskRes);
              if (mi >= 0 && mj >= 0 && mi < mw && mj < mh) sum += mask[mj * mw + mi];
            }
          }
          if (sum > 0) {
            v[y * N + x] -= (depth * sum) / 4;
            if (sum / 4 > this.#cover[y * N + x]) this.#cover[y * N + x] = sum / 4;
          }
        }
      }
    }
  }

  /** The potential (the glass, or the star alone without it; the pointer's barrier, the snap back's well) and its half step. */
  #buildPotential(pull: number): void {
    const o = this.#options;
    const v = this.#v;
    v.fill(0);
    this.#cover.fill(0);
    const e0 = (o.k * o.k) / 2;
    const sx = (this.#ox - this.#gx) / Cell;
    const sy = (this.#oy - this.#gy) / Cell;
    if (this.#glass && this.#mask) this.#glassPotential(o.prism * e0);
    else {
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
    this.#radiate();
    this.#sim += Dt;
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
  #current(dt: number): void {
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
        // What shows (none under the glass) sets the scale and the births (by |ψ|²).
        const shown = p * (1 - Math.min(1, this.#cover[i] * CoverHide));
        if (shown > most) most = shown;
        sum += shown;
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
    // The picture's scale eases toward the densest part, either way (a packet coming in never dims the rest at once).
    this.#rhoMax = Math.max(RhoFloor, this.#rhoMax + (most - this.#rhoMax) * Math.min(1, dt * ScaleEase));
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
    this.#flow += ((glass ? glass.assembled : 1) - this.#flow) * Math.min(1, dt * FlowEase);
    const whole = !glass || glass.whole;
    // The pieces moving (or leaning) move the glass's well with them.
    if (whole !== this.#wasWhole || pull > 0 || this.#wasPull > 0 || (glass && !glass.resting)) this.#potential = true;
    this.#wasWhole = whole;
    this.#wasPull = pull;
    if (this.#potential) {
      this.#buildPotential(pull);
      this.#potential = false;
    }

    // The light comes in: packets at the left, toward the star, each its own wavenumber (faster at the start).
    const sx = (this.#ox - this.#gx) / Cell;
    const sy = (this.#oy - this.#gy) / Cell;
    this.#emit += dt * o.rate;
    while (this.#emit >= 1) {
      this.#emit -= 1;
      const y = sy + (this.#random() * 2 - 1) * 2;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      this.#launch(Source, y, Math.atan2(sy - y, sx - Source), k, 1, true);
    }
    // The glass back together: a little of the light it held, out along its arms.
    while (this.#bursts >= 1) {
      this.#bursts -= 1;
      const arm = (Math.floor(this.#random() * 8) * Math.PI) / 4;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      this.#launch(sx + Math.cos(arm) * 4, sy + Math.sin(arm) * 4, arm, k, 0.8, false);
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
    this.#current(dt);

    // The specks ride the current (a midpoint step), fading as they go; gone when faded or off the grid.
    const toPx = this.#rate() * Cell;
    const gone = Lifetimes * Math.max(0.1, o.decay);
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
      if (this.#age[i] > gone || nx < x0 || ny < y0 || nx >= x1 || ny >= y1) {
        this.#kill(i);
        i--;
      }
    }
    // Births: |ψ|² is the chance of a speck. As many a second as asked at the wave's usual amount, more with more of it,
    // each drawn from its density.
    const amount = this.#cdf[Cells - 1];
    // The usual amount is learned while the flow is full; until then (the first assembly) births follow the flow.
    if (this.#flow > 0.9) {
      this.#amount = this.#amount ? this.#amount + (amount - this.#amount) * Math.min(1, dt * AmountEase) : amount;
    }
    this.#births += dt * o.spawn * (this.#amount > 1e-6 ? amount / this.#amount : this.#flow);
    while (this.#births >= 1) {
      this.#births -= 1;
      if (this.#alive >= Capacity || !this.#spawn(this.#alive)) continue;
      this.#alive++;
    }
  }

  #kill(i: number): void {
    const last = --this.#alive;
    this.#x[i] = this.#x[last];
    this.#y[i] = this.#y[last];
    this.#age[i] = this.#age[last];
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
    const glow = Math.min(1, o.glow) * (this.#dark ? DarkGlow : 1) * 255;
    const rgb = this.#rgb;
    const tr = this.#tintR;
    const tg = this.#tintG;
    const tb = this.#tintB;
    if (!this.#tinted) {
      tr.fill(rgb[0]);
      tg.fill(rgb[1]);
      tb.fill(rgb[2]);
      this.#tinted = true;
    }
    const ease = Math.min(1, this.#frameDt * TintEase);
    const lo = o.k * (1 - o.spectrum);
    const span = Math.max(1e-3, 2 * o.k * o.spectrum);
    for (let i = 0; i < Cells; i++) {
      const p = rho[i] * inv;
      const j = i * 4;
      if (p < 1e-4) {
        data[j + 3] = 0;
        continue;
      }
      // Its light's color: the stream's while straight on, else between the logo's colors by its wavenumber; trusted
      // only where the wave is dense (where it is thin its phase is noise).
      let r = rgb[0];
      let g = rgb[1];
      let b = rgb[2];
      if (o.colors) {
        const u = this.#vx[i];
        const w = this.#vy[i];
        const turn = Math.abs(Math.atan2(w, u));
        let colorful = Math.min(1, Math.max(0, (turn - Straight * 0.5) / Straight)) * Math.min(1, p * 12);
        colorful = colorful * colorful * (3 - 2 * colorful);
        if (colorful > 0) {
          const f = Math.min(2.999, Math.max(0, ((Math.sqrt(u * u + w * w) - lo) / span) * 4 - 0.5));
          const k = f | 0;
          const t = f - k;
          const a = (1 + k) * 3;
          const z = (2 + k) * 3;
          r += (rgb[a] * (1 - t) + rgb[z] * t - r) * colorful;
          g += (rgb[a + 1] * (1 - t) + rgb[z + 1] * t - g) * colorful;
          b += (rgb[a + 2] * (1 - t) + rgb[z + 2] * t - b) * colorful;
        }
      }
      tr[i] += (r - tr[i]) * ease;
      tg[i] += (g - tg[i]) * ease;
      tb[i] += (b - tb[i]) * ease;
      data[j] = tr[i];
      data[j + 1] = tg[i];
      data[j + 2] = tb[i];
      const front = 1 - o.ripples * 0.5 * (1 - this.#re[i] / Math.sqrt(rho[i]));
      const clear = 1 - Math.min(1, this.#cover[i] * CoverHide);
      data[j + 3] = Math.min(255, glow * Math.min(1, p ** 0.5) * front * edgeShow[i] * clear);
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
      // On a dark canvas the wave is light added to it (a veil of grey would sink into the black).
      if (this.#dark) ctx.globalCompositeOperation = 'lighter';
      ctx.drawImage(this.#picture, this.#gx, this.#gy, N * Cell, N * Cell);
      ctx.globalCompositeOperation = 'source-over';
    }

    // The specks, bucketed by color and alpha level (counting sort): one path each.
    const counts = this.#counts;
    counts.fill(0);
    const decay = Math.max(0.1, o.decay);
    for (let i = 0; i < this.#alive; i++) {
      const age = this.#age[i];
      // In over its birth, then halving every `decay` seconds (not the wave under it: it is a particle now).
      const x = this.#x[i];
      const y = this.#y[i];
      const clear = 1 - Math.min(1, this.#at(this.#cover, x, y) * CoverHide);
      const a = Math.min(1, age / Birth) * 0.5 ** (age / decay) * this.#at(edgeShow, x, y) * clear;
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
