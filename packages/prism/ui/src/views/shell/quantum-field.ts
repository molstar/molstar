/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The empty canvas's splash as a quantum system (beside prism-field.ts, the same surface; shell.md §1): the light is
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
import { GlassFootprint } from './glass-footprint.js';
import { cancelFrameCallback, requestFrameCallback } from '../../reactive/frame-state.js';
import { AmbientPace, type AmbientPaceOptions } from './ambient-pace.js';

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
  rate: 3.4,
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
const Dt = 0.15;
const MaxSteps = 4;
/** The absorbing edge: cells wide, and how much it takes a step at the very edge. */
const Edge = 14;
const Absorb = 0.18;
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
 * as dense as a unit one. At most `Emitters` at once; the light comes in `Source` cells from the left edge.
 */
const Emitters = 12;
const EmitTime = 20;
const EmitAlong = 7;
const EmitAcross = 8;
const EmitStrength = 0.07;
/**
 * The light comes in as a broad front, not a beam: packets start anywhere within this many cells of the star's line
 * (about the word's height) and all go straight right, so no one spot holds the light; the glass bends it.
 */
const FrontHalf = 14;
/** A packet's amplitude in the front, and in the snap's burst (soft: no spike at the star). */
const FrontAmplitude = 0.8;
/**
 * The start: the flow already under way, faint: this many packets (half size `SeedSize` cells, amplitude
 * `SeedAmplitude`, wavenumbers of the spectrum) spread along the beam's whole length and across its band, going right
 * within `SeedSpread` radians, so the picture flows from the first frame and the front strengthens it as the glass
 * assembles; it comes in over `SeedFade` s.
 */
const Seeds = 10;
const SeedSize = 10;
const SeedAmplitude = 0.5;
const SeedSpread = 0.25;
const SeedFade = 1.5;
const BurstAmplitude = 0.45;
const Source = Edge + 6;
/** How much of the vignette's radius the wave fades over toward its rim. */
const VignetteSoft = 0.45;
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
/** The same as `var(…)` values, made once (the colors are read every 1.5 s). */
const ColorRefs = ColorVars.map((v) => `var(${v})`);

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
const twiddleRe = new Float64Array(N / 2);
const twiddleIm = new Float64Array(N / 2);
for (let j = 0; j < N / 2; j++) {
  twiddleRe[j] = Math.cos((-2 * Math.PI * j) / N);
  twiddleIm[j] = Math.sin((-2 * Math.PI * j) / N);
}

// The split step needs no ordered spectrum: the forward transform is decimation in frequency (natural order in,
// bit-reversed out), the inverse decimation in time (bit-reversed in, natural out), and the kinetic table is kept in
// the bit-reversed order, so no permutation is ever done. Rows run along a row; columns run a whole row at a time
// (contiguous, the butterfly's twiddle fixed across it).

/** Forward 2D FFT, decimation in frequency: natural order in, bit-reversed out (both axes). */
function forward2(re: Float64Array, im: Float64Array): void {
  for (let o = 0; o < Cells; o += N) {
    for (let size = N; size >= 2; size >>= 1) {
      const half = size >> 1;
      const step = N / size;
      for (let start = 0; start < N; start += size) {
        for (let k = 0; k < half; k++) {
          const wr = twiddleRe[k * step];
          const wi = twiddleIm[k * step];
          const a = o + start + k;
          const b = a + half;
          const ar = re[a];
          const ai = im[a];
          const dr = ar - re[b];
          const di = ai - im[b];
          re[a] = ar + re[b];
          im[a] = ai + im[b];
          re[b] = dr * wr - di * wi;
          im[b] = dr * wi + di * wr;
        }
      }
    }
  }
  for (let size = N; size >= 2; size >>= 1) {
    const half = size >> 1;
    const step = N / size;
    for (let start = 0; start < N; start += size) {
      for (let k = 0; k < half; k++) {
        const wr = twiddleRe[k * step];
        const wi = twiddleIm[k * step];
        const ra = (start + k) * N;
        const rb = ra + half * N;
        for (let x = 0; x < N; x++) {
          const a = ra + x;
          const b = rb + x;
          const ar = re[a];
          const ai = im[a];
          const dr = ar - re[b];
          const di = ai - im[b];
          re[a] = ar + re[b];
          im[a] = ai + im[b];
          re[b] = dr * wr - di * wi;
          im[b] = dr * wi + di * wr;
        }
      }
    }
  }
}

/** Inverse 2D FFT (unscaled), decimation in time: bit-reversed in (both axes), natural order out. */
function inverse2(re: Float64Array, im: Float64Array): void {
  for (let size = 2; size <= N; size <<= 1) {
    const half = size >> 1;
    const step = N / size;
    for (let start = 0; start < N; start += size) {
      for (let k = 0; k < half; k++) {
        const wr = twiddleRe[k * step];
        const wi = -twiddleIm[k * step];
        const ra = (start + k) * N;
        const rb = ra + half * N;
        for (let x = 0; x < N; x++) {
          const a = ra + x;
          const b = rb + x;
          const br = re[b] * wr - im[b] * wi;
          const bi = re[b] * wi + im[b] * wr;
          re[b] = re[a] - br;
          im[b] = im[a] - bi;
          re[a] += br;
          im[a] += bi;
        }
      }
    }
  }
  for (let o = 0; o < Cells; o += N) {
    for (let size = 2; size <= N; size <<= 1) {
      const half = size >> 1;
      const step = N / size;
      for (let start = 0; start < N; start += size) {
        for (let k = 0; k < half; k++) {
          const wr = twiddleRe[k * step];
          const wi = -twiddleIm[k * step];
          const a = o + start + k;
          const b = a + half;
          const br = re[b] * wr - im[b] * wi;
          const bi = re[b] * wi + im[b] * wr;
          re[b] = re[a] - br;
          im[b] = im[a] - bi;
          re[a] += br;
          im[a] += bi;
        }
      }
    }
  }
}

/** Bit-reverses both axes in place (the tests' bridge to the natural order; the evolution never needs it). */
function permute2(re: Float64Array, im: Float64Array): void {
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const i = y * N + x;
      const j = reversed[y] * N + reversed[x];
      if (j <= i) continue;
      let t = re[i];
      re[i] = re[j];
      re[j] = t;
      t = im[i];
      im[i] = im[j];
      im[j] = t;
    }
  }
}

/** The 2D FFT of the N × N grid in natural order (inverse: unscaled). For the tests; the evolution skips the orders. */
export function fft2(re: Float64Array, im: Float64Array, inverse: boolean): void {
  if (inverse) {
    permute2(re, im);
    inverse2(re, im);
  } else {
    forward2(re, im);
    permute2(re, im);
  }
}

/** atan2 to about 1e-5 radians, much cheaper than `Math.atan2` (the current's phase differences, 4 a cell a frame). */
function atan2(y: number, x: number): number {
  const ax = Math.abs(x);
  const ay = Math.abs(y);
  const big = ax > ay ? ax : ay;
  if (big === 0) return 0;
  const a = (ax < ay ? ax : ay) / big;
  const s = a * a;
  let r = ((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;
  if (ay > ax) r = 1.5707963267948966 - r;
  if (x < 0) r = 3.141592653589793 - r;
  return y < 0 ? -r : r;
}

// The kinetic step in momentum space, e^{−i k² Dt / 2}, with the inverse FFT's 1 / N² folded in, in the forward
// transform's bit-reversed order; and the absorbing edge's share a step (1 inside), and its square root (a half step).
const kineticRe = new Float64Array(Cells);
const kineticIm = new Float64Array(Cells);
const edgeKeep = new Float32Array(Cells);
const edgeHalf = new Float32Array(Cells);
/** How much of the wave shows (1 inside, easing to 0 across the absorbing edge): the grid's square never shows. */
const edgeShow = new Float32Array(Cells);
for (let y = 0; y < N; y++) {
  const fy = reversed[y];
  const ky = (2 * Math.PI * (fy < N / 2 ? fy : fy - N)) / N;
  const ey = Math.max(0, Edge - Math.min(y, N - 1 - y)) / Edge;
  for (let x = 0; x < N; x++) {
    const fx = reversed[x];
    const kx = (2 * Math.PI * (fx < N / 2 ? fx : fx - N)) / N;
    const phase = (-(kx * kx + ky * ky) * Dt) / 2;
    kineticRe[y * N + x] = Math.cos(phase) / Cells;
    kineticIm[y * N + x] = Math.sin(phase) / Cells;
    const ex = Math.max(0, Edge - Math.min(x, N - 1 - x)) / Edge;
    const e = Math.max(ex, ey);
    edgeKeep[y * N + x] = 1 - Absorb * e * e;
    edgeHalf[y * N + x] = Math.sqrt(edgeKeep[y * N + x]);
    // A round vignette, clear of the square absorbing edge: the grid's shape never shows.
    const r = Math.hypot(x + 0.5 - N / 2, y + 0.5 - N / 2) / (N / 2 - Edge);
    const inside = Math.min(1, Math.max(0, (1 - r) / VignetteSoft));
    edgeShow[y * N + x] = inside * inside * (3 - 2 * inside);
  }
}
/** The cells inside the vignette (off the grid's border): the only ones the picture, the current and births need. */
const inner = (() => {
  const list: number[] = [];
  for (let y = 1; y < N - 1; y++) for (let x = 1; x < N - 1; x++) if (edgeShow[y * N + x] > 0) list.push(y * N + x);
  return Int32Array.from(list);
})();

export class QuantumField {
  #options: QuantumFieldOptions;
  readonly #canvas: HTMLCanvasElement;
  readonly #ctx: CanvasRenderingContext2D | null;
  /** The canvas's computed style, kept (it is live). */
  #style: CSSStyleDeclaration | undefined;
  #random: () => number;
  readonly #still: boolean;

  // The wave, the potential and its half step (with the edge), the current's velocity (cells a time unit) and the density.
  readonly #re = new Float64Array(Cells);
  readonly #im = new Float64Array(Cells);
  readonly #v = new Float32Array(Cells);
  readonly #halfRe = new Float64Array(Cells);
  readonly #halfIm = new Float64Array(Cells);
  readonly #fullRe = new Float64Array(Cells);
  readonly #fullIm = new Float64Array(Cells);
  // The emitters' separable sources (a row and a column factor each), reused.
  readonly #rowRe = new Float32Array(N);
  readonly #rowIm = new Float32Array(N);
  readonly #colRe = new Float32Array(N);
  readonly #colIm = new Float32Array(N);
  readonly #vx = new Float32Array(Cells);
  readonly #vy = new Float32Array(Cells);
  readonly #rho = new Float32Array(Cells);
  /** The density's running sum, to draw where a speck starts. */
  readonly #cdf = new Float32Array(inner.length);
  #rhoMax = RhoFloor;
  /** How strong the flow from the left is now (0..1, the glass's assembly, eased). */
  #flow = 0;
  /** Whether the start's waves are in (`#seed`). */
  #seeded = false;
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
  /** Full rate while awake, faded out and no frames at rest (ambient-pace.ts). */
  readonly #ambient: AmbientPace;
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
  /** Where the glass is (`glass-footprint.ts`): the well, and where the light is out of sight. */
  #footprint: GlassFootprint | undefined;
  /** How much of each cell the glass covers now (0..1, from its pieces; `CoverHide`). */
  readonly #cover = new Float32Array(Cells);
  #wasPull = 0;
  #bursts = 0;

  constructor(
    canvas: HTMLCanvasElement,
    options: Partial<QuantumFieldOptions> = {},
    init: { still?: boolean; pace?: AmbientPaceOptions } = {},
  ) {
    this.#canvas = canvas;
    this.#ctx = canvas.getContext('2d');
    this.#options = { ...QuantumFieldDefaults, ...options };
    this.#random = generator(this.#options.seed);
    this.#still = !!init.still;
    this.#ambient = new AmbientPace(init.pace);
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
    // A change (a story's control) is seen: awake.
    this.#awake();
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
      this.#seeded = false;
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
    // A new size: awake (the field fills its new room smoothly; a resized canvas is blank, the word redrawn).
    this.#awake();
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

  /**
   * Draws the last state now, at once: a resized canvas is blank, and the next frame comes after this one is painted
   * (from a resize observer, past the frame's tick), so without it every step of a resize shows one blank frame.
   */
  redraw(): void {
    if (this.#still || !this.#running || !this.#width) return;
    this.#draw(this.#ambient.asleep ? 0 : this.#ambient.level);
  }

  readColors(): void {
    if (typeof getComputedStyle !== 'function') return;
    const el = this.#canvas;
    const before = el.style.color;
    // Read every 1.5 s: the gradient (and the glass's atlas made with it) is made again only when a color changed.
    let changed = !this.#gradient;
    for (let i = 0; i < Colors; i++) {
      el.style.color = ColorRefs[i];
      this.#style ??= getComputedStyle(el);
      const v = this.#style.color;
      if (!v || v === this.#colors[i]) continue;
      this.#colors[i] = v;
      changed = true;
      const m = /rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)/.exec(v);
      if (m) {
        this.#rgb[i * 3] = +m[1];
        this.#rgb[i * 3 + 1] = +m[2];
        this.#rgb[i * 3 + 2] = +m[3];
      }
    }
    el.style.color = before;
    this.#dark = 0.299 * this.#rgb[0] + 0.587 * this.#rgb[1] + 0.114 * this.#rgb[2] > 128;
    if (this.#ctx && changed) {
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
      this.#bursts += 3;
    };
    this.#glass = glass;
    this.#footprint = new GlassFootprint();
    glass.scatter();
    return true;
  }

  shatterAt(x: number, y: number): boolean {
    const glass = this.#glass;
    if (!glass || !this.#options.interactive || this.#still) return false;
    const u = (x - this.#wordX) / this.#wordScale;
    const v = (y - this.#wordY) / this.#wordScale;
    if (!glass.hits(u, v)) return false;
    this.#awake();
    glass.shatter(u, v);
    this.#kickAt(x, y, this.#options.kick * 1.6, 90);
    this.#potential = true;
    return true;
  }

  pointerMove(x: number, y: number, _now: number): void {
    if (!this.#options.interactive || this.#still) return;
    this.#awake();
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
    // Its start at full rate (the word gathers, the wave comes in), until it rests (ambient-pace.ts).
    this.#ambient.input();
    if (this.#still) this.#drawStill();
    else this.#schedule();
  }

  stop(): void {
    this.#running = false;
    if (this.#frame) cancelFrameCallback(this.#frame);
    this.#frame = 0;
  }

  setVisible(visible: boolean): void {
    this.#visible = visible;
    if (visible) this.#schedule();
  }

  dispose(): void {
    this.stop();
  }

  /** The pace's probe: its steps (the rest guard counts them). */
  get steps(): number {
    return this.#ambient.steps;
  }

  /** Whether it has faded out and asks for no frame (a probe; input wakes it). */
  get asleep(): boolean {
    return this.#ambient.asleep;
  }

  /** Input (or a new size): awake, and the frames again if it was asleep. */
  #awake(): void {
    this.#ambient.input();
    this.#schedule();
  }

  #schedule(): void {
    if (this.#still || !this.#running || !this.#visible || this.#frame || this.#ambient.asleep) return;
    this.#last = 0;
    this.#ambient.reset();
    this.#frame = requestFrameCallback(this.#tick);
  }

  readonly #tick = (now: number) => {
    this.#frame = 0;
    if (!this.#running || !this.#visible) return;
    // At rest the wave and the specks fade out and the frames stop; the glass (shattering, gathering, settling) keeps it
    // awake. The dt below is clamped, so the wave resumes where it was rather than catching up.
    const level = this.#ambient.advance(!!this.#glass?.animating);
    if (this.#ambient.asleep) {
      // Faded out: the word alone, the last frame; nothing is asked for until input wakes it.
      this.#draw(0);
      return;
    }
    this.#ambient.stepped();
    const dt = this.#last ? Math.min(0.1, (now - this.#last) / 1000) : 1 / 60;
    this.#last = now;
    if (now - this.#colorsRead > 1500) {
      this.#colorsRead = now;
      this.readColors();
    }
    this.#frameDt = dt;
    this.#step(dt);
    this.#draw(level);
    this.#frame = requestFrameCallback(this.#tick);
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
    this.#seeded = false;
    this.#time = 0;
    this.#random = generator(this.#options.seed);
    for (let t = 0; t < 6; t += 1 / 20) this.#step(1 / 20);
    this.#draw(1);
  }

  /** Time units a second: a packet of the mean wavenumber goes `speed` px a second (its group velocity is k). */
  #rate(): number {
    const o = this.#options;
    return o.speed / (Math.max(0.2, o.k) * Cell);
  }

  /** An emitter at grid (cx, cy) radiating a packet going `angle` with wavenumber `k` (none if all are busy). */
  /** The start's broad, faint waves, all over the grid (inside its absorbing edge). */
  #seed(): void {
    const o = this.#options;
    const re = this.#re;
    const im = this.#im;
    const s2 = 2 * SeedSize * SeedSize;
    for (let n = 0; n < Seeds; n++) {
      // Along the beam's whole length (spread evenly, a little jitter), across its band, going right like the front.
      const sy = (this.#oy - this.#gy) / Cell;
      const cx = Source + ((n + 0.3 + 0.4 * this.#random()) / Seeds) * (N - Edge - Source);
      const cy = sy + (this.#random() * 2 - 1) * FrontHalf * 1.4;
      const angle = (this.#random() * 2 - 1) * SeedSpread;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      const kx = k * Math.cos(angle);
      const ky = k * Math.sin(angle);
      const phase0 = this.#random() * Math.PI * 2;
      for (let y = 0; y < N; y++) {
        const dy = y - cy;
        for (let x = 0; x < N; x++) {
          const dx = x - cx;
          const g = SeedAmplitude * Math.exp(-(dx * dx + dy * dy) / s2);
          if (g < 1e-4) continue;
          const phase = kx * dx + ky * dy + phase0;
          const i = y * N + x;
          re[i] += g * Math.cos(phase);
          im[i] += g * Math.sin(phase);
        }
      }
    }
  }

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
      // Separable: a Gaussian along the axes (long along the way for the front, which goes right; round for a burst
      // along an arm) times the plane wave e^{i(kx·dx + ky·dy)}, so a row factor times a column factor.
      const straight = Math.abs(sn) < 1e-3;
      const sx = straight ? EmitAlong : (EmitAlong + EmitAcross) / 2;
      const sy = straight ? EmitAcross : sx;
      const kx = k * c;
      const ky = k * sn;
      const x0 = Math.max(1, Math.floor(cx - 3 * sx));
      const x1 = Math.min(N - 2, Math.ceil(cx + 3 * sx));
      const y0 = Math.max(1, Math.floor(cy - 3 * sy));
      const y1 = Math.min(N - 2, Math.ceil(cy + 3 * sy));
      const rowRe = this.#rowRe;
      const rowIm = this.#rowIm;
      for (let x = x0; x <= x1; x++) {
        const dx = x - cx;
        const g = Math.exp(-(dx * dx) / (2 * sx * sx));
        rowRe[x] = g * Math.cos(kx * dx);
        rowIm[x] = g * Math.sin(kx * dx);
      }
      // Its strength and phase now: e^{−iωt}, ω = k²/2, folded into the column factor.
      const t = this.#sim - this.#emT0[e];
      const wt = (-k * k * t) / 2;
      const colRe = this.#colRe;
      const colIm = this.#colIm;
      for (let y = y0; y <= y1; y++) {
        const dy = y - cy;
        const g = strength * Math.exp(-(dy * dy) / (2 * sy * sy));
        colRe[y] = g * Math.cos(ky * dy + wt);
        colIm[y] = g * Math.sin(ky * dy + wt);
      }
      for (let y = y0; y <= y1; y++) {
        const pr = colRe[y];
        const pi = colIm[y];
        if (Math.abs(pr) + Math.abs(pi) < 1e-7) continue;
        const o = y * N;
        for (let x = x0; x <= x1; x++) {
          re[o + x] += pr * rowRe[x] - pi * rowIm[x];
          im[o + x] += pr * rowIm[x] + pi * rowRe[x];
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

  /** The potential (the glass, or the star alone without it; the pointer's barrier, the snap back's well) and its half step. */
  #buildPotential(pull: number): void {
    const o = this.#options;
    const v = this.#v;
    v.fill(0);
    const e0 = (o.k * o.k) / 2;
    const sx = (this.#ox - this.#gx) / Cell;
    const sy = (this.#oy - this.#gy) / Cell;
    if (this.#glass && this.#footprint?.ready) {
      // The glass: a well as deep as the glass covers each cell, carried by its pieces.
      this.#footprint.rasterize(
        this.#glass,
        this.#wordX,
        this.#wordY,
        this.#wordScale,
        this.#gx,
        this.#gy,
        Cell,
        N,
        N,
        this.#cover,
      );
      const depth = o.prism * e0;
      for (let i = 0; i < Cells; i++) v[i] -= depth * this.#cover[i];
    } else {
      this.#cover.fill(0);
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
    // Its half step and full step (two halves back to back, between the steps of a frame); most cells have none.
    const hr = this.#halfRe;
    const hi = this.#halfIm;
    const fr = this.#fullRe;
    const fi = this.#fullIm;
    for (let i = 0; i < Cells; i++) {
      if (v[i] === 0) {
        hr[i] = edgeHalf[i];
        hi[i] = 0;
        fr[i] = edgeKeep[i];
        fi[i] = 0;
        continue;
      }
      const phase = (-v[i] * Dt) / 2;
      const keep = edgeHalf[i];
      const c = Math.cos(phase);
      const sn = Math.sin(phase);
      hr[i] = c * keep;
      hi[i] = sn * keep;
      fr[i] = (c * c - sn * sn) * edgeKeep[i];
      fi[i] = 2 * c * sn * edgeKeep[i];
    }
  }

  /** ψ times a potential step (`#halfRe`/`#halfIm` or `#fullRe`/`#fullIm`). */
  #potentialStep(pr: Float64Array, pi: Float64Array): void {
    const re = this.#re;
    const im = this.#im;
    for (let i = 0; i < Cells; i++) {
      const a = re[i];
      const b = im[i];
      re[i] = a * pr[i] - b * pi[i];
      im[i] = a * pi[i] + b * pr[i];
    }
  }

  /**
   * `steps` steps of the evolution (split-step Fourier, Strang): a half step of the potential, then for each step the
   * kinetic step in momentum space and a full potential step (the two halves between steps at once), the last a half.
   * The emitters radiate before each.
   */
  #evolve(steps: number): void {
    const re = this.#re;
    const im = this.#im;
    this.#radiate();
    this.#sim += Dt;
    this.#potentialStep(this.#halfRe, this.#halfIm);
    for (let n = 0; n < steps; n++) {
      forward2(re, im);
      for (let i = 0; i < Cells; i++) {
        const a = re[i];
        const b = im[i];
        re[i] = a * kineticRe[i] - b * kineticIm[i];
        im[i] = a * kineticIm[i] + b * kineticRe[i];
      }
      inverse2(re, im);
      if (n === steps - 1) {
        this.#potentialStep(this.#halfRe, this.#halfIm);
      } else {
        this.#potentialStep(this.#fullRe, this.#fullIm);
        this.#radiate();
        this.#sim += Dt;
      }
    }
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
    const cdf = this.#cdf;
    const cover = this.#cover;
    // Only the vignette's cells: outside it nothing shows, nothing is born and the velocity stays 0.
    for (let n = 0; n < inner.length; n++) {
      const i = inner[n];
      const a = re[i];
      const b = im[i];
      const p = a * a + b * b;
      rho[i] = p;
      // What shows (none under the glass) sets the scale and the births (by |ψ|²).
      const c = cover[i] * CoverHide;
      const shown = c >= 1 ? 0 : p * (1 - c);
      if (shown > most) most = shown;
      sum += shown;
      cdf[n] = sum;
      if (p < 1e-9) {
        vx[i] = 0;
        vy[i] = 0;
        continue;
      }
      // ∇ of the phase from the phase differences to the neighbors (arg ψ(x+1)ψ*(x), and back), averaged: exact for a
      // plane wave, where a central difference of ψ would give sin k (and the colors and speeds with it).
      const r1 = re[i + 1];
      const i1 = im[i + 1];
      const r0 = re[i - 1];
      const i0 = im[i - 1];
      let u = (atan2(i1 * a - r1 * b, r1 * a + i1 * b) + atan2(b * r0 - a * i0, a * r0 + b * i0)) / 2;
      const rd = re[i + N];
      const id = im[i + N];
      const ru = re[i - N];
      const iu = im[i - N];
      let w = (atan2(id * a - rd * b, rd * a + id * b) + atan2(b * ru - a * iu, a * ru + b * iu)) / 2;
      const s2 = u * u + w * w;
      if (s2 > cap2) {
        const k = cap / Math.sqrt(s2);
        u *= k;
        w *= k;
      }
      vx[i] = u;
      vy[i] = w;
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
    if (s < 1e-3 || Math.abs(atan2(w, u)) < Straight) return 0;
    const lo = o.k * (1 - o.spectrum);
    const band = (s - lo) / Math.max(1e-3, 2 * o.k * o.spectrum);
    return 1 + Math.min(3, Math.max(0, Math.floor(band * 4)));
  }

  /** A speck where the wave is (drawn from its density); false if there is next to none. */
  #spawn(i: number): boolean {
    const cdf = this.#cdf;
    const total = cdf[cdf.length - 1];
    if (total < 1e-3) return false;
    const pick = this.#random() * total;
    let lo = 0;
    let hi = cdf.length - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (cdf[mid] < pick) lo = mid + 1;
      else hi = mid;
    }
    const cell = inner[lo];
    const gx = cell % N;
    const gy = (cell / N) | 0;
    this.#x[i] = this.#gx + (gx + this.#random()) * Cell;
    this.#y[i] = this.#gy + (gy + this.#random()) * Cell;
    this.#age[i] = 0;
    this.#shimmer[i] = this.#random() * Math.PI * 2;
    this.#color[i] = this.#colorFor(this.#vx[cell], this.#vy[cell]);
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
    if (!this.#seeded) {
      this.#seeded = true;
      this.#seed();
    }
    this.#emit += dt * o.rate;
    while (this.#emit >= 1) {
      this.#emit -= 1;
      const y = sy + (this.#random() * 2 - 1) * FrontHalf;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      this.#launch(Source, y, 0, k, FrontAmplitude, true);
    }
    // The glass back together: a little of the light it held, out along its arms.
    while (this.#bursts >= 1) {
      this.#bursts -= 1;
      const arm = (Math.floor(this.#random() * 8) * Math.PI) / 4;
      const k = o.k * (1 + o.spectrum * (this.#random() * 2 - 1));
      this.#launch(sx + Math.cos(arm) * 4, sy + Math.sin(arm) * 4, arm, k, BurstAmplitude, false);
    }

    // The evolution: as many steps as the time that passed takes.
    this.#clock += dt * this.#rate();
    let steps = 0;
    while (this.#clock >= Dt && steps < MaxSteps) {
      this.#clock -= Dt;
      steps++;
    }
    if (steps === MaxSteps) this.#clock = 0;
    if (steps) this.#evolve(steps);
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
      // Gone when faded, or out of the vignette (nothing shows there, and the current is not kept there).
      if (
        this.#age[i] > gone ||
        nx < x0 ||
        ny < y0 ||
        nx >= x1 ||
        ny >= y1 ||
        edgeShow[((((ny - y0) / Cell) | 0) * N + (((nx - x0) / Cell) | 0)) | 0] === 0
      ) {
        this.#kill(i);
        i--;
      }
    }
    // Births: |ψ|² is the chance of a speck. As many a second as asked at the wave's usual amount, more with more of it,
    // each drawn from its density.
    const amount = this.#cdf[this.#cdf.length - 1];
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
    const glow = Math.min(1, o.glow) * (this.#dark ? DarkGlow : 1) * 255 * Math.min(1, this.#time / SeedFade);
    const rgb = this.#rgb;
    const tr = this.#tintR;
    const tg = this.#tintG;
    const tb = this.#tintB;
    if (!this.#tinted) {
      // From the middle of the logo's colors (brand-2 and -3); each cell finds its own as the wave comes.
      const m = o.colors ? 2 * 3 : 0;
      const n = o.colors ? 3 * 3 : 0;
      tr.fill((rgb[m] + rgb[n]) / 2);
      tg.fill((rgb[m + 1] + rgb[n + 1]) / 2);
      tb.fill((rgb[m + 2] + rgb[n + 2]) / 2);
      this.#tinted = true;
    }
    const ease = Math.min(1, this.#frameDt * TintEase);
    const lo = o.k * (1 - o.spectrum);
    const span = Math.max(1e-3, 2 * o.k * o.spectrum);
    // Only the vignette's cells (the rest of the image stays clear from the start).
    for (let n = 0; n < inner.length; n++) {
      const i = inner[n];
      const p = rho[i] * inv;
      const j = i * 4;
      if (p < 1e-4) {
        data[j + 3] = 0;
        continue;
      }
      // Its light's color: the logo's, by its wavenumber (brand-1 the longest waves to brand-4 the shortest, blended),
      // followed only as far as the wave is dense (where it is thin its phase is noise: the cell keeps its color).
      let r = rgb[0];
      let g = rgb[1];
      let b = rgb[2];
      let follow = 1;
      if (o.colors) {
        const u = this.#vx[i];
        const w = this.#vy[i];
        const f = Math.min(2.999, Math.max(0, ((Math.sqrt(u * u + w * w) - lo) / span) * 4 - 0.5));
        const k = f | 0;
        const t = f - k;
        const a = (1 + k) * 3;
        const z = (2 + k) * 3;
        r = rgb[a] * (1 - t) + rgb[z] * t;
        g = rgb[a + 1] * (1 - t) + rgb[z + 1] * t;
        b = rgb[a + 2] * (1 - t) + rgb[z + 2] * t;
        follow = Math.min(1, p * 12);
      }
      const e = ease * follow;
      tr[i] += (r - tr[i]) * e;
      tg[i] += (g - tg[i]) * e;
      tb[i] += (b - tb[i]) * e;
      data[j] = tr[i];
      data[j + 1] = tg[i];
      data[j + 2] = tb[i];
      const front = 1 - o.ripples * 0.5 * (1 - this.#re[i] / Math.sqrt(rho[i]));
      const clear = 1 - Math.min(1, this.#cover[i] * CoverHide);
      data[j + 3] = Math.min(255, glow * Math.min(1, Math.sqrt(p)) * front * edgeShow[i] * clear);
    }
    this.#pictureCtx.putImageData(image, 0, 0);
  }

  /** `level`: how much of the wave and the specks shows (0..1; the word is always drawn whole). */
  #draw(level: number): void {
    const ctx = this.#ctx;
    if (!ctx) return;
    const o = this.#options;
    const dpr = this.#dpr;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, this.#canvas.width, this.#canvas.height);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    if (level > 0 && o.glow > 0 && this.#picture) {
      this.#paint();
      ctx.imageSmoothingEnabled = true;
      ctx.globalAlpha = level;
      // On a dark canvas the wave is light added to it (a veil of grey would sink into the black).
      if (this.#dark) ctx.globalCompositeOperation = 'lighter';
      ctx.drawImage(this.#picture, this.#gx, this.#gy, N * Cell, N * Cell);
      ctx.globalCompositeOperation = 'source-over';
      ctx.globalAlpha = 1;
    }

    // The specks, bucketed by color and alpha level (counting sort): one path each.
    const counts = this.#counts;
    counts.fill(0);
    // Halving every `decay` seconds.
    const fade = -Math.LN2 / Math.max(0.1, o.decay);
    for (let i = 0; level > 0 && i < this.#alive; i++) {
      const age = this.#age[i];
      // In over its birth, then halving every `decay` seconds (not the wave under it: it is a particle now).
      const x = this.#x[i];
      const y = this.#y[i];
      const clear = 1 - Math.min(1, this.#at(this.#cover, x, y) * CoverHide);
      const a = Math.min(1, age / Birth) * Math.exp(age * fade) * this.#at(edgeShow, x, y) * clear;
      const lv = Math.max(0, Math.min(Levels - 1, Math.round(a * (Levels - 1))));
      const c = o.colors ? this.#color[i] : 0;
      const b = c * Levels + lv;
      this.#bucketOf[i] = b;
      counts[b]++;
    }
    let sum = 0;
    for (let b = 0; b < counts.length; b++) {
      this.#starts[b] = sum;
      sum += counts[b];
    }
    counts.fill(0);
    for (let i = 0; level > 0 && i < this.#alive; i++) {
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
        ctx.globalAlpha = ((c === 0 ? 0.5 : 1) * o.alpha * l * level) / (Levels - 1);
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
      this.#dpr,
    );
  }
}
