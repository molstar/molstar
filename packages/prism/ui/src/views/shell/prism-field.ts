/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The empty canvas's splash (shell.md, "Appear when relevant"): light through the Mol* star. A sparse, faint
// stream of specks drifts in from the left behind the wordmark (white light, in the canvas's subtle text color) and, where
// it reaches the star, leaves it along the star's arms (the asterisk's eight at 45°, less the one the light comes in
// along), each ray a small spectrum of the logo's colors (brand-1 to brand-4, each bent a little its own way across the
// ray), spreading and fading out. After Stories Untold's thread field (its beads that scatter in and out), without the
// threads and much quieter.
//
// Circles, playful: dots in the stream, bubbles of varied sizes in the rays that swell and shrink and twinkle each at its
// own pace. It moves on its own: the stream rushes in at the start (meeting the wordmark's colors as they wipe in), the
// light comes in waves, the rays breathe. The light goes through air (`fluid.ts`, Navier–Stokes on a coarse grid), stirred
// by a few slow eddies of its own, and each speck moves at its own velocity plus the air's where it is (a classic
// Runge–Kutta step through the flow). And it answers the pointer (mouse, pen or touch) without taking it from the canvas:
// specks near it are pushed aside, its stroke drives the air (so the light swirls round behind it, in eddies that keep
// turning a while), the arms facing it take more of the light, moving it brings more light, and a press releases a burst
// along every arm. A speck pushed off the star's line passes by uncolored.
//
// A 2D canvas: everything per frame is preallocated (memory: render paths stay allocation-free); it draws nothing while
// hidden or stopped, and one still frame under reduced motion.

import { WordmarkGeometry } from '../../components/icon/wordmark.js';
import { Fluid } from './fluid.js';
import { Glass } from './glass.js';
import { GlassFootprint } from './glass-footprint.js';
import { cancelFrameCallback, requestFrameCallback } from '../../reactive/frame-state.js';
import { AmbientPace, type AmbientPaceOptions } from './ambient-pace.js';

export interface PrismFieldOptions {
  /** Specks emitted per second. */
  rate: number;
  /** Their speed, px per second. */
  speed: number;
  /** How far before the split the stream starts (px; it fades in over its first third). */
  beam: number;
  /** The stream's half height at its start (px); it narrows into the split. */
  width: number;
  /** A ray's half angle (degrees): across it the first color bends one way by it, the last the other. */
  spread: number;
  /** How far past the split a colored speck goes before it is gone (px). */
  reach: number;
  /** The most a speck shows (0..1); the stream shows at a third of it. */
  alpha: number;
  /** A speck's radius (px): a dot of the stream; a colored one is a little larger. */
  size: number;
  /** How much their sizes vary (0: all alike; 1: now and then a bubble three times as large). */
  variety: number;
  /** The colored ones swell and shrink as they go (0..1), each at its own pace. */
  bob: number;
  /** How much they twinkle (0..1). */
  twinkle: number;
  /** The rays in the logo's colors; in the stream's color when off. */
  colors: boolean;
  /** The light comes in waves (0: steady, 1: from nothing to twice the rate). */
  pulse: number;
  /**
   * Seconds the stream rushes in at the start (faster and denser, easing to its pace), so the first light reaches the star
   * as the wordmark fades in; 0: it starts at its pace.
   */
  intro: number;
  /** The air stirs on its own: slow eddies the light drifts in (0: still air). */
  curl: number;
  /** The rays breathe: their spread swells and narrows by this much (0..1). */
  breathe: number;
  /** It answers the pointer. */
  interactive: boolean;
  /** How far the pointer reaches (px). */
  radius: number;
  /** How hard it pushes specks aside (px/s² at its center). */
  push: number;
  /** How much its stroke drives the air (0..1): the light swirls on behind it. */
  wake: number;
  /** How long the air's eddies keep turning (0: they die down soon; 1: they hold their spin). */
  swirl: number;
  /** How much more light goes out along the arms facing it (0..1). */
  follow: number;
  /** Specks a press releases from the split. */
  burst: number;
  seed: number;
}

export const PrismFieldDefaults: PrismFieldOptions = {
  rate: 24,
  speed: 26,
  beam: 260,
  width: 10,
  spread: 7,
  reach: 320,
  alpha: 0.55,
  size: 1.4,
  variety: 0.7,
  bob: 0.35,
  twinkle: 0.45,
  colors: true,
  pulse: 0.55,
  intro: 2.4,
  curl: 0.6,
  breathe: 0.3,
  interactive: true,
  radius: 90,
  push: 520,
  wake: 0.5,
  swirl: 0.5,
  follow: 0.5,
  burst: 36,
  seed: 3,
};

/** Specks alive at most. */
const Capacity = 1536;
/** Alpha levels a frame batches specks into (per color), so a frame fills a few paths. */
const Levels = 12;
/** The stream's color and the four of the logo. */
const Colors = 5;
/** Seconds a still frame (reduced motion) simulates first, so the rays are full. */
const StillSeconds = 30;
/** Half the split's height (px): a speck crossing its line farther from its point passes by. */
const SplitHalf = 22;
/** How far an uncolored speck goes past the split's line before it is gone (px). */
const PassBy = 140;
/** Seconds the warmth of a stroke lasts, and px of stroke to full warmth. */
const WarmFade = 1.6;
const WarmTravel = 260;
/**
 * The star's arms (radians, canvas y down): the asterisk's eight at 45°, less the left one, which runs into the l and
 * along which the light comes in. The light leaves the star along them.
 */
const Arms = [0, 1, 2, 3, 5, 6, 7].map((k) => (k * Math.PI) / 4);
/** An arm's weight toward the pointer (cosine of its angle to it, raised: the facing arms take most). */
const armWeights = new Float32Array(Arms.length);

/** How far a shatter's shockwave reaches (px), and how hard it throws the specks (px/s at its middle). */
const ShockRadius = 220;
const ShockSpeed = 260;

/**
 * The glass's gravity as it snaps back: how far it reaches (px from the word's middle) and how hard it pulls (px/s² at
 * the middle, at full pull).
 */
const GravityRadius = 360;
const Gravity = 900;

/** The air: px past the light's reach it covers, how fast its motion dies down (1/s), its confinement at full `swirl`. */
const AirMargin = 60;
const AirDrag = 0.5;
const AirConfinement = 5;
/** The fastest stroke the air takes on (px/s). */
const StrokeCap = 700;
/** The air's own eddies: how many, how wide (px), how hard each turns it at full `curl` (px/s²). */
const Eddies = 3;
const EddyRadius = 150;
const EddyForce = 30;
/** A shatter's gust: two eddies turning opposite ways on either side of the press (px apart, px wide, px/s). */
const GustApart = 36;
const GustRadius = 110;
const GustSpeed = 160;

/**
 * The light out of sight under the glass (`glass-footprint.ts`), so the word stays clear, whole or in pieces: a grid of
 * `CoverCell` px round the word, `CoverMargin` px past it each way (room for the flying pieces); a speck under glass
 * shows `1 − cover × CoverHide` of itself.
 */
const CoverCell = 4;
const CoverMargin = 200;
const CoverHide = 1.6;

/** How much faster (and denser) the stream is at the very start of its intro. */
const IntroBoost = 8;
/** Seconds a speck's speed takes to settle back after a push or a burst. */
const Settle = 0.8;

/** The theme's colors the field draws in, as CSS custom properties on (or above) its canvas. */
const ColorVars = [
  '--msp-color-text-subtle',
  '--msp-color-brand-1',
  '--msp-color-brand-2',
  '--msp-color-brand-3',
  '--msp-color-brand-4',
];
/** The same as `var(…)` values, made once (the colors are read every 1.5 s). */
const ColorRefs = ColorVars.map((v) => `var(${v})`);

/** A small seeded generator (mulberry32): the same field for the same seed. */
function generator(seed: number): () => number {
  let a = seed >>> 0 || 1;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export class PrismField {
  #options: PrismFieldOptions;
  readonly #canvas: HTMLCanvasElement;
  readonly #ctx: CanvasRenderingContext2D | null;
  /** The canvas's computed style, kept (it is live). */
  #style: CSSStyleDeclaration | undefined;
  #random: () => number;

  // The specks: position (px, canvas CSS pixels), velocity (px/s), the speed they settle back to, where they split
  // (px), color (0: the stream; 1–4: the logo's), a shimmer phase; `alive` packs the live ones first.
  readonly #x = new Float32Array(Capacity);
  readonly #y = new Float32Array(Capacity);
  readonly #vx = new Float32Array(Capacity);
  readonly #vy = new Float32Array(Capacity);
  readonly #base = new Float32Array(Capacity);
  readonly #sx = new Float32Array(Capacity);
  readonly #sy = new Float32Array(Capacity);
  readonly #shimmer = new Float32Array(Capacity);
  /** Its size (a factor of `size`) and its own pace (twinkle, bob, curl). */
  readonly #scale = new Float32Array(Capacity);
  readonly #pace = new Float32Array(Capacity);
  readonly #color = new Uint8Array(Capacity);
  #alive = 0;
  // A frame's buckets: per color and alpha level, the specks in it.
  readonly #bucket = new Uint16Array(Capacity);
  readonly #bucketOf = new Uint16Array(Capacity);
  readonly #counts = new Uint16Array(Colors * Levels);
  readonly #starts = new Uint16Array(Colors * Levels);

  readonly #colors: string[] = ['#888', '#7a4f2f', '#e0ae2c', '#6f9a80', '#3f86ad'];
  #width = 0;
  #height = 0;
  #dpr = 1;
  /** The split's point (canvas CSS pixels): where the stream splits. */
  #ox = 0;
  #oy = 0;
  #emit = 0;
  #time = 0;
  #last = 0;
  #frame = 0;
  /** Full rate while awake, faded out and no frames at rest (ambient-pace.ts). */
  readonly #ambient: AmbientPace;
  #running = false;
  #visible = true;
  #colorsRead = Number.NEGATIVE_INFINITY;
  readonly #still: boolean;

  // The pointer: where it is (canvas px), its stroke's velocity (px/s, smoothed), whether it is over the field, when it
  // last moved (ms), the warmth its strokes left, a burst waiting.
  #px = 0;
  #py = 0;
  #pvx = 0;
  #pvy = 0;
  #pointer = false;
  #moved = 0;
  #warmth = 0;
  /** Where the pointer is from the star (radians), and how much it draws the light (0..1, eased). */
  #toward = 0;
  #pull = 0;

  // The wordmark as glass (`glass.ts`): where its drawing origin is (canvas px), px a unit, the logo's gradient for its
  // line (in drawing units, made with the colors).
  #glass: Glass | undefined;
  #wordX = 0;
  #wordY = 0;
  #wordScale = 1;
  #gradient: CanvasGradient | undefined;
  #burst = 0;
  // Where the glass is: its footprint on the cover grid (from (coverX, coverY), canvas px), redrawn while it moves.
  #footprint: GlassFootprint | undefined;
  #cover = new Float32Array(0);
  #coverX = 0;
  #coverY = 0;
  #coverW = 0;
  #coverH = 0;
  #coverDirty = true;

  /** The air the light goes through (none when still), and its eddies' phases. */
  readonly #air: Fluid | undefined;
  readonly #eddy = new Float32Array(Eddies * 2);

  constructor(
    canvas: HTMLCanvasElement,
    options: Partial<PrismFieldOptions> = {},
    init: { still?: boolean; pace?: AmbientPaceOptions } = {},
  ) {
    this.#canvas = canvas;
    this.#ctx = canvas.getContext('2d');
    this.#options = { ...PrismFieldDefaults, ...options };
    this.#random = generator(this.#options.seed);
    this.#still = !!init.still;
    this.#ambient = new AmbientPace(init.pace);
    if (!this.#still) {
      this.#air = new Fluid();
      for (let i = 0; i < this.#eddy.length; i++) this.#eddy[i] = this.#random() * Math.PI * 2;
    }
  }

  setOptions(options: Partial<PrismFieldOptions>): void {
    // A change (a story's control) is seen: awake.
    this.#awake();
    const seed = this.#options.seed;
    this.#options = { ...this.#options, ...options };
    if (this.#options.seed !== seed) {
      this.#random = generator(this.#options.seed);
      this.#alive = 0;
    }
    this.#layoutAir();
    if (this.#still) this.#drawStill();
  }

  /** Where the stream splits, in the canvas's CSS pixels (the star's middle). */
  setOrigin(x: number, y: number): void {
    this.#ox = x;
    this.#oy = y;
    this.#layoutAir();
    if (this.#still) this.#drawStill();
  }

  /** The air covers where the light goes: the stream's start to the rays' reach round the star, within the canvas. */
  #layoutAir(): void {
    const air = this.#air;
    if (!air) return;
    const o = this.#options;
    const x0 = Math.max(0, this.#ox - Math.max(o.beam, o.reach) - AirMargin);
    const x1 = Math.min(this.#width, this.#ox + o.reach + AirMargin);
    const y0 = Math.max(0, this.#oy - o.reach - AirMargin);
    const y1 = Math.min(this.#height, this.#oy + o.reach + AirMargin);
    air.layout(x0, y0, Math.max(0, x1 - x0), Math.max(0, y1 - y0));
  }

  /** Reads the canvas's size (and the screen's pixel ratio) and the theme's colors. */
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
    this.#layoutAir();
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

  /** The theme's colors (it may change under the canvas: the chrome follows the canvas background). */
  readColors(): void {
    if (typeof getComputedStyle !== 'function') return;
    // Through the element's `color`, so a `color-mix()` token comes back as a color a canvas takes.
    const el = this.#canvas;
    const before = el.style.color;
    // Read every 1.5 s: the gradient (and the glass's atlas made with it) is made again only when a color changed.
    let changed = !this.#gradient;
    for (let i = 0; i < Colors; i++) {
      el.style.color = ColorRefs[i];
      this.#style ??= getComputedStyle(el);
      const v = this.#style.color;
      if (v && v !== this.#colors[i]) {
        this.#colors[i] = v;
        changed = true;
      }
    }
    el.style.color = before;
    // The wordmark's gradient, as its outline has it (wordmark.tsx), in drawing units.
    if (this.#ctx && changed) {
      const g = this.#ctx.createLinearGradient(0, 0, WordmarkGeometry.width - 2, 0);
      g.addColorStop(0, this.#colors[1]);
      g.addColorStop(0.35, this.#colors[2]);
      g.addColorStop(0.6, this.#colors[3]);
      g.addColorStop(1, this.#colors[4]);
      this.#gradient = g;
    }
  }

  /**
   * The wordmark the light goes through, as glass: its drawing origin (canvas CSS px) and px a unit. The first time, it
   * starts in pieces and assembles. True if the canvas draws the word (the page's wordmark then only places it); false
   * when still (reduced motion) or where there is nothing to draw it with (Node, jsdom).
   */
  setWord(x: number, y: number, scale: number): boolean {
    this.#wordX = x;
    this.#wordY = y;
    this.#wordScale = scale;
    this.#layoutCover();
    if (this.#glass) return true;
    if (this.#still || typeof Path2D !== 'function') return false;
    const glass = new Glass(this.#options.seed + 8);
    // Back together: the star splits the light again, with a pop of it.
    glass.onSnap = () => {
      this.#burst += this.#options.burst * 0.6;
    };
    this.#glass = glass;
    this.#footprint = new GlassFootprint();
    glass.scatter();
    return true;
  }

  /** The cover grid round the word where it is now (grown when it must, else reused). */
  #layoutCover(): void {
    const w = Math.ceil((WordmarkGeometry.width * this.#wordScale + 2 * CoverMargin) / CoverCell);
    const h = Math.ceil((WordmarkGeometry.height * this.#wordScale + 2 * CoverMargin) / CoverCell);
    if (this.#cover.length < w * h) this.#cover = new Float32Array(w * h);
    this.#coverX = this.#wordX - CoverMargin;
    this.#coverY = this.#wordY - CoverMargin;
    this.#coverW = w;
    this.#coverH = h;
    this.#coverDirty = true;
  }

  /** How much of a speck at (x, y) shows (0 under the glass, 1 clear of it). */
  #clear(x: number, y: number): number {
    const cx = (x - this.#coverX) / CoverCell;
    const cy = (y - this.#coverY) / CoverCell;
    if (cx < 0 || cy < 0 || cx >= this.#coverW || cy >= this.#coverH) return 1;
    return 1 - Math.min(1, this.#cover[(cy | 0) * this.#coverW + (cx | 0)] * CoverHide);
  }

  /**
   * A press at (x, y) (canvas CSS px) on the wordmark: it shatters, the light it held thrown out in a shockwave and a
   * burst; true if it was on the word.
   */
  shatterAt(x: number, y: number): boolean {
    const glass = this.#glass;
    if (!glass || !this.#options.interactive || this.#still) return false;
    const u = (x - this.#wordX) / this.#wordScale;
    const v = (y - this.#wordY) / this.#wordScale;
    // Whole, or already in pieces (shattered further from where they are).
    if (!glass.hits(u, v)) return false;
    this.#awake();
    glass.shatter(u, v);
    const r2 = ShockRadius * ShockRadius;
    for (let i = 0; i < this.#alive; i++) {
      const dx = this.#x[i] - x;
      const dy = this.#y[i] - y;
      const d2 = dx * dx + dy * dy;
      if (d2 >= r2) continue;
      const d = Math.sqrt(d2) + 1e-3;
      const f = 1 - d / ShockRadius;
      this.#vx[i] += (dx / d) * ShockSpeed * f;
      this.#vy[i] += (dy / d) * ShockSpeed * f;
    }
    // The shards stir the air: a gust out of the press, two eddies either side of it turning opposite ways.
    const air = this.#air;
    if (air) {
      const a = this.#random() * Math.PI * 2;
      const ax = Math.cos(a) * GustApart;
      const ay = Math.sin(a) * GustApart;
      air.swirl(x + ax, y + ay, GustSpeed, GustRadius);
      air.swirl(x - ax, y - ay, -GustSpeed, GustRadius);
    }
    this.#burst += this.#options.burst;
    return true;
  }

  /** The pointer moved, at (x, y) in the canvas's CSS pixels; `now` in ms. */
  pointerMove(x: number, y: number, now: number): void {
    if (!this.#options.interactive || this.#still) return;
    this.#awake();
    if (this.#pointer && this.#moved) {
      const dt = Math.max(0.008, (now - this.#moved) / 1000);
      // Smoothed: a stroke, not a jitter.
      this.#pvx += ((x - this.#px) / dt - this.#pvx) * 0.35;
      this.#pvy += ((y - this.#py) / dt - this.#pvy) * 0.35;
      this.#warmth = Math.min(1, this.#warmth + Math.hypot(x - this.#px, y - this.#py) / WarmTravel);
    }
    this.#px = x;
    this.#py = y;
    this.#moved = now;
    this.#pointer = true;
    this.#glass?.setPointer((x - this.#wordX) / this.#wordScale, (y - this.#wordY) / this.#wordScale);
  }

  /** The pointer left (or a touch ended). */
  pointerLeave(): void {
    this.#pointer = false;
    this.#glass?.clearPointer();
    this.#pvx = 0;
    this.#pvy = 0;
  }

  /** A press (a click, a tap): a burst of color from the split. */
  pointerDown(x: number, y: number, now: number): void {
    if (!this.#options.interactive || this.#still) return;
    this.pointerMove(x, y, now);
    this.#burst += this.#options.burst;
  }

  /** Runs while visible; nothing is drawn or scheduled while stopped or hidden. */
  start(): void {
    this.#running = true;
    // Its start at full rate (the word gathers, the light comes in), until it rests (ambient-pace.ts).
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
    // At rest the particles fade out and the frames stop; the glass (shattering, gathering, settling) keeps it awake.
    const level = this.#ambient.advance(!!this.#glass?.animating);
    if (this.#ambient.asleep) {
      // Faded out: the word alone, the last frame; nothing is asked for until input wakes it.
      this.#draw(0);
      return;
    }
    this.#ambient.stepped();
    // A long pause (a hidden tab) resumes where it was rather than catching up.
    const dt = this.#last ? Math.min(0.1, (now - this.#last) / 1000) : 1 / 60;
    this.#last = now;
    if (now - this.#colorsRead > 1500) {
      this.#colorsRead = now;
      this.readColors();
    }
    // A pointer that stopped moving has no stroke.
    if (now - this.#moved > 80) {
      this.#pvx *= 0.85;
      this.#pvy *= 0.85;
    }
    this.#step(dt);
    this.#draw(level);
    this.#frame = requestFrameCallback(this.#tick);
  };

  #drawStill(): void {
    if (!this.#width) return;
    this.#alive = 0;
    this.#time = 0;
    this.#random = generator(this.#options.seed);
    for (let t = 0; t < StillSeconds; t += 1 / 30) this.#step(1 / 30);
    this.#draw(1);
  }

  /** A speck of the stream, at its start, toward the split's point. */
  #spawnStream(): void {
    if (this.#alive >= Capacity) return;
    const o = this.#options;
    const i = this.#alive++;
    this.#x[i] = this.#ox - o.beam;
    this.#y[i] = this.#oy + (this.#random() * 2 - 1) * o.width;
    const dy = this.#oy - this.#y[i];
    const len = Math.hypot(o.beam, dy) || 1;
    const v = o.speed * (0.8 + this.#random() * 0.4);
    this.#vx[i] = (v * o.beam) / len;
    this.#vy[i] = (v * dy) / len;
    this.#base[i] = v;
    this.#color[i] = 0;
    this.#shimmer[i] = this.#random() * Math.PI * 2;
    // A speck of the stream is a small dot; it gets its own size when it splits.
    this.#scale[i] = 0.55 + this.#random() * 0.35;
    this.#pace[i] = 0.6 + this.#random() * 0.9;
  }

  /** A speck of a burst: colored at the star, fast at first, along an arm and a little wider. */
  #spawnBurst(spread: number): void {
    if (this.#alive >= Capacity) return;
    const o = this.#options;
    const i = this.#alive++;
    const c = 1 + Math.floor(this.#random() * 4);
    const angle = this.#fanAngle(c, spread * 1.6);
    const v = o.speed * (2.2 + this.#random() * 1.6);
    this.#x[i] = this.#ox;
    this.#sx[i] = this.#ox;
    this.#y[i] = this.#oy;
    this.#sy[i] = this.#oy;
    this.#vx[i] = v * Math.cos(angle);
    this.#vy[i] = v * Math.sin(angle);
    this.#base[i] = o.speed * (0.8 + this.#random() * 0.4);
    this.#color[i] = c;
    this.#shimmer[i] = this.#random() * Math.PI * 2;
    this.#scale[i] = this.#bubble();
    this.#pace[i] = 0.6 + this.#random() * 0.9;
  }

  /** A colored speck's size: most about the same, now and then a bubble (up to three times as large). */
  #bubble(): number {
    const r = this.#random();
    return 1 + this.#options.variety * 2 * r * r * r;
  }

  /**
   * The way a speck of color `c` goes: along one of the star's arms (those facing the pointer more likely), bent across
   * the arm by its color (a small spectrum on each), a little scattered.
   */
  #fanAngle(c: number, spread: number): number {
    let total = 0;
    for (let k = 0; k < Arms.length; k++) {
      const facing = Math.max(0, Math.cos(Arms[k] - this.#toward));
      armWeights[k] = 1 + this.#pull * 6 * facing * facing * facing;
      total += armWeights[k];
    }
    let pick = this.#random() * total;
    let arm = Arms[Arms.length - 1];
    for (let k = 0; k < Arms.length; k++) {
      pick -= armWeights[k];
      if (pick <= 0) {
        arm = Arms[k];
        break;
      }
    }
    return arm + spread * ((c - 1) / 1.5 - 1) + (this.#random() * 2 - 1) * spread * 0.22;
  }

  #kill(i: number): void {
    const last = --this.#alive;
    this.#x[i] = this.#x[last];
    this.#y[i] = this.#y[last];
    this.#vx[i] = this.#vx[last];
    this.#vy[i] = this.#vy[last];
    this.#base[i] = this.#base[last];
    this.#sx[i] = this.#sx[last];
    this.#sy[i] = this.#sy[last];
    this.#color[i] = this.#color[last];
    this.#shimmer[i] = this.#shimmer[last];
    this.#scale[i] = this.#scale[last];
    this.#pace[i] = this.#pace[last];
  }

  #step(dt: number): void {
    const o = this.#options;
    this.#time += dt;
    const t = this.#time;
    this.#warmth = Math.max(0, this.#warmth - dt / WarmFade);
    this.#glass?.step(dt);
    // Where the glass is, while its pieces move (or the word moved).
    const glass = this.#glass;
    if (glass && this.#footprint?.ready && this.#coverW && (this.#coverDirty || !glass.resting)) {
      this.#footprint.rasterize(
        glass,
        this.#wordX,
        this.#wordY,
        this.#wordScale,
        this.#coverX,
        this.#coverY,
        CoverCell,
        this.#coverW,
        this.#coverH,
        this.#cover,
      );
      this.#coverDirty = false;
    }
    // The star splits the light only while the glass is whole.
    const prism = !this.#glass || this.#glass.whole;
    // The intro: the stream rushes in at full speed for its first half, then eases (smoothly) to its pace.
    const u = o.intro > 0 ? Math.min(1, Math.max(0, (t / o.intro - 0.5) * 2)) : 1;
    const boost = 1 + (IntroBoost - 1) * (1 - u * u * (3 - 2 * u));

    // The rays: breathing; the pointer draws the light into the arms facing it (eased in and out).
    const spread = ((o.spread * Math.PI) / 180) * (1 + o.breathe * Math.sin(t * 0.37));
    const drawn = this.#pointer && o.interactive ? o.follow : 0;
    if (drawn) this.#toward = Math.atan2(this.#py - this.#oy, this.#px - this.#ox);
    this.#pull += (drawn - this.#pull) * Math.min(1, dt * 1.5);

    // The light in waves (two slow ones beating), more while the pointer warms it.
    const wave = 0.5 + 0.35 * Math.sin(t * 0.83) + 0.15 * Math.sin(t * 2.1 + 1.3);
    this.#emit += dt * o.rate * (1 - o.pulse + o.pulse * 2 * wave) * (1 + 2 * this.#warmth) * boost;
    while (this.#emit >= 1) {
      this.#emit -= 1;
      this.#spawnStream();
    }
    while (this.#burst >= 1) {
      this.#burst -= 1;
      this.#spawnBurst(spread);
    }

    const pointer = this.#pointer && o.interactive;
    const air = this.#air;
    if (air) {
      // The air's own eddies wander slowly round the star, turning alternate ways.
      if (o.curl > 0) {
        for (let e = 0; e < Eddies; e++) {
          const ex = this.#ox + o.reach * 0.55 * Math.sin(t * (0.05 + 0.03 * e) + this.#eddy[2 * e]);
          const ey = this.#oy + o.reach * 0.45 * Math.cos(t * (0.04 + 0.025 * e) + this.#eddy[2 * e + 1]);
          air.swirl(ex, ey, (e % 2 ? -1 : 1) * o.curl * EddyForce * dt, EddyRadius);
        }
      }
      // A stroke drives the air under it along.
      if (pointer && o.wake > 0 && (this.#pvx || this.#pvy)) {
        // A flick no harder than a brisk stroke.
        const k = o.wake * Math.min(1, StrokeCap / (Math.sqrt(this.#pvx * this.#pvx + this.#pvy * this.#pvy) || 1));
        air.drive(this.#px, this.#py, this.#pvx * k, this.#pvy * k, o.radius * 0.7, Math.min(1, dt * 12));
      }
      air.drag = AirDrag;
      air.confinement = o.swirl * AirConfinement;
      air.step(dt);
    }

    // The glass snapping back pulls the light toward the word's middle for a moment (its pop throws it out again).
    const pull = this.#glass ? this.#glass.pull : 0;
    const gx = this.#wordX + (WordmarkGeometry.width / 2 - 1) * this.#wordScale;
    const gy = this.#wordY + (WordmarkGeometry.height / 2 - 1) * this.#wordScale;
    const g2 = GravityRadius * GravityRadius;
    const r = o.radius;
    const r2 = r * r;
    const settle = Math.min(1, dt / Settle);
    for (let i = 0; i < this.#alive; i++) {
      let vx = this.#vx[i];
      let vy = this.#vy[i];
      const x = this.#x[i];
      const y = this.#y[i];

      if (pointer) {
        const dx = x - this.#px;
        const dy = y - this.#py;
        const d2 = dx * dx + dy * dy;
        if (d2 < r2) {
          const d = Math.sqrt(d2) + 1e-3;
          const f = 1 - d / r;
          const f2 = f * f;
          // Aside (its stroke moves the air, which carries them).
          vx += (dx / d) * o.push * f2 * dt;
          vy += (dy / d) * o.push * f2 * dt;
        }
      }

      // The speed settles back to the speck's own (its direction stays where it was pushed).
      const v = Math.sqrt(vx * vx + vy * vy) || 1;
      const k = 1 + (this.#base[i] / v - 1) * settle;
      vx *= k;
      vy *= k;

      if (pull > 0) {
        const dx = gx - x;
        const dy = gy - y;
        const d2 = dx * dx + dy * dy;
        if (d2 < g2) {
          const d = Math.sqrt(d2) + 1e-3;
          const f = Gravity * pull * (1 - d / GravityRadius) * dt;
          vx += (dx / d) * f;
          vy += (dy / d) * f;
        }
      }

      // On through the air: its own velocity (the stream at the intro's pace) plus the air's where it is, a classic
      // Runge–Kutta step along the flow.
      const pace = this.#color[i] === 0 ? boost : 1;
      const ox = vx * pace;
      const oy = vy * pace;
      let nx = x + ox * dt;
      let ny = y + oy * dt;
      if (air?.on) {
        const h = dt / 2;
        air.sample(x, y);
        const k1x = air.su;
        const k1y = air.sv;
        air.sample(x + (ox + k1x) * h, y + (oy + k1y) * h);
        const k2x = air.su;
        const k2y = air.sv;
        air.sample(x + (ox + k2x) * h, y + (oy + k2y) * h);
        const k3x = air.su;
        const k3y = air.sv;
        air.sample(x + (ox + k3x) * dt, y + (oy + k3y) * dt);
        nx += ((k1x + 2 * k2x + 2 * k3x + air.su) * dt) / 6;
        ny += ((k1y + 2 * k2y + 2 * k3y + air.sv) * dt) / 6;
      }
      this.#vx[i] = vx;
      this.#vy[i] = vy;
      this.#x[i] = nx;
      this.#y[i] = ny;

      if (this.#color[i] === 0) {
        if (prism && x < this.#ox && nx >= this.#ox && Math.abs(ny - this.#oy) <= SplitHalf) {
          // Through the star: one of the logo's colors, bent its own way, split where it crossed.
          const c = 1 + Math.floor(this.#random() * 4);
          this.#color[i] = c;
          const angle = this.#fanAngle(c, spread);
          // In the intro it leaves the star with some of the rush, and settles to its pace: the first colors pop out.
          const speed = Math.hypot(vx, vy) * Math.min(3, boost);
          this.#vx[i] = speed * Math.cos(angle);
          this.#vy[i] = speed * Math.sin(angle);
          this.#sx[i] = nx;
          this.#sy[i] = ny;
          this.#scale[i] = this.#bubble();
        } else if (nx > this.#ox + PassBy || nx < this.#ox - o.beam * 1.5 || Math.abs(ny - this.#oy) > o.reach) {
          // Pushed off the split's line: it passes by, and goes.
          this.#kill(i);
          i--;
        }
      } else {
        const dx = nx - this.#sx[i];
        const dy = ny - this.#sy[i];
        if (dx * dx + dy * dy > o.reach * o.reach) {
          this.#kill(i);
          i--;
        }
      }
    }
  }

  /** A speck's alpha level (0: not drawn). */
  #level(i: number): number {
    const o = this.#options;
    let a: number;
    if (this.#color[i] === 0) {
      // The stream: in over its first third, at a third of the most; past the split's line (it missed), out.
      const x = this.#x[i];
      const along = 1 - (this.#ox - x) / o.beam;
      a = (o.alpha / 3) * Math.min(1, Math.max(0, along * 3));
      if (x > this.#ox) a *= Math.max(0, 1 - (x - this.#ox) / PassBy);
    } else {
      const d = Math.hypot(this.#x[i] - this.#sx[i], this.#y[i] - this.#sy[i]) / o.reach;
      // Brightest just past the split, gone at the reach.
      a = o.alpha * Math.min(1, d * 12) * (1 - d) ** 1.6;
    }
    // Out of sight under the glass.
    a *= this.#clear(this.#x[i], this.#y[i]);
    // A twinkle, each at its own pace.
    a *=
      1 - o.twinkle + o.twinkle * (0.5 + 0.5 * Math.sin(this.#time * (1.4 + 2.6 * this.#pace[i]) + this.#shimmer[i]));
    return Math.max(0, Math.min(Levels - 1, Math.round((a / Math.max(o.alpha, 1e-6)) * (Levels - 1))));
  }

  /** `level`: how much of the light shows (0..1; the word is always drawn whole). */
  #draw(level: number): void {
    const ctx = this.#ctx;
    if (!ctx) return;
    const o = this.#options;
    const dpr = this.#dpr;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, this.#canvas.width, this.#canvas.height);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    // Bucket by color and level (counting sort), so each bucket is one path.
    const counts = this.#counts;
    counts.fill(0);
    for (let i = 0; level > 0 && i < this.#alive; i++) {
      const c = o.colors ? this.#color[i] : 0;
      const b = c * Levels + this.#level(i);
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

    // Circles: a dot of the stream, a colored one larger and swelling and shrinking at its own pace.
    const t = this.#time;
    for (let c = 0; c < Colors; c++) {
      ctx.fillStyle = this.#colors[c];
      for (let l = 1; l < Levels; l++) {
        const b = c * Levels + l;
        const n = counts[b];
        if (!n) continue;
        ctx.globalAlpha = (o.alpha * l * level) / (Levels - 1);
        ctx.beginPath();
        const s = this.#starts[b];
        for (let k = 0; k < n; k++) {
          const i = this.#bucket[s + k];
          const x = this.#x[i];
          const y = this.#y[i];
          let r = o.size * this.#scale[i];
          if (this.#color[i] !== 0) r *= 1 + o.bob * Math.sin(t * (1.5 + 2 * this.#pace[i]) + this.#shimmer[i]);
          if (r < 0.2) continue;
          ctx.moveTo(x + r, y);
          ctx.arc(x, y, r, 0, Math.PI * 2);
        }
        ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
    // The wordmark's glass, over the light.
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
