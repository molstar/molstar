/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The empty canvas's splash (design/shell.md, "Appear when relevant"): light through Prism's mark. A sparse, faint stream
// of specks drifts in from the left behind the wordmark (white light, in the canvas's subtle text color); where it
// reaches the mark's prism it splits into a slow fan of the logo's colors (brand-1 to brand-4, as a spectrum: each color
// bends its own way), which spreads and fades out. After Stories Untold's thread field (its beads that scatter in and
// out), without the threads and much quieter. A 2D canvas: everything per frame is preallocated (memory: render paths
// stay allocation-free); it draws nothing while hidden or stopped, and one still frame under reduced motion.

export interface PrismFieldOptions {
  /** Specks emitted per second. */
  rate: number;
  /** Their speed, px per second. */
  speed: number;
  /** How far before the prism the stream starts (px; it fades in over its first third). */
  beam: number;
  /** The stream's half height at its start (px); it narrows into the prism. */
  width: number;
  /** The fan's half angle (degrees): the first color bends up by it, the last down. */
  spread: number;
  /** How far past the prism a colored speck goes before it is gone (px). */
  reach: number;
  /** The most a speck shows (0..1); the stream shows at a third of it. */
  alpha: number;
  /** A speck's radius (px). */
  size: number;
  /** A colored speck's streak behind it (px; 0: a dot). */
  tail: number;
  /** The fan in the logo's colors; in the stream's color when off. */
  colors: boolean;
  seed: number;
}

export const PrismFieldDefaults: PrismFieldOptions = {
  rate: 9,
  speed: 22,
  beam: 260,
  width: 10,
  spread: 16,
  reach: 300,
  alpha: 0.55,
  size: 1,
  tail: 5,
  colors: true,
  seed: 3,
};

/** Specks alive at most. */
const Capacity = 1024;
/** Alpha levels a frame batches specks into (per color), so a frame strokes a few paths. */
const Levels = 12;
/** The stream's color and the four of the logo. */
const Colors = 5;
/** Seconds a still frame (reduced motion) simulates first, so the fan is full. */
const StillSeconds = 30;

/** The theme's colors the field draws in, as CSS custom properties on (or above) its canvas. */
const ColorVars = [
  '--msp-color-text-subtle',
  '--msp-color-brand-1',
  '--msp-color-brand-2',
  '--msp-color-brand-3',
  '--msp-color-brand-4',
];

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
  #random: () => number;

  // The specks: position (px, canvas CSS pixels), velocity (px/s), age (s), color (0: the stream; 1–4: the logo's),
  // a shimmer phase; `alive` packs the live ones first.
  readonly #x = new Float32Array(Capacity);
  readonly #y = new Float32Array(Capacity);
  readonly #vx = new Float32Array(Capacity);
  readonly #vy = new Float32Array(Capacity);
  readonly #age = new Float32Array(Capacity);
  readonly #shimmer = new Float32Array(Capacity);
  readonly #color = new Uint8Array(Capacity);
  #alive = 0;
  // A frame's buckets: per color and alpha level, the specks in it.
  readonly #bucket = new Uint16Array(Capacity);
  readonly #bucketOf = new Uint16Array(Capacity);
  readonly #counts = new Uint16Array(Colors * Levels);
  readonly #starts = new Uint16Array(Colors * Levels);

  #colors: string[] = ['#888', '#7a4f2f', '#e0ae2c', '#6f9a80', '#3f86ad'];
  #width = 0;
  #height = 0;
  #dpr = 1;
  /** The prism's point (canvas CSS pixels): where the stream splits. */
  #ox = 0;
  #oy = 0;
  #emit = 0;
  #time = 0;
  #last = 0;
  #frame = 0;
  #running = false;
  #visible = true;
  #colorsRead = -Infinity;
  readonly #still: boolean;

  constructor(canvas: HTMLCanvasElement, options: Partial<PrismFieldOptions> = {}, init: { still?: boolean } = {}) {
    this.#canvas = canvas;
    this.#ctx = canvas.getContext('2d');
    this.#options = { ...PrismFieldDefaults, ...options };
    this.#random = generator(this.#options.seed);
    this.#still = !!init.still;
  }

  setOptions(options: Partial<PrismFieldOptions>): void {
    const seed = this.#options.seed;
    this.#options = { ...this.#options, ...options };
    if (this.#options.seed !== seed) {
      this.#random = generator(this.#options.seed);
      this.#alive = 0;
    }
    if (this.#still) this.#drawStill();
  }

  /** Where the stream splits, in the canvas's CSS pixels (the mark's prism). */
  setOrigin(x: number, y: number): void {
    this.#ox = x;
    this.#oy = y;
    if (this.#still) this.#drawStill();
  }

  /** Reads the canvas's size (and the screen's pixel ratio) and the theme's colors. */
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

  /** The theme's colors (it may change under the canvas: the chrome follows the canvas background). */
  readColors(): void {
    if (typeof getComputedStyle !== 'function') return;
    // Through the element's `color`, so a `color-mix()` token comes back as a color a canvas takes.
    const el = this.#canvas;
    const before = el.style.color;
    for (let i = 0; i < Colors; i++) {
      el.style.color = `var(${ColorVars[i]})`;
      const v = getComputedStyle(el).color;
      if (v) this.#colors[i] = v;
    }
    el.style.color = before;
  }

  /** Runs while visible; nothing is drawn or scheduled while stopped or hidden. */
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
    // A long pause (a hidden tab) resumes where it was rather than catching up.
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
    this.#alive = 0;
    this.#random = generator(this.#options.seed);
    for (let t = 0; t < StillSeconds; t += 1 / 30) this.#step(1 / 30);
    this.#draw();
  }

  #spawn(): void {
    if (this.#alive >= Capacity) return;
    const o = this.#options;
    const i = this.#alive++;
    this.#x[i] = this.#ox - o.beam;
    this.#y[i] = this.#oy + (this.#random() * 2 - 1) * o.width;
    // Toward the prism's point: the stream narrows into it.
    const dy = this.#oy - this.#y[i];
    const len = Math.hypot(o.beam, dy) || 1;
    const v = o.speed * (0.85 + this.#random() * 0.3);
    this.#vx[i] = (v * o.beam) / len;
    this.#vy[i] = (v * dy) / len;
    this.#age[i] = 0;
    this.#color[i] = 0;
    this.#shimmer[i] = this.#random() * Math.PI * 2;
  }

  #kill(i: number): void {
    const last = --this.#alive;
    this.#x[i] = this.#x[last];
    this.#y[i] = this.#y[last];
    this.#vx[i] = this.#vx[last];
    this.#vy[i] = this.#vy[last];
    this.#age[i] = this.#age[last];
    this.#color[i] = this.#color[last];
    this.#shimmer[i] = this.#shimmer[last];
  }

  #step(dt: number): void {
    const o = this.#options;
    this.#time += dt;
    this.#emit += dt * o.rate;
    while (this.#emit >= 1) {
      this.#emit -= 1;
      this.#spawn();
    }
    const spread = (o.spread * Math.PI) / 180;
    for (let i = 0; i < this.#alive; i++) {
      this.#age[i] += dt;
      this.#x[i] += this.#vx[i] * dt;
      this.#y[i] += this.#vy[i] * dt;
      if (this.#color[i] === 0) {
        // Through the prism: one of the logo's colors, bent its own way (the first up, the last down), a little
        // scattered.
        if (this.#x[i] >= this.#ox) {
          const c = 1 + Math.floor(this.#random() * 4);
          this.#color[i] = c;
          const angle = spread * ((c - 1) / 1.5 - 1) + (this.#random() * 2 - 1) * spread * 0.22;
          const v = Math.hypot(this.#vx[i], this.#vy[i]);
          this.#vx[i] = v * Math.cos(angle);
          this.#vy[i] = v * Math.sin(angle);
          this.#x[i] = this.#ox;
          this.#y[i] = this.#oy;
        }
      } else {
        const dx = this.#x[i] - this.#ox;
        const dy = this.#y[i] - this.#oy;
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
      // The stream: in over its first third, at a third of the most.
      const along = 1 - (this.#ox - this.#x[i]) / o.beam;
      a = (o.alpha / 3) * Math.min(1, Math.max(0, along * 3));
    } else {
      const d = Math.hypot(this.#x[i] - this.#ox, this.#y[i] - this.#oy) / o.reach;
      // Brightest just past the prism, gone at the reach.
      a = o.alpha * Math.min(1, d * 12) * (1 - d) ** 1.6;
    }
    a *= 0.7 + 0.3 * Math.sin(this.#time * 1.7 + this.#shimmer[i]);
    return Math.max(0, Math.min(Levels - 1, Math.round((a / Math.max(o.alpha, 1e-6)) * (Levels - 1))));
  }

  #draw(): void {
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
    for (let i = 0; i < this.#alive; i++) {
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
    for (let b = 0; b < counts.length; b++) counts[b] = 0;
    for (let i = 0; i < this.#alive; i++) {
      const b = this.#bucketOf[i];
      this.#bucket[this.#starts[b] + counts[b]++] = i;
    }

    ctx.lineCap = 'round';
    ctx.lineWidth = o.size * 2;
    for (let c = 0; c < Colors; c++) {
      ctx.strokeStyle = this.#colors[c];
      ctx.fillStyle = this.#colors[c];
      for (let l = 1; l < Levels; l++) {
        const b = c * Levels + l;
        const n = counts[b];
        if (!n) continue;
        ctx.globalAlpha = (o.alpha * l) / (Levels - 1);
        ctx.beginPath();
        const s = this.#starts[b];
        for (let k = 0; k < n; k++) {
          const i = this.#bucket[s + k];
          const x = this.#x[i];
          const y = this.#y[i];
          const tail = this.#color[i] === 0 ? o.tail * 0.6 : o.tail;
          if (tail > 0) {
            const v = Math.hypot(this.#vx[i], this.#vy[i]) || 1;
            ctx.moveTo(x - (this.#vx[i] / v) * tail, y - (this.#vy[i] / v) * tail);
            ctx.lineTo(x, y);
          } else {
            ctx.moveTo(x + 0.01, y);
            ctx.lineTo(x, y);
          }
        }
        ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
  }
}
