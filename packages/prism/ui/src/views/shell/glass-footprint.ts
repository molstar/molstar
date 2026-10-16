/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Where the wordmark's glass is, on a grid (the splashes', prism-field.ts and quantum-field.ts): the word filled (its
// letters and star as solid glass, a little rounded), sampled once, then carried piece by piece where each piece of the
// `Glass` is now. The quantum splash takes its potential from it; both keep their light out of sight under it, so the
// word stays clear, whole or in pieces. Allocation-free once made (memory: render paths).

import { WordmarkGeometry } from '../../components/icon/wordmark.js';
import type { Glass } from './glass.js';

/** Samples a drawing unit, and how much the filled word is softened (units). */
const MaskRes = 4;
const MaskSoft = 0.8;

export class GlassFootprint {
  /** The filled word (0..1, `MaskRes` samples a unit from (−1, −1)); none where there is no canvas (Node, jsdom). */
  readonly #mask: Float32Array | undefined;
  readonly #width: number = 0;
  readonly #height: number = 0;
  readonly #pose = new Float32Array(5);

  constructor() {
    if (typeof document === 'undefined' || typeof Path2D !== 'function') return;
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
    this.#width = w;
    this.#height = h;
  }

  /** Whether there is a footprint to give (a canvas made the mask). */
  get ready(): boolean {
    return !!this.#mask;
  }

  /**
   * How much glass covers each cell of a grid, into `out` (cleared first): `nx` × `ny` cells of `cell` px from (gx, gy),
   * canvas CSS px; the word's drawing origin at (wordX, wordY), `scale` px a unit. A cell takes, over 2 × 2 samples,
   * the word where the piece over it came from (the most of any piece over it).
   */
  rasterize(
    glass: Glass,
    wordX: number,
    wordY: number,
    scale: number,
    gx: number,
    gy: number,
    cell: number,
    nx: number,
    ny: number,
    out: Float32Array,
  ): void {
    out.fill(0, 0, nx * ny);
    const mask = this.#mask;
    if (!mask) return;
    const pose = this.#pose;
    const corners = glass.corners;
    const mw = this.#width;
    const mh = this.#height;
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
      let x0 = Number.POSITIVE_INFINITY;
      let x1 = Number.NEGATIVE_INFINITY;
      let y0 = Number.POSITIVE_INFINITY;
      let y1 = Number.NEGATIVE_INFINITY;
      for (let k = 0; k < 3; k++) {
        const lx = corners[i * 6 + 2 * k];
        const ly = corners[i * 6 + 2 * k + 1];
        const px = (wordX + (mx + lx * c - ly * sn) * scale - gx) / cell;
        const py = (wordY + (my + lx * sn + ly * c) * scale - gy) / cell;
        if (px < x0) x0 = px;
        if (px > x1) x1 = px;
        if (py < y0) y0 = py;
        if (py > y1) y1 = py;
      }
      const cx0 = Math.max(0, Math.floor(x0));
      const cx1 = Math.min(nx - 1, Math.ceil(x1));
      const cy0 = Math.max(0, Math.floor(y0));
      const cy1 = Math.min(ny - 1, Math.ceil(y1));
      for (let y = cy0; y <= cy1; y++) {
        for (let x = cx0; x <= cx1; x++) {
          let sum = 0;
          for (let sy = 0; sy < 2; sy++) {
            for (let sx = 0; sx < 2; sx++) {
              // The sample in the piece's own frame (drawing units about its middle).
              const du = (gx + (x + 0.25 + 0.5 * sx) * cell - wordX) / scale - mx;
              const dv = (gy + (y + 0.25 + 0.5 * sy) * cell - wordY) / scale - my;
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
          const k = y * nx + x;
          if (sum / 4 > out[k]) out[k] = sum / 4;
        }
      }
    }
  }
}
