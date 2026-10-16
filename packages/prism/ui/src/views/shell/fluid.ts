/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The air the splash's light goes through (prism-field.ts): incompressible Navier–Stokes on a coarse grid, after Stam's
// stable fluids. Each step: advection (semi-Lagrangian, its backtrace a midpoint Runge–Kutta step), viscosity as a drag,
// vorticity confinement (the eddies keep their spin), and a pressure projection (Gauss–Seidel, warm started) that keeps
// it divergence-free, so a push goes round as a swirl rather than piling up. Velocities sit at the cells' middles (px/s,
// in the canvas's CSS pixels); outside the grid the edge's. Arrays are allocated only when the grid grows (memory: render
// paths stay allocation-free).

/** A cell's side (px), and the cells at most (the grid coarsens past them). */
const Cell = 12;
const MaxCells = 4800;
/** Pressure iterations a step. */
const Iterations = 10;

export class Fluid {
  /** The velocity where `sample` last looked (px/s). */
  su = 0;
  sv = 0;
  /** Viscosity, as a drag (1/s): how fast motion dies down. */
  drag = 0.6;
  /** Vorticity confinement: how much the eddies keep their spin. */
  confinement = 0;

  #x0 = 0;
  #y0 = 0;
  #nx = 0;
  #ny = 0;
  #h = Cell;
  #u = new Float32Array(0);
  #v = new Float32Array(0);
  #u0 = new Float32Array(0);
  #v0 = new Float32Array(0);
  #p = new Float32Array(0);
  #w = new Float32Array(0);
  // A lookup's cell and weights (`#locate`).
  #k = 0;
  #fx = 0;
  #fy = 0;

  get on(): boolean {
    return this.#nx >= 2 && this.#ny >= 2;
  }

  /** The box it covers (canvas CSS px); the field is kept if the grid's size stays, else it starts still. */
  layout(x: number, y: number, width: number, height: number): void {
    this.#x0 = x;
    this.#y0 = y;
    let h = Cell;
    if ((width / h) * (height / h) > MaxCells) h = Math.sqrt((width * height) / MaxCells);
    const nx = Math.max(0, Math.ceil(width / h));
    const ny = Math.max(0, Math.ceil(height / h));
    if (nx === this.#nx && ny === this.#ny && h === this.#h) return;
    this.#nx = nx;
    this.#ny = ny;
    this.#h = h;
    const n = nx * ny;
    if (this.#u.length < n) {
      this.#u = new Float32Array(n);
      this.#v = new Float32Array(n);
      this.#u0 = new Float32Array(n);
      this.#v0 = new Float32Array(n);
      this.#p = new Float32Array(n);
      this.#w = new Float32Array(n);
    }
    this.clear();
  }

  /** Still again. */
  clear(): void {
    this.#u.fill(0);
    this.#v.fill(0);
    this.#p.fill(0);
  }

  /** The velocity at (x, y) (canvas px) into `su`, `sv`. */
  sample(x: number, y: number): void {
    if (!this.on) {
      this.su = 0;
      this.sv = 0;
      return;
    }
    this.#locate((x - this.#x0) / this.#h - 0.5, (y - this.#y0) / this.#h - 0.5);
    this.su = this.#lerp(this.#u);
    this.sv = this.#lerp(this.#v);
  }

  /** Draws the air within `radius` of (x, y) toward the velocity (tx, ty) (px/s), by up to `amount` (0..1) at the middle. */
  drive(x: number, y: number, tx: number, ty: number, radius: number, amount: number): void {
    if (!this.on) return;
    const nx = this.#nx;
    const h = this.#h;
    const r2 = radius * radius;
    const i0 = Math.max(0, Math.floor((x - radius - this.#x0) / h));
    const i1 = Math.min(nx - 1, Math.ceil((x + radius - this.#x0) / h));
    const j0 = Math.max(0, Math.floor((y - radius - this.#y0) / h));
    const j1 = Math.min(this.#ny - 1, Math.ceil((y + radius - this.#y0) / h));
    for (let j = j0; j <= j1; j++) {
      const dy = this.#y0 + (j + 0.5) * h - y;
      for (let i = i0; i <= i1; i++) {
        const dx = this.#x0 + (i + 0.5) * h - x;
        const d2 = dx * dx + dy * dy;
        if (d2 >= r2) continue;
        const f = 1 - Math.sqrt(d2) / radius;
        const a = amount * f * f;
        const k = j * nx + i;
        this.#u[k] += (tx - this.#u[k]) * a;
        this.#v[k] += (ty - this.#v[k]) * a;
      }
    }
  }

  /**
   * Spins the air round (x, y): adds up to `speed` (px/s; positive: clockwise on the screen) across it, the most a third
   * of the way out, none at the middle or past `radius`.
   */
  swirl(x: number, y: number, speed: number, radius: number): void {
    if (!this.on) return;
    const nx = this.#nx;
    const h = this.#h;
    const r2 = radius * radius;
    const i0 = Math.max(0, Math.floor((x - radius - this.#x0) / h));
    const i1 = Math.min(nx - 1, Math.ceil((x + radius - this.#x0) / h));
    const j0 = Math.max(0, Math.floor((y - radius - this.#y0) / h));
    const j1 = Math.min(this.#ny - 1, Math.ceil((y + radius - this.#y0) / h));
    for (let j = j0; j <= j1; j++) {
      const dy = this.#y0 + (j + 0.5) * h - y;
      for (let i = i0; i <= i1; i++) {
        const dx = this.#x0 + (i + 0.5) * h - x;
        const d2 = dx * dx + dy * dy;
        if (d2 >= r2 || d2 < 1e-6) continue;
        const d = Math.sqrt(d2);
        const q = d / radius;
        // Peaks at a third of the way out (6.75 q (1 - q)² is 1 there).
        const s = (speed * 6.75 * q * (1 - q) * (1 - q)) / d;
        const k = j * nx + i;
        this.#u[k] -= dy * s;
        this.#v[k] += dx * s;
      }
    }
  }

  /** The air moves on by `dt` (s). */
  step(dt: number): void {
    if (!this.on) return;
    this.#advect(dt);
    if (this.confinement > 0) this.#confine(dt);
    this.#project();
  }

  /** Each cell takes the velocity from where its air came from (traced back along the flow by a midpoint step). */
  #advect(dt: number): void {
    const nx = this.#nx;
    const ny = this.#ny;
    const n = nx * ny;
    const u = this.#u;
    const v = this.#v;
    const u0 = this.#u0;
    const v0 = this.#v0;
    u0.set(u.subarray(0, n));
    v0.set(v.subarray(0, n));
    const s = dt / this.#h;
    const fade = Math.exp(-this.drag * dt);
    for (let j = 0, k = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++, k++) {
        this.#locate(i - 0.5 * s * u0[k], j - 0.5 * s * v0[k]);
        const um = this.#lerp(u0);
        const vm = this.#lerp(v0);
        this.#locate(i - s * um, j - s * vm);
        u[k] = this.#lerp(u0) * fade;
        v[k] = this.#lerp(v0) * fade;
      }
    }
  }

  /** Vorticity confinement: a push along the eddies' edges, toward their cores, so they keep turning. */
  #confine(dt: number): void {
    const nx = this.#nx;
    const ny = this.#ny;
    const u = this.#u;
    const v = this.#v;
    const w = this.#w;
    const mag = this.#u0;
    const h2 = 2 * this.#h;
    for (let j = 0, k = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++, k++) {
        const r = i < nx - 1 ? k + 1 : k;
        const l = i > 0 ? k - 1 : k;
        const d = j < ny - 1 ? k + nx : k;
        const t = j > 0 ? k - nx : k;
        w[k] = (v[r] - v[l] - (u[d] - u[t])) / h2;
        mag[k] = Math.abs(w[k]);
      }
    }
    const eps = this.confinement * this.#h * dt;
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1, k = j * nx + 1; i < nx - 1; i++, k++) {
        const gx = mag[k + 1] - mag[k - 1];
        const gy = mag[k + nx] - mag[k - nx];
        const len = Math.sqrt(gx * gx + gy * gy) + 1e-5;
        u[k] += eps * (gy / len) * w[k];
        v[k] -= eps * (gx / len) * w[k];
      }
    }
  }

  /** Takes out the divergence: solves for the pressure (open edges: none past them) and subtracts its gradient. */
  #project(): void {
    const nx = this.#nx;
    const ny = this.#ny;
    const u = this.#u;
    const v = this.#v;
    const p = this.#p;
    const div = this.#w;
    const h = this.#h;
    const h2 = 2 * h;
    for (let j = 0, k = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++, k++) {
        const r = i < nx - 1 ? u[k + 1] : u[k];
        const l = i > 0 ? u[k - 1] : u[k];
        const d = j < ny - 1 ? v[k + nx] : v[k];
        const t = j > 0 ? v[k - nx] : v[k];
        div[k] = ((r - l + d - t) / h2) * h * h;
      }
    }
    for (let it = 0; it < Iterations; it++) {
      for (let j = 0, k = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++, k++) {
          const r = i < nx - 1 ? p[k + 1] : 0;
          const l = i > 0 ? p[k - 1] : 0;
          const d = j < ny - 1 ? p[k + nx] : 0;
          const t = j > 0 ? p[k - nx] : 0;
          p[k] = (r + l + d + t - div[k]) * 0.25;
        }
      }
    }
    for (let j = 0, k = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++, k++) {
        const r = i < nx - 1 ? p[k + 1] : 0;
        const l = i > 0 ? p[k - 1] : 0;
        const d = j < ny - 1 ? p[k + nx] : 0;
        const t = j > 0 ? p[k - nx] : 0;
        u[k] -= (r - l) / h2;
        v[k] -= (d - t) / h2;
      }
    }
  }

  /** A lookup at grid coordinates (cells, from the first cell's middle; clamped to the grid): its cell and weights. */
  #locate(gx: number, gy: number): void {
    const nx = this.#nx;
    const ny = this.#ny;
    if (gx < 0) gx = 0;
    else if (gx > nx - 1) gx = nx - 1;
    if (gy < 0) gy = 0;
    else if (gy > ny - 1) gy = ny - 1;
    let i = gx | 0;
    let j = gy | 0;
    if (i > nx - 2) i = nx - 2;
    if (j > ny - 2) j = ny - 2;
    this.#fx = gx - i;
    this.#fy = gy - j;
    this.#k = j * nx + i;
  }

  #lerp(f: Float32Array): number {
    const k = this.#k;
    const nx = this.#nx;
    const fx = this.#fx;
    const fy = this.#fy;
    return (f[k] * (1 - fx) + f[k + 1] * fx) * (1 - fy) + (f[k + nx] * (1 - fx) + f[k + nx + 1] * fx) * fy;
  }
}
