/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ribbons' geometry (_spec/README.md, Stems, Ribbons, Stubs): one frame for a ring and a row. A point is given by
// its place along the band (`t`: an angle on a ring, a page x on a row) and its depth inward from the band's inner edge
// (`d`, CSS px: toward the centre on a ring, down on a row). Cubics are flattened at a fixed number of segments, and the
// polygons drawn are the polygons hit (one cover: the closed hit equals the drawn run).

/** Where the ribbons are drawn: a ring's hole, or a row's band. */
export class RibbonFrame {
  polar = true;
  /** A ring: its centre and the radius of the band's inner edge (the mouths); a row: the band's top and its height. */
  cx = 0;
  cy = 0;
  r = 0;
  /** The room inward: the hole's radius on a ring, the band's height on a row. */
  room = 0;
  /** A row's band: its left and right edge (the room a ribbon arcs over). */
  x0 = 0;
  x1 = 0;
  /** `point`'s result. */
  px = 0;
  py = 0;

  /** The page point at place `t` and depth `d`, into `px`, `py`. */
  point(t: number, d: number): void {
    if (this.polar) {
      const r = Math.max(0, this.r - d);
      this.px = this.cx + r * Math.sin(t);
      this.py = this.cy - r * Math.cos(t);
    } else {
      this.px = t;
      this.py = this.r + d;
    }
  }

  /** Place units a CSS px at depth `d` (radians a px on a ring, 1 on a row). */
  perPx(d: number): number {
    return this.polar ? 1 / Math.max(1e-6, this.r - d) : 1;
  }
}

/** A growable list of polygons: their points (x, y) end to end, each one's start and length. */
export class Polygons {
  n = 0;
  start = new Int32Array(64);
  len = new Int32Array(64);
  /** Bounding boxes: x0, y0, x1, y1 a polygon. */
  box = new Float64Array(256);
  xy = new Float64Array(4096);
  m = 0;

  reset(): void {
    this.n = 0;
    this.m = 0;
  }

  begin(): number {
    const i = this.n;
    if (i >= this.start.length) {
      this.start = grow(this.start, i + 1);
      this.len = grow(this.len, i + 1);
      this.box = grow(this.box, 4 * (i + 1));
    }
    this.start[i] = this.m;
    this.len[i] = 0;
    this.box[4 * i] = Infinity;
    this.box[4 * i + 1] = Infinity;
    this.box[4 * i + 2] = -Infinity;
    this.box[4 * i + 3] = -Infinity;
    return i;
  }

  add(x: number, y: number): void {
    const i = this.n;
    if (2 * this.m + 2 > this.xy.length) this.xy = grow(this.xy, 2 * this.m + 2);
    this.xy[2 * this.m] = x;
    this.xy[2 * this.m + 1] = y;
    this.m++;
    this.len[i]++;
    const b = this.box;
    if (x < b[4 * i]) b[4 * i] = x;
    if (y < b[4 * i + 1]) b[4 * i + 1] = y;
    if (x > b[4 * i + 2]) b[4 * i + 2] = x;
    if (y > b[4 * i + 3]) b[4 * i + 3] = y;
  }

  end(): number {
    return this.n++;
  }

  /** Whether page point (x, y) is inside polygon `i` (even-odd). */
  contains(i: number, x: number, y: number): boolean {
    const b = this.box;
    if (x < b[4 * i] || x > b[4 * i + 2] || y < b[4 * i + 1] || y > b[4 * i + 3]) return false;
    const s = this.start[i];
    const n = this.len[i];
    const p = this.xy;
    let inside = false;
    for (let k = 0, j = n - 1; k < n; j = k++) {
      const xk = p[2 * (s + k)];
      const yk = p[2 * (s + k) + 1];
      const xj = p[2 * (s + j)];
      const yj = p[2 * (s + j) + 1];
      if (yk > y !== yj > y && x < ((xj - xk) * (y - yk)) / (yj - yk) + xk) inside = !inside;
    }
    return inside;
  }
}

export function grow<T extends Int32Array | Float64Array | Uint32Array | Uint8Array | Int8Array>(
  a: T,
  need: number,
): T {
  if (need <= a.length) return a;
  const out = new (a.constructor as { new (n: number): T })(Math.max(need, 2 * a.length));
  out.set(a);
  return out;
}

/**
 * Adds a cubic from (t0, d0) to (t3, d3) with controls (t1, d1) and (t2, d2), in place units and depths, flattened
 * into `segments` lines (its first point left out: the polygon has it).
 */
export function addCubic(
  f: RibbonFrame,
  out: Polygons,
  segments: number,
  t0: number,
  d0: number,
  t1: number,
  d1: number,
  t2: number,
  d2: number,
  t3: number,
  d3: number,
): void {
  // the controls on the page: a cubic in the page's frame (a ring's chord bends through its centre)
  f.point(t0, d0);
  const x0 = f.px;
  const y0 = f.py;
  f.point(t1, d1);
  const x1 = f.px;
  const y1 = f.py;
  f.point(t2, d2);
  const x2 = f.px;
  const y2 = f.py;
  f.point(t3, d3);
  const x3 = f.px;
  const y3 = f.py;
  for (let k = 1; k <= segments; k++) {
    const u = k / segments;
    const v = 1 - u;
    const a = v * v * v;
    const b = 3 * v * v * u;
    const c = 3 * v * u * u;
    const e = u * u * u;
    out.add(a * x0 + b * x1 + c * x2 + e * x3, a * y0 + b * y1 + c * y2 + e * y3);
  }
}

/** Adds an arc at depth `d` from place `ta` to `tb` (lines a few px apart; its first point left out). */
export function addArc(f: RibbonFrame, out: Polygons, ta: number, tb: number, d: number): void {
  const px = Math.abs(tb - ta) / f.perPx(d);
  const n = f.polar ? Math.max(1, Math.min(64, Math.ceil(px / 4))) : 1;
  for (let k = 1; k <= n; k++) {
    f.point(ta + ((tb - ta) * k) / n, d);
    out.add(f.px, f.py);
  }
}

/**
 * Whether a ribbon is crossed (_spec/README.md, Ribbons, The twist): `s = strand × the display direction of A × that
 * of B`; −1 crossed (low edge to low edge), +1 or a strand of 0 uncrossed (A's low edge to B's high edge).
 */
export function crossed(strand: number, dirA: number, dirB: number): boolean {
  return strand * dirA * dirB < 0;
}

/**
 * A ribbon's four corners (_spec/README.md, Ribbons, The twist), into `out`: its first side from `out[0]` (at A) to
 * `out[1]` (at B), its second from `out[2]` (at B) back to `out[3]` (at A). A slice's low edge is its first in display
 * order (clockwise on a ring); uncrossed (`s` +1), A's low edge joins B's high (two arcs facing each other across the
 * hole have opposite tangents: the flat ribbon), crossed (`s` −1) low joins low. The blocks' directions are in `s`.
 */
export function ribbonCorners(
  a0: number,
  a1: number,
  dirA: number,
  b0: number,
  b1: number,
  dirB: number,
  strand: number,
  out: Float64Array,
): void {
  // the edges in display order (a slice is a share of the count, not positions); the blocks' directions are in `s`
  const aLow = a0;
  const aHigh = a1;
  const bLow = b0;
  const bHigh = b1;
  const x = crossed(strand, dirA, dirB);
  out[0] = aLow;
  out[1] = x ? bLow : bHigh;
  out[2] = x ? bHigh : bLow;
  out[3] = aHigh;
}

/** Twice the signed area of a triangle: the side of line (a, b) point c is on. */
function orient(ax: number, ay: number, bx: number, by: number, cx: number, cy: number): number {
  return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

/** Whether segments (p1, p2) and (p3, p4) cross (properly). */
export function segmentsCross(
  x1: number,
  y1: number,
  x2: number,
  y2: number,
  x3: number,
  y3: number,
  x4: number,
  y4: number,
): boolean {
  const a = orient(x1, y1, x2, y2, x3, y3);
  const b = orient(x1, y1, x2, y2, x4, y4);
  const c = orient(x3, y3, x4, y4, x1, y1);
  const d = orient(x3, y3, x4, y4, x2, y2);
  return a * b < 0 && c * d < 0;
}
