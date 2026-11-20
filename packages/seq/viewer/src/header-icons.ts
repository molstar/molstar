/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type HeaderIcon, PathCommand } from '@molstar/seqx-core/render';
import { lucideEyeOff } from '@molstar/ui/components/icon/_generated/lucide/eye-off';
import { lucidePalette } from '@molstar/ui/components/icon/_generated/lucide/palette';

// A header entry's icon as path data (_spec/README.md "Behaviors", Headers: icons as paths): its markup (Lucide's, in a
// 24-unit box: the entry's own `iconMarkup`, else a built-in's icon, imported by itself from `@molstar/ui`) read once
// into move, line, cubic and close commands, what core's draw lists stroke. The elements read: `path`, `line`,
// `circle`, `polyline`, `polygon`, `rect`, `ellipse`; another makes the icon unreadable (undefined: the slot draws its
// label's first letter). No map of every icon: a host's entry brings its own markup.

/** The built-in entries' icons by name. */
const BuiltIn: Readonly<Record<string, string>> = { 'eye-off': lucideEyeOff, palette: lucidePalette };

const Cache = new Map<string, HeaderIcon | null>();

/** An icon's paths (its markup, else a built-in's by name), read once; undefined when no path reads it. */
export function headerIcon(name: string | undefined, markup?: string): HeaderIcon | undefined {
  const text = markup ?? (name === undefined ? undefined : BuiltIn[name]);
  if (text === undefined) return undefined;
  let out = Cache.get(text);
  if (out === undefined) {
    out = readMarkup(text) ?? null;
    Cache.set(text, out);
  }
  return out ?? undefined;
}

/** Path data being built. */
class Out {
  readonly commands: number[] = [];
  readonly points: number[] = [];
  move(x: number, y: number) {
    this.commands.push(PathCommand.move);
    this.points.push(x, y);
  }
  line(x: number, y: number) {
    this.commands.push(PathCommand.line);
    this.points.push(x, y);
  }
  cubic(a: number, b: number, c: number, d: number, x: number, y: number) {
    this.commands.push(PathCommand.cubic);
    this.points.push(a, b, c, d, x, y);
  }
  close() {
    this.commands.push(PathCommand.close);
  }
}

/** The markup's elements read into one path; undefined at an element it does not read. */
export function readMarkup(markup: string): HeaderIcon | undefined {
  const out = new Out();
  for (const m of markup.matchAll(/<(\w+)([^>]*)\/?>/g)) {
    const [, tag, attrs] = m;
    const a: Record<string, string> = {};
    for (const p of attrs.matchAll(/([\w-]+)="([^"]*)"/g)) a[p[1]] = p[2];
    const n = (k: string, d = 0) => (a[k] === undefined ? d : Number.parseFloat(a[k]));
    switch (tag) {
      case 'path':
        if (!readPath(a.d ?? '', out)) return undefined;
        break;
      case 'line':
        out.move(n('x1'), n('y1'));
        out.line(n('x2'), n('y2'));
        break;
      case 'polyline':
      case 'polygon': {
        const v = numbers(a.points ?? '');
        for (let i = 0; i + 1 < v.length; i += 2) (i ? out.line : out.move).call(out, v[i], v[i + 1]);
        if (tag === 'polygon') out.close();
        break;
      }
      case 'circle':
        ellipse(out, n('cx'), n('cy'), n('r'), n('r'));
        break;
      case 'ellipse':
        ellipse(out, n('cx'), n('cy'), n('rx'), n('ry'));
        break;
      case 'rect': {
        const x = n('x');
        const y = n('y');
        const w = n('width');
        const h = n('height');
        const r = Math.min(n('rx', n('ry')), w / 2, h / 2);
        if (!(r > 0)) {
          out.move(x, y);
          out.line(x + w, y);
          out.line(x + w, y + h);
          out.line(x, y + h);
          out.close();
          break;
        }
        const k = r * Kappa;
        out.move(x + r, y);
        out.line(x + w - r, y);
        out.cubic(x + w - r + k, y, x + w, y + r - k, x + w, y + r);
        out.line(x + w, y + h - r);
        out.cubic(x + w, y + h - r + k, x + w - r + k, y + h, x + w - r, y + h);
        out.line(x + r, y + h);
        out.cubic(x + r - k, y + h, x, y + h - r + k, x, y + h - r);
        out.line(x, y + r);
        out.cubic(x, y + r - k, x + r - k, y, x + r, y);
        out.close();
        break;
      }
      default:
        return undefined;
    }
  }
  return { commands: Int8Array.from(out.commands), points: Float32Array.from(out.points) };
}

/** A quarter circle's control distance, over its radius. */
const Kappa = 0.5522847498;

function ellipse(out: Out, cx: number, cy: number, rx: number, ry: number) {
  const kx = rx * Kappa;
  const ky = ry * Kappa;
  out.move(cx + rx, cy);
  out.cubic(cx + rx, cy + ky, cx + kx, cy + ry, cx, cy + ry);
  out.cubic(cx - kx, cy + ry, cx - rx, cy + ky, cx - rx, cy);
  out.cubic(cx - rx, cy - ky, cx - kx, cy - ry, cx, cy - ry);
  out.cubic(cx + kx, cy - ry, cx + rx, cy - ky, cx + rx, cy);
  out.close();
}

/** A list of numbers as SVG writes them (`1.5-2`, `.5.5`, exponents). */
function numbers(s: string): number[] {
  return [...s.matchAll(/[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g)].map((m) => Number(m[0]));
}

/** Path data's commands read into `out`: every SVG command, quadratics and arcs as cubics; false when it is malformed. */
function readPath(d: string, out: Out): boolean {
  let x = 0;
  let y = 0;
  let sx = 0;
  let sy = 0;
  // the last control point (a smooth curve's reflection)
  let cx = 0;
  let cy = 0;
  let last = '';
  for (const m of d.matchAll(/([MmLlHhVvCcSsQqTtAaZz])([^MmLlHhVvCcSsQqTtAaZz]*)/g)) {
    const c = m[1];
    const rel = c === c.toLowerCase();
    const C = c.toUpperCase();
    // an arc's flags may be written run together (`0 011 1`): read one character each
    const v = C === 'A' ? arcNumbers(m[2]) : numbers(m[2]);
    if (!v) return false;
    const ox = rel ? x : 0;
    const oy = rel ? y : 0;
    const per = { M: 2, L: 2, H: 1, V: 1, C: 6, S: 4, Q: 4, T: 2, A: 7, Z: 0 }[C] as number;
    if (C === 'Z') {
      out.close();
      x = sx;
      y = sy;
      last = C;
      continue;
    }
    if (!v.length || v.length % per) return false;
    for (let i = 0; i < v.length; i += per) {
      const ax = rel ? x : 0;
      const ay = rel ? y : 0;
      switch (C) {
        case 'M':
          if (i === 0) {
            x = v[i] + ox;
            y = v[i + 1] + oy;
            out.move(x, y);
            sx = x;
            sy = y;
          } else {
            x = v[i] + ax;
            y = v[i + 1] + ay;
            out.line(x, y);
          }
          break;
        case 'L':
          x = v[i] + ax;
          y = v[i + 1] + ay;
          out.line(x, y);
          break;
        case 'H':
          x = v[i] + ax;
          out.line(x, y);
          break;
        case 'V':
          y = v[i] + ay;
          out.line(x, y);
          break;
        case 'C':
          out.cubic(v[i] + ax, v[i + 1] + ay, (cx = v[i + 2] + ax), (cy = v[i + 3] + ay), v[i + 4] + ax, v[i + 5] + ay);
          x = v[i + 4] + ax;
          y = v[i + 5] + ay;
          break;
        case 'S': {
          const smooth = last === 'C' || last === 'S';
          const x1 = smooth ? 2 * x - cx : x;
          const y1 = smooth ? 2 * y - cy : y;
          out.cubic(x1, y1, (cx = v[i] + ax), (cy = v[i + 1] + ay), v[i + 2] + ax, v[i + 3] + ay);
          x = v[i + 2] + ax;
          y = v[i + 3] + ay;
          break;
        }
        case 'Q':
        case 'T': {
          const smooth = C === 'T' && (last === 'Q' || last === 'T');
          const qx = C === 'Q' ? v[i] + ax : smooth ? 2 * x - cx : x;
          const qy = C === 'Q' ? v[i + 1] + ay : smooth ? 2 * y - cy : y;
          const ex = v[i + (C === 'Q' ? 2 : 0)] + ax;
          const ey = v[i + (C === 'Q' ? 3 : 1)] + ay;
          out.cubic(
            x + (2 / 3) * (qx - x),
            y + (2 / 3) * (qy - y),
            ex + (2 / 3) * (qx - ex),
            ey + (2 / 3) * (qy - ey),
            ex,
            ey,
          );
          cx = qx;
          cy = qy;
          x = ex;
          y = ey;
          break;
        }
        case 'A': {
          const ex = v[i + 5] + ax;
          const ey = v[i + 6] + ay;
          arc(out, x, y, v[i], v[i + 1], v[i + 2], v[i + 3] !== 0, v[i + 4] !== 0, ex, ey);
          x = ex;
          y = ey;
          break;
        }
      }
      last = C === 'M' ? 'L' : C;
    }
    last = C;
  }
  return true;
}

/** An arc's numbers: radii, rotation, then its two flags one character each, then its end. */
function arcNumbers(s: string): number[] | undefined {
  const out: number[] = [];
  let rest = s;
  const num = /^[\s,]*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)/;
  const flag = /^[\s,]*([01])/;
  while (rest.trim()) {
    for (let k = 0; k < 7; k++) {
      const m = (k === 3 || k === 4 ? flag : num).exec(rest);
      if (!m) return undefined;
      out.push(Number(m[1]));
      rest = rest.slice(m[0].length);
    }
  }
  return out;
}

/** An SVG arc as cubics (the endpoint parametrization to the centre's, then a cubic per quarter at most). */
function arc(
  out: Out,
  x0: number,
  y0: number,
  rx: number,
  ry: number,
  angle: number,
  large: boolean,
  sweep: boolean,
  x: number,
  y: number,
) {
  rx = Math.abs(rx);
  ry = Math.abs(ry);
  if (!rx || !ry || (x0 === x && y0 === y)) {
    out.line(x, y);
    return;
  }
  const phi = (angle * Math.PI) / 180;
  const cos = Math.cos(phi);
  const sin = Math.sin(phi);
  const dx = (x0 - x) / 2;
  const dy = (y0 - y) / 2;
  const x1 = cos * dx + sin * dy;
  const y1 = -sin * dx + cos * dy;
  const l = (x1 * x1) / (rx * rx) + (y1 * y1) / (ry * ry);
  if (l > 1) {
    rx *= Math.sqrt(l);
    ry *= Math.sqrt(l);
  }
  const num = rx * rx * ry * ry - rx * rx * y1 * y1 - ry * ry * x1 * x1;
  const den = rx * rx * y1 * y1 + ry * ry * x1 * x1;
  const k = (large === sweep ? -1 : 1) * Math.sqrt(Math.max(0, num / den));
  const cxp = (k * rx * y1) / ry;
  const cyp = (-k * ry * x1) / rx;
  const cx = cos * cxp - sin * cyp + (x0 + x) / 2;
  const cy = sin * cxp + cos * cyp + (y0 + y) / 2;
  const ang = (ux: number, uy: number, vx: number, vy: number) => {
    const a = Math.atan2(ux * vy - uy * vx, ux * vx + uy * vy);
    return a;
  };
  const t0 = ang(1, 0, (x1 - cxp) / rx, (y1 - cyp) / ry);
  let dt = ang((x1 - cxp) / rx, (y1 - cyp) / ry, (-x1 - cxp) / rx, (-y1 - cyp) / ry);
  if (!sweep && dt > 0) dt -= 2 * Math.PI;
  else if (sweep && dt < 0) dt += 2 * Math.PI;
  const n = Math.max(1, Math.ceil(Math.abs(dt) / (Math.PI / 2) - 1e-9));
  const step = dt / n;
  const e = (4 / 3) * Math.tan(step / 4);
  const at = (t: number) => {
    const ct = Math.cos(t);
    const st = Math.sin(t);
    return [
      cx + rx * ct * cos - ry * st * sin,
      cy + rx * ct * sin + ry * st * cos,
      -rx * st * cos - ry * ct * sin,
      -rx * st * sin + ry * ct * cos,
    ];
  };
  let t = t0;
  let [px, py, pdx, pdy] = at(t);
  for (let i = 0; i < n; i++) {
    t += step;
    const [qx, qy, qdx, qdy] = at(t);
    out.cubic(px + e * pdx, py + e * pdy, qx - e * qdx, qy - e * qdy, i === n - 1 ? x : qx, i === n - 1 ? y : qy);
    [px, py, pdx, pdy] = [qx, qy, qdx, qdy];
  }
}
