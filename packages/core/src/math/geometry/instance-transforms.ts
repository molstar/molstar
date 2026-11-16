/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../linear-algebra/3d/mat4.js';
import type { Vec3 } from '../linear-algebra/3d/vec3.js';
import { GridPositions, type Positions } from './grid-positions.js';

/**
 * The transforms of a set of instances as their readers take them: a matrix an instance, or a
 * pose an instance. Readers that go over all instances (bounds, grids, culling cells) branch on `kind` once and run a
 * loop of their own for each, here as the kernels below; none reads them a frame
 */
export type InstanceTransforms = InstanceMatrices | InstancePoses;

/** A matrix an instance: 16 floats, column-major, rigid with a uniform scale (as instance transforms are) */
export interface InstanceMatrices {
  readonly kind: 'matrices';
  readonly count: number;
  readonly matrices: Float32Array;
  /** the first float of the first instance's matrix */
  readonly offset?: number;
}

/**
 * A pose an instance: the instance's element (a particle) placed by its position, rotation (a unit quaternion) and
 * scale, around `center` (subtracted first), then `matrix` (the representation's, rigid with a uniform scale):
 * world = matrix · (position + scale · rotate(rotation, local - center)). Instance `i`'s element is `order[i]` (its
 * id, in the grid's order) mapped by `indices` (or `first` + id where the elements are a range)
 */
export interface InstancePoses {
  readonly kind: 'poses';
  readonly count: number;
  /** x, y, z an element */
  /** x, y, z an element: float32 or on a grid (`GridPositions`, decoded where read) */
  readonly positions: Positions;
  /**
   * x, y, z, w an element, a quaternion: unit float32, or snorm16 (`Int16Array`, a unit component times 32767); read
   * as it is, the rotation by `2 / |q|²` in place of 2 (exact for either); none: no rotation
   */
  readonly rotations?: Float32Array | Int16Array;
  /** an element's scale, `scale` without; one not above 0 is read as 1 (a particle without a radius) */
  readonly scales?: Float32Array;
  /** every element's scale where `scales` is none (a target's particles of one radius), 1 without; not above 0 read as 1 */
  readonly scale?: number;
  /** the instances' ids in their order (the grid's), `i` without */
  readonly order?: ArrayLike<number>;
  /** an id's element, `first + id` without */
  readonly indices?: ArrayLike<number>;
  readonly first?: number;
  readonly center: Vec3;
  readonly matrix?: Mat4;
}

/** Instance `i`'s element (`InstancePoses`) */
function elementOf(t: InstancePoses, i: number) {
  const id = t.order ? t.order[i] : i;
  return t.indices ? t.indices[id] : (t.first ?? 0) + id;
}

/** The scale of every element where `scales` is none (`InstancePoses.scale`) */
function uniformScaleOf(t: InstancePoses) {
  const s = t.scale;
  return s !== undefined && s > 0 ? s : 1;
}

/** Element `e`'s scale (`InstancePoses.scales`, `scale`) */
function scaleOf(t: InstancePoses, e: number) {
  if (!t.scales) return uniformScaleOf(t);
  const s = t.scales[e];
  return s > 0 ? s : 1;
}

const tmp = [0, 0, 0];
const pos = [0, 0, 0];

/**
 * Pose `e`'s placement of `local` without `matrix` into `out` (x, y, z): position + scale · rotate(q, local - center),
 * the rotation `v + f w (q × v) + f q × (q × v)` with `f = 2 / |q|²` (`rotations`), no matrix built
 */
function placePose(t: InstancePoses, e: number, lx: number, ly: number, lz: number, out: number[]) {
  const { positions: p, rotations: r, scales, center } = t;
  const vx = lx - center[0];
  const vy = ly - center[1];
  const vz = lz - center[2];
  const s = scales ? scaleOf(t, e) : uniformScaleOf(t);
  GridPositions.get(p, e, pos);
  if (!r) {
    out[0] = pos[0] + s * vx;
    out[1] = pos[1] + s * vy;
    out[2] = pos[2] + s * vz;
    return;
  }
  const q = e * 4;
  const qx = r[q],
    qy = r[q + 1],
    qz = r[q + 2],
    qw = r[q + 3];
  const f = rotationFactor(qx, qy, qz, qw);
  // t = f (q × v), v' = v + w t + q × t
  const tx = f * (qy * vz - qz * vy);
  const ty = f * (qz * vx - qx * vz);
  const tz = f * (qx * vy - qy * vx);
  out[0] = pos[0] + s * (vx + qw * tx + (qy * tz - qz * ty));
  out[1] = pos[1] + s * (vy + qw * ty + (qz * tx - qx * tz));
  out[2] = pos[2] + s * (vz + qw * tz + (qx * ty - qy * tx));
}

/** A quaternion's rotation factor, `2 / |q|²` (2 for a unit one; 0 for none, the identity) */
function rotationFactor(x: number, y: number, z: number, w: number) {
  const n = x * x + y * y + z * z + w * w;
  return n > 0 ? 2 / n : 0;
}

/** `out` (x, y, z) by `m` (affine) in place */
function applyMatrix(m: Mat4, out: number[]) {
  const x = out[0],
    y = out[1],
    z = out[2];
  out[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
  out[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
  out[2] = m[2] * x + m[6] * y + m[10] * z + m[14];
}

export namespace InstanceTransforms {
  export function ofMatrices(matrices: Float32Array, count: number, offset = 0): InstanceMatrices {
    return { kind: 'matrices', count, matrices, offset };
  }

  /**
   * Poses with every field set (`undefined` where not given), in one order: one shape for all poses objects, so that
   * the readers' property loads stay monomorphic. Make and change poses by it (`{ ...poses, order }` keeps the shape)
   */
  export function ofPoses(p: Omit<InstancePoses, 'kind'>): InstancePoses {
    return {
      kind: 'poses',
      count: p.count,
      positions: p.positions,
      rotations: p.rotations,
      scales: p.scales,
      scale: p.scale,
      order: p.order,
      indices: p.indices,
      first: p.first,
      center: p.center,
      matrix: p.matrix,
    };
  }

  /** Instance `i`'s placement of `local` into `out` */
  export function getPoint(t: InstanceTransforms, i: number, local: Vec3, out: Vec3): Vec3 {
    if (t.kind === 'matrices') {
      const m = t.matrices;
      const o = (t.offset ?? 0) + i * 16;
      const [x, y, z] = local;
      out[0] = m[o] * x + m[o + 4] * y + m[o + 8] * z + m[o + 12];
      out[1] = m[o + 1] * x + m[o + 5] * y + m[o + 9] * z + m[o + 13];
      out[2] = m[o + 2] * x + m[o + 6] * y + m[o + 10] * z + m[o + 14];
      return out;
    }
    placePose(t, elementOf(t, i), local[0], local[1], local[2], tmp);
    if (t.matrix) applyMatrix(t.matrix, tmp);
    out[0] = tmp[0];
    out[1] = tmp[1];
    out[2] = tmp[2];
    return out;
  }

  /** Instance `i`'s direction of `local` (a normal: the linear part, normalized) into `out` */
  export function getDirection(t: InstanceTransforms, i: number, local: Vec3, out: Vec3): Vec3 {
    const [x, y, z] = local;
    let ox: number;
    let oy: number;
    let oz: number;
    if (t.kind === 'matrices') {
      const m = t.matrices;
      const o = (t.offset ?? 0) + i * 16;
      ox = m[o] * x + m[o + 4] * y + m[o + 8] * z;
      oy = m[o + 1] * x + m[o + 5] * y + m[o + 9] * z;
      oz = m[o + 2] * x + m[o + 6] * y + m[o + 10] * z;
    } else {
      // rotated by q alone: the scale is uniform, the center a translation
      const r = t.rotations;
      const q = elementOf(t, i) * 4;
      const qx = r ? r[q] : 0,
        qy = r ? r[q + 1] : 0,
        qz = r ? r[q + 2] : 0,
        qw = r ? r[q + 3] : 1;
      const f = rotationFactor(qx, qy, qz, qw);
      const tx = f * (qy * z - qz * y);
      const ty = f * (qz * x - qx * z);
      const tz = f * (qx * y - qy * x);
      ox = x + qw * tx + (qy * tz - qz * ty);
      oy = y + qw * ty + (qz * tx - qx * tz);
      oz = z + qw * tz + (qx * ty - qy * tx);
      const m = t.matrix;
      if (m) {
        const ax = ox,
          ay = oy;
        ox = m[0] * ax + m[4] * ay + m[8] * oz;
        oy = m[1] * ax + m[5] * ay + m[9] * oz;
        oz = m[2] * ax + m[6] * ay + m[10] * oz;
      }
    }
    const len = Math.hypot(ox, oy, oz);
    const k = len > 0 ? 1 / len : 1;
    out[0] = ox * k;
    out[1] = oy * k;
    out[2] = oz * k;
    return out;
  }

  /** Instance `i`'s matrix into `out` (column-major): what a reader of matrices takes (export, few instances) */
  export function getMatrix(t: InstanceTransforms, i: number, out: Mat4): Mat4 {
    if (t.kind === 'matrices') return Mat4.fromArray(out, t.matrices, (t.offset ?? 0) + i * 16);
    const e = elementOf(t, i);
    const r = t.rotations;
    const q = e * 4;
    const x = r ? r[q] : 0,
      y = r ? r[q + 1] : 0,
      z = r ? r[q + 2] : 0,
      w = r ? r[q + 3] : 1;
    const s = scaleOf(t, e);
    const f = rotationFactor(x, y, z, w);
    // the rotation's columns times the scale (`Mat4.fromQuat`, by `f`)
    out[0] = s * (1 - f * (y * y + z * z));
    out[1] = s * f * (x * y + w * z);
    out[2] = s * f * (x * z - w * y);
    out[3] = 0;
    out[4] = s * f * (x * y - w * z);
    out[5] = s * (1 - f * (x * x + z * z));
    out[6] = s * f * (y * z + w * x);
    out[7] = 0;
    out[8] = s * f * (x * z + w * y);
    out[9] = s * f * (y * z - w * x);
    out[10] = s * (1 - f * (x * x + y * y));
    out[11] = 0;
    // the translation: position - (scale · rotation) · center
    GridPositions.get(t.positions, e, pos);
    const [cx, cy, cz] = t.center;
    out[12] = pos[0] - (out[0] * cx + out[4] * cy + out[8] * cz);
    out[13] = pos[1] - (out[1] * cx + out[5] * cy + out[9] * cz);
    out[14] = pos[2] - (out[2] * cx + out[6] * cy + out[10] * cz);
    out[15] = 1;
    return t.matrix ? Mat4.mul(out, t.matrix, out) : out;
  }

  /**
   * The first three rows of the matrices of instances `first` to `first + count` (12 floats each: an affine matrix
   * without its last row, as the GPU takes them) into `out` from `outOffset`, an instance each `outStride` floats
   */
  export function writeAffineRows(
    t: InstanceTransforms,
    first: number,
    count: number,
    out: Float32Array,
    outOffset = 0,
    outStride = 12,
  ) {
    if (t.kind === 'matrices') {
      const m = t.matrices;
      // unrolled: millions of instances (mesoscale)
      for (let i = 0, s = (t.offset ?? 0) + first * 16, o = outOffset; i < count; ++i, s += 16, o += outStride) {
        out[o] = m[s];
        out[o + 1] = m[s + 4];
        out[o + 2] = m[s + 8];
        out[o + 3] = m[s + 12];
        out[o + 4] = m[s + 1];
        out[o + 5] = m[s + 5];
        out[o + 6] = m[s + 9];
        out[o + 7] = m[s + 13];
        out[o + 8] = m[s + 2];
        out[o + 9] = m[s + 6];
        out[o + 10] = m[s + 10];
        out[o + 11] = m[s + 14];
      }
      return out;
    }
    // inlined as `writePoints`: A = [s R | position - s R center], then the matrix's M A
    const { positions: p, rotations: r, scales, order, indices, matrix: m } = t;
    // the positions' kind once: float32, else the grid's, decoded inline (`GridPositions`)
    const pf = p instanceof Float32Array ? p : undefined;
    const g = p instanceof Float32Array ? undefined : p;
    const gd = g ? g.data : pf!;
    const gs = g ? g.step : 0;
    const [gox, goy, goz] = g ? g.origin : [0, 0, 0];
    const uniform = uniformScaleOf(t);
    const elementFirst = t.first ?? 0;
    const [cx, cy, cz] = t.center;
    for (let k = 0, o = outOffset; k < count; ++k, o += outStride) {
      const i = first + k;
      const id = order ? order[i] : i;
      const e = indices ? indices[id] : elementFirst + id;
      let sc = scales ? scales[e] : uniform;
      if (!(sc > 0)) sc = 1;
      let x = 0,
        y = 0,
        z = 0,
        w = 1;
      if (r) {
        const q = e * 4;
        x = r[q];
        y = r[q + 1];
        z = r[q + 2];
        w = r[q + 3];
      }
      // a_ij: row i, column j of s R (`Mat4.fromQuat`, by `f = 2 / |q|²`)
      const nq = x * x + y * y + z * z + w * w;
      const f = nq > 0 ? 2 / nq : 0;
      const a00 = sc * (1 - f * (y * y + z * z)),
        a10 = sc * f * (x * y + w * z),
        a20 = sc * f * (x * z - w * y);
      const a01 = sc * f * (x * y - w * z),
        a11 = sc * (1 - f * (x * x + z * z)),
        a21 = sc * f * (y * z + w * x);
      const a02 = sc * f * (x * z + w * y),
        a12 = sc * f * (y * z - w * x),
        a22 = sc * (1 - f * (x * x + y * y));
      let ex: number, ey: number, ez: number;
      if (pf) {
        ex = pf[e * 3];
        ey = pf[e * 3 + 1];
        ez = pf[e * 3 + 2];
      } else {
        const a = gd[e * 2],
          b = gd[e * 2 + 1];
        ex = gox + gs * (a & 0x1fffff);
        ey = goy + gs * ((a >>> 21) | ((b & 0x3ff) << 11));
        ez = goz + gs * ((b >>> 10) & 0x1fffff);
      }
      const t0 = ex - (a00 * cx + a01 * cy + a02 * cz);
      const t1 = ey - (a10 * cx + a11 * cy + a12 * cz);
      const t2 = ez - (a20 * cx + a21 * cy + a22 * cz);
      if (!m) {
        out[o] = a00;
        out[o + 1] = a01;
        out[o + 2] = a02;
        out[o + 3] = t0;
        out[o + 4] = a10;
        out[o + 5] = a11;
        out[o + 6] = a12;
        out[o + 7] = t1;
        out[o + 8] = a20;
        out[o + 9] = a21;
        out[o + 10] = a22;
        out[o + 11] = t2;
        continue;
      }
      // M A: row i of M (m[i], m[4 + i], m[8 + i], m[12 + i]) times A's columns and translation
      for (let ri = 0; ri < 3; ++ri) {
        const m0 = m[ri],
          m1 = m[4 + ri],
          m2 = m[8 + ri],
          m3 = m[12 + ri];
        const ro = o + ri * 4;
        out[ro] = m0 * a00 + m1 * a10 + m2 * a20;
        out[ro + 1] = m0 * a01 + m1 * a11 + m2 * a21;
        out[ro + 2] = m0 * a02 + m1 * a12 + m2 * a22;
        out[ro + 3] = m0 * t0 + m1 * t1 + m2 * t2 + m3;
      }
    }
    return out;
  }

  /** Instance `i`'s largest scale on an axis (rigid instances: their uniform scale) */
  export function getMaxScale(t: InstanceTransforms, i: number): number {
    if (t.kind === 'matrices') {
      const m = t.matrices;
      const o = (t.offset ?? 0) + i * 16;
      return Math.sqrt(
        Math.max(
          m[o] * m[o] + m[o + 1] * m[o + 1] + m[o + 2] * m[o + 2],
          m[o + 4] * m[o + 4] + m[o + 5] * m[o + 5] + m[o + 6] * m[o + 6],
          m[o + 8] * m[o + 8] + m[o + 9] * m[o + 9] + m[o + 10] * m[o + 10],
        ),
      );
    }
    const s = scaleOf(t, elementOf(t, i));
    return t.matrix ? s * Mat4.getMaxScaleOnAxis(t.matrix) : s;
  }

  /** The largest scale on an axis of all instances (0 without instances) */
  export function getMaxScaleOfAll(t: InstanceTransforms): number {
    if (t.kind === 'matrices') {
      let max = 0;
      const m = t.matrices;
      for (let i = 0, o = t.offset ?? 0; i < t.count; ++i, o += 16) {
        max = Math.max(
          max,
          m[o] * m[o] + m[o + 1] * m[o + 1] + m[o + 2] * m[o + 2],
          m[o + 4] * m[o + 4] + m[o + 5] * m[o + 5] + m[o + 6] * m[o + 6],
          m[o + 8] * m[o + 8] + m[o + 9] * m[o + 9] + m[o + 10] * m[o + 10],
        );
      }
      return Math.sqrt(max);
    }
    if (t.count === 0) return 0;
    let max = uniformScaleOf(t);
    if (t.scales) {
      max = 0;
      // inlined as `writePoints`: the fields read once
      const { scales, order, indices } = t;
      const first = t.first ?? 0;
      for (let i = 0; i < t.count; ++i) {
        const id = order ? order[i] : i;
        const s = scales[indices ? indices[id] : first + id];
        const v = s > 0 ? s : 1;
        if (v > max) max = v;
      }
    }
    return t.matrix ? max * Mat4.getMaxScaleOnAxis(t.matrix) : max;
  }

  /**
   * The placements of `local` by every instance into `x`, `y`, `z` (from 0), returning their box (`min`, `max`): what
   * the grids, cells and bounds start from
   */
  export function writePoints(
    t: InstanceTransforms,
    local: Vec3,
    x: Float32Array | undefined,
    y: Float32Array | undefined,
    z: Float32Array | undefined,
    min: Vec3,
    max: Vec3,
  ) {
    let minX = Infinity,
      minY = Infinity,
      minZ = Infinity;
    let maxX = -Infinity,
      maxY = -Infinity,
      maxZ = -Infinity;
    const [lx, ly, lz] = local;
    if (t.kind === 'matrices') {
      const m = t.matrices;
      for (let i = 0, o = t.offset ?? 0; i < t.count; ++i, o += 16) {
        const px = m[o] * lx + m[o + 4] * ly + m[o + 8] * lz + m[o + 12];
        const py = m[o + 1] * lx + m[o + 5] * ly + m[o + 9] * lz + m[o + 13];
        const pz = m[o + 2] * lx + m[o + 6] * ly + m[o + 10] * lz + m[o + 14];
        if (x) {
          x[i] = px;
          y![i] = py;
          z![i] = pz;
        }
        if (px < minX) minX = px;
        if (py < minY) minY = py;
        if (pz < minZ) minZ = pz;
        if (px > maxX) maxX = px;
        if (py > maxY) maxY = py;
        if (pz > maxZ) maxZ = pz;
      }
    } else {
      // inlined, the fields read once and `local - center` once (`placePose` called an instance, its fields read each
      // time from objects of several shapes, cost MG_1184's setup about 3%)
      const { positions: p, rotations: r, scales, order, indices, matrix: m } = t;
      // the positions' kind once: float32, else the grid's, decoded inline (`GridPositions`)
      const pf = p instanceof Float32Array ? p : undefined;
      const g = p instanceof Float32Array ? undefined : p;
      const gd = g ? g.data : pf!;
      const gs = g ? g.step : 0;
      const [gox, goy, goz] = g ? g.origin : [0, 0, 0];
      const uniform = uniformScaleOf(t);
      const first = t.first ?? 0;
      const vx = lx - t.center[0],
        vy = ly - t.center[1],
        vz = lz - t.center[2];
      for (let i = 0; i < t.count; ++i) {
        const id = order ? order[i] : i;
        const e = indices ? indices[id] : first + id;
        let sc = scales ? scales[e] : uniform;
        if (!(sc > 0)) sc = 1;
        let rx = vx,
          ry = vy,
          rz = vz;
        if (r) {
          const q = e * 4;
          const qx = r[q],
            qy = r[q + 1],
            qz = r[q + 2],
            qw = r[q + 3];
          // by `f = 2 / |q|²` (`rotations`)
          const nq = qx * qx + qy * qy + qz * qz + qw * qw;
          const f = nq > 0 ? 2 / nq : 0;
          const tx = f * (qy * vz - qz * vy);
          const ty = f * (qz * vx - qx * vz);
          const tz = f * (qx * vy - qy * vx);
          rx = vx + qw * tx + (qy * tz - qz * ty);
          ry = vy + qw * ty + (qz * tx - qx * tz);
          rz = vz + qw * tz + (qx * ty - qy * tx);
        }
        let px: number, py: number, pz: number;
        if (pf) {
          px = pf[e * 3];
          py = pf[e * 3 + 1];
          pz = pf[e * 3 + 2];
        } else {
          const a = gd[e * 2],
            b = gd[e * 2 + 1];
          px = gox + gs * (a & 0x1fffff);
          py = goy + gs * ((a >>> 21) | ((b & 0x3ff) << 11));
          pz = goz + gs * ((b >>> 10) & 0x1fffff);
        }
        px += sc * rx;
        py += sc * ry;
        pz += sc * rz;
        if (m) {
          const ax = px,
            ay = py;
          px = m[0] * ax + m[4] * ay + m[8] * pz + m[12];
          py = m[1] * ax + m[5] * ay + m[9] * pz + m[13];
          pz = m[2] * ax + m[6] * ay + m[10] * pz + m[14];
        }
        if (x) {
          x[i] = px;
          y![i] = py;
          z![i] = pz;
        }
        if (px < minX) minX = px;
        if (py < minY) minY = py;
        if (pz < minZ) minZ = pz;
        if (px > maxX) maxX = px;
        if (py > maxY) maxY = py;
        if (pz > maxZ) maxZ = pz;
      }
    }
    min[0] = minX;
    min[1] = minY;
    min[2] = minZ;
    max[0] = maxX;
    max[1] = maxY;
    max[2] = maxZ;
  }

  /** The largest squared distance of the placements of `local` by every instance from `from` */
  export function getMaxDistanceSq(t: InstanceTransforms, local: Vec3, from: Vec3): number {
    const [lx, ly, lz] = local;
    const [bx, by, bz] = from;
    let maxDistSq = 0;
    if (t.kind === 'matrices') {
      const m = t.matrices;
      for (let i = 0, o = t.offset ?? 0; i < t.count; ++i, o += 16) {
        const dx = m[o] * lx + m[o + 4] * ly + m[o + 8] * lz + m[o + 12] - bx;
        const dy = m[o + 1] * lx + m[o + 5] * ly + m[o + 9] * lz + m[o + 13] - by;
        const dz = m[o + 2] * lx + m[o + 6] * ly + m[o + 10] * lz + m[o + 14] - bz;
        const d = dx * dx + dy * dy + dz * dz;
        if (d > maxDistSq) maxDistSq = d;
      }
    } else {
      // inlined, the fields read once and `local - center` once (`placePose` called an instance, its fields read each
      // time from objects of several shapes, cost MG_1184's setup about 3%)
      const { positions: p, rotations: r, scales, order, indices, matrix: m } = t;
      // the positions' kind once: float32, else the grid's, decoded inline (`GridPositions`)
      const pf = p instanceof Float32Array ? p : undefined;
      const g = p instanceof Float32Array ? undefined : p;
      const gd = g ? g.data : pf!;
      const gs = g ? g.step : 0;
      const [gox, goy, goz] = g ? g.origin : [0, 0, 0];
      const uniform = uniformScaleOf(t);
      const first = t.first ?? 0;
      const vx = lx - t.center[0],
        vy = ly - t.center[1],
        vz = lz - t.center[2];
      for (let i = 0; i < t.count; ++i) {
        const id = order ? order[i] : i;
        const e = indices ? indices[id] : first + id;
        let sc = scales ? scales[e] : uniform;
        if (!(sc > 0)) sc = 1;
        let rx = vx,
          ry = vy,
          rz = vz;
        if (r) {
          const q = e * 4;
          const qx = r[q],
            qy = r[q + 1],
            qz = r[q + 2],
            qw = r[q + 3];
          // by `f = 2 / |q|²` (`rotations`)
          const nq = qx * qx + qy * qy + qz * qz + qw * qw;
          const f = nq > 0 ? 2 / nq : 0;
          const tx = f * (qy * vz - qz * vy);
          const ty = f * (qz * vx - qx * vz);
          const tz = f * (qx * vy - qy * vx);
          rx = vx + qw * tx + (qy * tz - qz * ty);
          ry = vy + qw * ty + (qz * tx - qx * tz);
          rz = vz + qw * tz + (qx * ty - qy * tx);
        }
        let px: number, py: number, pz: number;
        if (pf) {
          px = pf[e * 3];
          py = pf[e * 3 + 1];
          pz = pf[e * 3 + 2];
        } else {
          const a = gd[e * 2],
            b = gd[e * 2 + 1];
          px = gox + gs * (a & 0x1fffff);
          py = goy + gs * ((a >>> 21) | ((b & 0x3ff) << 11));
          pz = goz + gs * ((b >>> 10) & 0x1fffff);
        }
        px += sc * rx;
        py += sc * ry;
        pz += sc * rz;
        if (m) {
          const ax = px,
            ay = py;
          px = m[0] * ax + m[4] * ay + m[8] * pz + m[12];
          py = m[1] * ax + m[5] * ay + m[9] * pz + m[13];
          pz = m[2] * ax + m[6] * ay + m[10] * pz + m[14];
        }
        const dx = px - bx,
          dy = py - by,
          dz = pz - bz;
        const d = dx * dx + dy * dy + dz * dz;
        if (d > maxDistSq) maxDistSq = d;
      }
    }
    return maxDistSq;
  }
}
