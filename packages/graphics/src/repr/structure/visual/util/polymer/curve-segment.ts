/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { NumberArray } from '@molstar/core/util/type-helpers';
import { lerp, spline } from '@molstar/core/math/interpolate';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3set = Vec3.set;
const v3orthogonalize = Vec3.orthogonalize;

export interface CurveSegmentState {
  curvePoints: NumberArray;
  tangentVectors: NumberArray;
  normalVectors: NumberArray;
  binormalVectors: NumberArray;
  widthValues: NumberArray;
  heightValues: NumberArray;
  linearSegments: number;
}

export interface CurveSegmentControls {
  secStrucFirst: boolean;
  secStrucLast: boolean;
  p0: Vec3;
  p1: Vec3;
  p2: Vec3;
  p3: Vec3;
  p4: Vec3;
  d12: Vec3;
  d23: Vec3;
}

export function createCurveSegmentState(linearSegments: number): CurveSegmentState {
  const n = linearSegments + 1;
  const pn = n * 3;
  return {
    curvePoints: new Float32Array(pn),
    tangentVectors: new Float32Array(pn),
    normalVectors: new Float32Array(pn),
    binormalVectors: new Float32Array(pn),
    widthValues: new Float32Array(n),
    heightValues: new Float32Array(n),
    linearSegments,
  };
}

export function interpolateCurveSegment(
  state: CurveSegmentState,
  controls: CurveSegmentControls,
  tension: number,
  shift: number,
) {
  interpolatePointsAndTangents(state, controls, tension, shift);
  interpolateNormals(state, controls);
}

// The per-point loops below work on scalar locals instead of `Vec3` calls: the `Vec3` helpers are called app-wide
// with many array kinds, so their keyed loads are megamorphic, also where V8 inlines them. In
// `interpolatePointsAndTangents` each expression keeps the operation order of the `Vec3` function it replaces (named in
// the comments), so the output is bit-identical.

export function interpolatePointsAndTangents(
  state: CurveSegmentState,
  controls: CurveSegmentControls,
  tension: number,
  shift: number,
) {
  const { curvePoints, tangentVectors, linearSegments } = state;
  const { p0, p1, p2, p3, p4, secStrucFirst, secStrucLast } = controls;

  const p0x = p0[0],
    p0y = p0[1],
    p0z = p0[2];
  const p1x = p1[0],
    p1y = p1[1],
    p1z = p1[2];
  const p2x = p2[0],
    p2y = p2[1],
    p2z = p2[2];
  const p3x = p3[0],
    p3y = p3[1],
    p3z = p3[2];
  const p4x = p4[0],
    p4y = p4[1],
    p4z = p4[2];

  const shift1 = 1 - shift;

  const tensionBeg = secStrucFirst ? 0.5 : tension;
  const tensionEnd = secStrucLast ? 0.5 : tension;

  for (let j = 0; j <= linearSegments; ++j) {
    const t = (j * 1.0) / linearSegments;
    let cx: number, cy: number, cz: number, ax: number, ay: number, az: number, bx: number, by: number, bz: number;
    // `Vec3.spline` of the curve point and two points around it
    if (t < shift1) {
      const te = lerp(tensionBeg, tension, t);
      const tc = t + shift;
      cx = spline(p0x, p1x, p2x, p3x, tc, te);
      cy = spline(p0y, p1y, p2y, p3y, tc, te);
      cz = spline(p0z, p1z, p2z, p3z, tc, te);
      const ta = t + shift + 0.01;
      ax = spline(p0x, p1x, p2x, p3x, ta, tensionBeg);
      ay = spline(p0y, p1y, p2y, p3y, ta, tensionBeg);
      az = spline(p0z, p1z, p2z, p3z, ta, tensionBeg);
      const tb = t + shift - 0.01;
      bx = spline(p0x, p1x, p2x, p3x, tb, tensionBeg);
      by = spline(p0y, p1y, p2y, p3y, tb, tensionBeg);
      bz = spline(p0z, p1z, p2z, p3z, tb, tensionBeg);
    } else {
      const te = lerp(tension, tensionEnd, t);
      const tc = t - shift1;
      cx = spline(p1x, p2x, p3x, p4x, tc, te);
      cy = spline(p1y, p2y, p3y, p4y, tc, te);
      cz = spline(p1z, p2z, p3z, p4z, tc, te);
      const ta = t - shift1 + 0.01;
      ax = spline(p1x, p2x, p3x, p4x, ta, te);
      ay = spline(p1y, p2y, p3y, p4y, ta, te);
      az = spline(p1z, p2z, p3z, p4z, ta, te);
      const tb = t - shift1 - 0.01;
      bx = spline(p1x, p2x, p3x, p4x, tb, te);
      by = spline(p1y, p2y, p3y, p4y, tb, te);
      bz = spline(p1z, p2z, p3z, p4z, tb, te);
    }
    const j3 = j * 3;
    curvePoints[j3] = cx;
    curvePoints[j3 + 1] = cy;
    curvePoints[j3 + 2] = cz;
    // `Vec3.normalize(Vec3.sub(a, b))`
    let tx = ax - bx,
      ty = ay - by,
      tz = az - bz;
    let len = tx * tx + ty * ty + tz * tz;
    if (len > 0) {
      len = 1 / Math.sqrt(len);
      tx *= len;
      ty *= len;
      tz *= len;
    }
    tangentVectors[j3] = tx;
    tangentVectors[j3 + 1] = ty;
    tangentVectors[j3 + 2] = tz;
  }
}

const tmpTangent = Vec3();
const tmpDirection = Vec3();
const tmpNormal = Vec3();
const firstTangentVec = Vec3();
const lastTangentVec = Vec3();
const firstNormalVec = Vec3();
const lastNormalVec = Vec3();

const HalfPi = Math.PI / 2;

/**
 * Populate normalVectors and binormalVectors with a rotation minimizing frame, propagated from firstDirection by the
 * double reflection method (Wang et al. 2008). The residual twist to reach lastDirection, taken modulo 180° since the
 * profiles are 2-fold symmetric, is distributed evenly along the segment.
 */
export function interpolateNormals(state: CurveSegmentState, controls: CurveSegmentControls) {
  const { curvePoints, tangentVectors, normalVectors, binormalVectors } = state;
  const { d12: firstDirection, d23: lastDirection } = controls;

  const n = curvePoints.length / 3;
  const n1 = n - 1;

  v3fromArray(firstTangentVec, tangentVectors, 0);
  v3fromArray(lastTangentVec, tangentVectors, n1 * 3);

  v3orthogonalize(firstNormalVec, firstTangentVec, firstDirection);
  v3orthogonalize(lastNormalVec, lastTangentVec, lastDirection);

  let nx = firstNormalVec[0],
    ny = firstNormalVec[1],
    nz = firstNormalVec[2];
  normalVectors[0] = nx;
  normalVectors[1] = ny;
  normalVectors[2] = nz;

  let qx = curvePoints[0],
    qy = curvePoints[1],
    qz = curvePoints[2];
  let tx = tangentVectors[0],
    ty = tangentVectors[1],
    tz = tangentVectors[2];

  for (let i = 1; i < n; ++i) {
    const i3 = i * 3;
    const ux = tx,
      uy = ty,
      uz = tz;
    tx = tangentVectors[i3];
    ty = tangentVectors[i3 + 1];
    tz = tangentVectors[i3 + 2];
    const px = curvePoints[i3],
      py = curvePoints[i3 + 1],
      pz = curvePoints[i3 + 2];
    const sx = px - qx,
      sy = py - qy,
      sz = pz - qz;
    qx = px;
    qy = py;
    qz = pz;

    let rx = nx,
      ry = ny,
      rz = nz;
    const c1 = sx * sx + sy * sy + sz * sz;
    if (c1 > 1e-12) {
      const k1 = -2 / c1;
      const kr = k1 * (sx * rx + sy * ry + sz * rz);
      rx += sx * kr;
      ry += sy * kr;
      rz += sz * kr;
      const ku = k1 * (sx * ux + sy * uy + sz * uz);
      const vx = tx - (ux + sx * ku),
        vy = ty - (uy + sy * ku),
        vz = tz - (uz + sz * ku);
      const c2 = vx * vx + vy * vy + vz * vz;
      if (c2 > 1e-12) {
        const k2 = (-2 / c2) * (vx * rx + vy * ry + vz * rz);
        rx += vx * k2;
        ry += vy * k2;
        rz += vz * k2;
      }
    }

    const cx = ty * rz - tz * ry,
      cy = tz * rx - tx * rz,
      cz = tx * ry - ty * rx;
    nx = cy * tz - cz * ty;
    ny = cz * tx - cx * tz;
    nz = cx * ty - cy * tx;
    let nl = nx * nx + ny * ny + nz * nz;
    if (nl > 0) {
      nl = 1 / Math.sqrt(nl);
      nx *= nl;
      ny *= nl;
      nz *= nl;
    } else {
      v3set(tmpTangent, tx, ty, tz);
      v3set(tmpDirection, rx, ry, rz);
      v3orthogonalize(tmpNormal, tmpTangent, tmpDirection);
      nx = tmpNormal[0];
      ny = tmpNormal[1];
      nz = tmpNormal[2];
    }
    normalVectors[i3] = nx;
    normalVectors[i3 + 1] = ny;
    normalVectors[i3 + 2] = nz;
  }

  const lx = lastNormalVec[0],
    ly = lastNormalVec[1],
    lz = lastNormalVec[2];
  let twist = Math.atan2(
    (ty * nz - tz * ny) * lx + (tz * nx - tx * nz) * ly + (tx * ny - ty * nx) * lz,
    nx * lx + ny * ly + nz * lz,
  );
  if (twist > HalfPi) twist -= Math.PI;
  else if (twist < -HalfPi) twist += Math.PI;

  const cosStep = Math.cos(twist / n1);
  const sinStep = Math.sin(twist / n1);
  let cos = 1,
    sin = 0;
  for (let i = 0; i < n; ++i) {
    const i3 = i * 3;
    const tx = tangentVectors[i3],
      ty = tangentVectors[i3 + 1],
      tz = tangentVectors[i3 + 2];
    const nx = normalVectors[i3],
      ny = normalVectors[i3 + 1],
      nz = normalVectors[i3 + 2];
    let bx = ty * nz - tz * ny,
      by = tz * nx - tx * nz,
      bz = tx * ny - ty * nx;
    let bl = bx * bx + by * by + bz * bz;
    if (bl > 0) {
      bl = 1 / Math.sqrt(bl);
      bx *= bl;
      by *= bl;
      bz *= bl;
      normalVectors[i3] = nx * cos + bx * sin;
      normalVectors[i3 + 1] = ny * cos + by * sin;
      normalVectors[i3 + 2] = nz * cos + bz * sin;
      binormalVectors[i3] = bx * cos - nx * sin;
      binormalVectors[i3 + 1] = by * cos - ny * sin;
      binormalVectors[i3 + 2] = bz * cos - nz * sin;
    } else {
      binormalVectors[i3] = 0;
      binormalVectors[i3 + 1] = 0;
      binormalVectors[i3 + 2] = 0;
    }
    const c = cos * cosStep - sin * sinStep;
    sin = sin * cosStep + cos * sinStep;
    cos = c;
  }
}

export function interpolateSizes(
  state: CurveSegmentState,
  w0: number,
  w1: number,
  w2: number,
  h0: number,
  h1: number,
  h2: number,
  shift: number,
) {
  const { widthValues, heightValues, linearSegments } = state;

  const shift1 = 1 - shift;

  for (let i = 0; i <= linearSegments; ++i) {
    const t = (i * 1.0) / linearSegments;
    if (t < shift1) {
      widthValues[i] = lerp(w0, w1, t + shift);
      heightValues[i] = lerp(h0, h1, t + shift);
    } else {
      widthValues[i] = lerp(w1, w2, t - shift1);
      heightValues[i] = lerp(h1, h2, t - shift1);
    }
  }
}
