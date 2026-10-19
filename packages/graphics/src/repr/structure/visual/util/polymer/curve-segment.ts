/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { NumberArray } from '@molstar/core/util/type-helpers';
import { clamp, lerp, smoothstep, spline } from '@molstar/core/math/interpolate';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3set = Vec3.set;
const v3orthogonalize = Vec3.orthogonalize;
const v3matchDirection = Vec3.matchDirection;

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
// with many array kinds, so their keyed loads are megamorphic, also where V8 inlines them. Each expression keeps the
// operation order of the `Vec3` function it replaces (named in the comments), so the output is bit-identical.

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

/**
 * Populate normalVectors by interpolating from firstDirection to lastDirection with
 * resulting vector perpendicular to tangentVectors and binormalVectors
 */
export function interpolateNormals(state: CurveSegmentState, controls: CurveSegmentControls) {
  const { curvePoints, tangentVectors, normalVectors, binormalVectors } = state;
  const { d12: firstDirection, d23: lastDirection } = controls;

  const n = curvePoints.length / 3;

  v3fromArray(firstTangentVec, tangentVectors, 0);
  v3fromArray(lastTangentVec, tangentVectors, (n - 1) * 3);

  v3orthogonalize(firstNormalVec, firstTangentVec, firstDirection);
  v3orthogonalize(lastNormalVec, lastTangentVec, lastDirection);
  v3matchDirection(lastNormalVec, lastNormalVec, firstNormalVec);

  const lx = lastNormalVec[0],
    ly = lastNormalVec[1],
    lz = lastNormalVec[2];
  // previous normal
  let px = firstNormalVec[0],
    py = firstNormalVec[1],
    pz = firstNormalVec[2];

  const n1 = n - 1;
  for (let i = 0; i < n; ++i) {
    const j = smoothstep(0, n1, i) * n1;
    const t = i === 0 ? 0 : 1 / (n - j);
    const i3 = i * 3;

    const tx = tangentVectors[i3],
      ty = tangentVectors[i3 + 1],
      tz = tangentVectors[i3 + 2];

    // `Vec3.slerp(prev, last, t)`
    const d = clamp(px * lx + py * ly + pz * lz, -1, 1);
    const theta = Math.acos(d) * t;
    const md = -d;
    let rx = lx + px * md,
      ry = ly + py * md,
      rz = lz + pz * md;
    let rl = rx * rx + ry * ry + rz * rz;
    if (rl > 0) {
      rl = 1 / Math.sqrt(rl);
      rx *= rl;
      ry *= rl;
      rz *= rl;
    }
    const cosTheta = Math.cos(theta);
    const sinTheta = Math.sin(theta);
    const sx = px * cosTheta + rx * sinTheta,
      sy = py * cosTheta + ry * sinTheta,
      sz = pz * cosTheta + rz * sinTheta;

    // `Vec3.orthogonalize(tangent, slerped)`: `normalize(cross(cross(tangent, slerped), tangent))`
    const cx = ty * sz - tz * sy,
      cy = tz * sx - tx * sz,
      cz = tx * sy - ty * sx;
    let nx = cy * tz - cz * ty,
      ny = cz * tx - cx * tz,
      nz = cx * ty - cy * tx;
    let nl = nx * nx + ny * ny + nz * nz;
    if (nl > 0) {
      nl = 1 / Math.sqrt(nl);
      nx *= nl;
      ny *= nl;
      nz *= nl;
    }
    if (nx === 0 && ny === 0 && nz === 0) {
      // the fallbacks for a parallel or zero input
      v3set(tmpTangent, tx, ty, tz);
      v3set(tmpDirection, sx, sy, sz);
      v3orthogonalize(tmpNormal, tmpTangent, tmpDirection);
      nx = tmpNormal[0];
      ny = tmpNormal[1];
      nz = tmpNormal[2];
    }
    normalVectors[i3] = nx;
    normalVectors[i3 + 1] = ny;
    normalVectors[i3 + 2] = nz;

    px = nx;
    py = ny;
    pz = nz;

    setBinormal(binormalVectors, i3, tx, ty, tz, nx, ny, nz);
  }

  for (let i = 1; i < n1; ++i) {
    const i3 = i * 3;
    // `Vec3.scale(Vec3.add(prev, Vec3.add(next, normal)), 1 / 3)`
    const nx = (normalVectors[i3 - 3] + (normalVectors[i3 + 3] + normalVectors[i3])) * (1 / 3);
    const ny = (normalVectors[i3 - 2] + (normalVectors[i3 + 4] + normalVectors[i3 + 1])) * (1 / 3);
    const nz = (normalVectors[i3 - 1] + (normalVectors[i3 + 5] + normalVectors[i3 + 2])) * (1 / 3);
    normalVectors[i3] = nx;
    normalVectors[i3 + 1] = ny;
    normalVectors[i3 + 2] = nz;

    setBinormal(binormalVectors, i3, tangentVectors[i3], tangentVectors[i3 + 1], tangentVectors[i3 + 2], nx, ny, nz);
  }
}

/** `Vec3.normalize(Vec3.cross(tangent, normal))` */
function setBinormal(
  binormalVectors: NumberArray,
  i3: number,
  tx: number,
  ty: number,
  tz: number,
  nx: number,
  ny: number,
  nz: number,
) {
  let bx = ty * nz - tz * ny,
    by = tz * nx - tx * nz,
    bz = tx * ny - ty * nx;
  let len = bx * bx + by * by + bz * bz;
  if (len > 0) {
    len = 1 / Math.sqrt(len);
    bx *= len;
    by *= len;
    bz *= len;
  }
  binormalVectors[i3] = bx;
  binormalVectors[i3 + 1] = by;
  binormalVectors[i3 + 2] = bz;
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
