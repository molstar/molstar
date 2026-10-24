/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Fred Ludlow <fred.ludlow@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * ported from NGL (https://github.com/arose/ngl), licensed under MIT
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RuntimeContext } from '@molstar/core/task';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { fillGridDim, type PositionData } from './common.js';
import type { Boundary } from './boundary.js';
import { GridLookup3D } from './lookup3d/grid.js';
import { Box3D } from './primitives/box3d.js';
import { Vec3 } from '../linear-algebra/3d/vec3.js';
import { Tensor } from '../linear-algebra/tensor.js';
import { Mat4 } from '../linear-algebra/3d/mat4.js';

type AnglesTables = { cosTable: Float32Array; sinTable: Float32Array };
function getAngleTables(probePositions: number): AnglesTables {
  let theta = 0.0;
  const step = (2 * Math.PI) / probePositions;

  const cosTable = new Float32Array(probePositions);
  const sinTable = new Float32Array(probePositions);
  for (let i = 0; i < probePositions; i++) {
    cosTable[i] = Math.cos(theta);
    sinTable[i] = Math.sin(theta);
    theta += step;
  }
  return { cosTable, sinTable };
}

//

export const MolecularSurfaceCalculationParams = {
  probeRadius: PD.Numeric(
    1.4,
    { min: 0, max: 10, step: 0.1 },
    { description: 'Radius of the probe tracing the molecular surface.' },
  ),
  resolution: PD.Numeric(0.5, { min: 0.01, max: 20, step: 0.01 }, { description: 'Grid resolution/cell spacing.' }),
  probePositions: PD.Numeric(
    36,
    { min: 12, max: 90, step: 1 },
    { description: 'Number of positions tested for probe target intersection.' },
  ),
};
export const DefaultMolecularSurfaceCalculationProps = PD.getDefaultValues(MolecularSurfaceCalculationParams);
export type MolecularSurfaceCalculationProps = typeof DefaultMolecularSurfaceCalculationProps;

export async function calcMolecularSurface(
  ctx: RuntimeContext,
  position: Required<PositionData>,
  boundary: Boundary,
  maxRadius: number,
  box: Box3D | null,
  props: MolecularSurfaceCalculationProps,
) {
  // Field generation method adapted from AstexViewer (Mike Hartshorn) by Fred Ludlow.
  // Other parts based heavily on NGL (Alexander Rose) EDT Surface class

  // Occlusion only needs a yes/no, so it scans compact copies, not `neighbours` through `OrderedSet.getAt`
  // and scattered loads. An atom that obscures a point on the sphere of atom `a` is within `r_a + r_c` of
  // `a`, so in its `find` result (up to the last bit at exact tangency, where the old cross-atom `lastClip`
  // could still test it); the scan order and the `lastClip` slot never change the outcome.

  /** Neighbours of the current atom, without the atom itself: resolved index, ordinal, centre, radius, radius². */
  let nbCapacity = 0;
  let nbCount = 0;
  let nbIndex = new Int32Array(0);
  let nbOrdinal = new Int32Array(0);
  let nbX = new Float64Array(0);
  let nbY = new Float64Array(0);
  let nbZ = new Float64Array(0);
  let nbR = new Float64Array(0);
  let nbRSq = new Float64Array(0);
  /** Neighbours that can reach the probe circle of the current torus (`projectTorus` keeps the count): centre, radius². */
  let tX = new Float64Array(0);
  let tY = new Float64Array(0);
  let tZ = new Float64Array(0);
  let tRSq = new Float64Array(0);
  /** Slot of the last obscuring atom in the list being scanned (very often the same in subsequent calls). */
  let lastClip = -1;

  function growNeighbours(count: number) {
    nbCapacity = Math.max(64, 2 * count);
    nbIndex = new Int32Array(nbCapacity);
    nbOrdinal = new Int32Array(nbCapacity);
    nbX = new Float64Array(nbCapacity);
    nbY = new Float64Array(nbCapacity);
    nbZ = new Float64Array(nbCapacity);
    nbR = new Float64Array(nbCapacity);
    nbRSq = new Float64Array(nbCapacity);
    tX = new Float64Array(nbCapacity);
    tY = new Float64Array(nbCapacity);
    tZ = new Float64Array(nbCapacity);
    tRSq = new Float64Array(nbCapacity);
  }

  /**
   * Copy the current `lookup3d.find` result, skipping the query atom `self` (a resolved index).
   */
  function loadNeighbours(self: number) {
    const { count, indices: ordinals } = neighbours;
    if (count > nbCapacity) growNeighbours(count);
    let c = 0;
    for (let j = 0; j < count; ++j) {
      const o = ordinals[j];
      const ai = OrderedSet.getAt(indices, o);
      if (ai === self) continue;
      const r = radius[ai];
      nbIndex[c] = ai;
      nbOrdinal[c] = o;
      nbX[c] = px[ai];
      nbY[c] = py[ai];
      nbZ[c] = pz[ai];
      nbR[c] = Math.abs(r); // the occlusion test squares r
      nbRSq[c] = r * r;
      ++c;
    }
    nbCount = c;
    lastClip = -1;
  }

  /**
   * Is the point at x,y,z inside any of the `count` spheres (centres `cx/cy/cz`, squared radii `cRSq`)?
   */
  function obscured(
    x: number,
    y: number,
    z: number,
    cx: Float64Array,
    cy: Float64Array,
    cz: Float64Array,
    cRSq: Float64Array,
    count: number,
  ) {
    if (lastClip !== -1) {
      const j = lastClip;
      const dx = cx[j] - x;
      const dy = cy[j] - y;
      const dz = cz[j] - z;
      if (dx * dx + dy * dy + dz * dz < cRSq[j]) return true;
      lastClip = -1;
    }

    for (let j = 0; j < count; ++j) {
      const dx = cx[j] - x;
      const dy = cy[j] - y;
      const dz = cz[j] - z;
      if (dx * dx + dy * dy + dz * dz < cRSq[j]) {
        lastClip = j;
        return true;
      }
    }

    return false;
  }

  /**
   * Is the projection `(dx, dy, dz) * ap + (vx, vy, vz)` onto the current atom's sphere obscured by a neighbour?
   */
  function projectionObscured(vx: number, vy: number, vz: number, dx: number, dy: number, dz: number, ap: number) {
    return obscured(dx * ap + vx, dy * ap + vy, dz * ap + vz, nbX, nbY, nbZ, nbRSq, nbCount);
  }

  /**
   * For each atom:
   *     Iterate over a subsection of the grid, for each point:
   *         If current value < 0.0, unvisited, set positive
   *
   *         If the delta distance would lower the grid value (raise it outside, in extended mode):
   *             Project this point onto surface of the atomic sphere
   *             If this projected point is not obscured by any other atom
   *                 set grid value to delta
   */
  function projectPointsRange(begI: number, endI: number) {
    for (let i = begI; i < endI; ++i) {
      const j = OrderedSet.getAt(indices, i);
      const vx = px[j],
        vy = py[j],
        vz = pz[j];
      const rad = radius[j];
      const rSq = rad * rad;
      const extended = ngPoints > 0;

      lookup3d.find(vx, vy, vz, rad);
      loadNeighbours(j);

      // Number of grid points, round this up...
      const ng = Math.ceil(rad * scaleFactor) + ngPoints;

      // Center of the atom, mapped to grid points (take floor)
      const iax = Math.floor(scaleFactor * (vx - minX));
      const iay = Math.floor(scaleFactor * (vy - minY));
      const iaz = Math.floor(scaleFactor * (vz - minZ));

      // Extents of grid to consider for this atom
      const begX = Math.max(0, iax - ng);
      const begY = Math.max(0, iay - ng);
      const begZ = Math.max(0, iaz - ng);

      // Add two to these points:
      // - iax are floor'd values so this ensures coverage
      // - these are loop limits (exclusive)
      const endX = Math.min(dimX, iax + ng + 2);
      const endY = Math.min(dimY, iay + ng + 2);
      const endZ = Math.min(dimZ, iaz + ng + 2);

      for (let xi = begX; xi < endX; ++xi) {
        const dx = gridx[xi] - vx;
        const xIdx = xi * iuv;
        for (let yi = begY; yi < endY; ++yi) {
          const dy = gridy[yi] - vy;
          const dxySq = dx * dx + dy * dy;
          const xyIdx = yi * iu + xIdx;
          for (let zi = begZ; zi < endZ; ++zi) {
            const dz = gridz[zi] - vz;
            const dSq = dxySq + dz * dz;

            if (dSq < rSq) {
              const idx = zi + xyIdx;
              // if unvisited, make positive
              if (data[idx] < 0.0) data[idx] *= -1;

              // Occlusion is tested only when the voxel would change (the outcome is the same).
              // Project on to the surface of the sphere:
              // sp is the projected point ( dx, dy, dz ) * ( ra / d )
              const d = Math.sqrt(dSq);
              const dd = rad - d;
              if (dd < data[idx] && !projectionObscured(vx, vy, vz, dx, dy, dz, rad / d)) {
                data[idx] = dd;
                idData[idx] = id[i];
              }
            } else if (extended) {
              const idx = zi + xyIdx;
              const d = Math.sqrt(dSq);
              const dd = rad - d;
              if (dd > data[idx] && !projectionObscured(vx, vy, vz, dx, dy, dz, rad / d)) data[idx] = dd;
            }
          }
        }
      }
    }
  }

  async function projectPoints() {
    for (let i = 0; i < n; i += updateChunk) {
      projectPointsRange(i, Math.min(i + updateChunk, n));

      if (ctx.shouldUpdate) {
        await ctx.update({ message: 'projecting points', current: i, max: n });
      }
    }
  }

  /**
   * `a` and `b` must be resolved indices, `idA` and `idB` their ids; the neighbours of `a` are loaded.
   *
   * Scalar locals for the a→b axis, the two circle normals and the centre, computed as `Vec3.normalize`, `cross` and
   * `scale` compute them (same operations, same order: the same field to the bit), and a branchless gather of the
   * neighbours that can reach the probe circle (perf-audit-2 02/P11: the Vec3 helpers boxed ~24 numbers per torus once
   * their call sites went generic, and the gather mispredicted on about half of ~44 neighbours).
   */
  function projectTorus(a: number, b: number, idA: number, idB: number) {
    const rA = radius[a];
    const rB = radius[b];
    const dx = px[b] - px[a];
    const dy = py[b] - py[a];
    const dz = pz[b] - pz[a];
    const dSq = dx * dx + dy * dy + dz * dz;

    // This check now redundant as already done in AVHash.withinRadii
    // if (dSq > ((rA + rB) * (rA + rB))) { return }

    const d = Math.sqrt(dSq);

    // Find angle between a->b vector and the circle
    // of their intersection by cosine rule
    const cosA = (rA * rA + d * d - rB * rB) / (2.0 * rA * d);

    // distance along a->b at intersection
    const dmp = rA * cosA;

    // the unit a->b vector
    let ax = dx,
      ay = dy,
      az = dz;
    let len = ax * ax + ay * ay + az * az;
    if (len > 0) {
      len = 1 / Math.sqrt(len);
      ax *= len;
      ay *= len;
      az *= len;
    }

    // Create normal to line
    let n1x = 1.0,
      n1y = 1.0,
      n1z = 1.0;
    if (ax !== 0) {
      n1x = (ay + az) / -ax;
    } else if (ay !== 0) {
      n1y = (ax + az) / -ay;
    } else if (az !== 0) {
      n1z = (ax + ay) / -az;
    }
    len = n1x * n1x + n1y * n1y + n1z * n1z;
    if (len > 0) {
      len = 1 / Math.sqrt(len);
      n1x *= len;
      n1y *= len;
      n1z *= len;
    }

    // Cross together for second normal vector
    let n2x = ay * n1z - az * n1y,
      n2y = az * n1x - ax * n1z,
      n2z = ax * n1y - ay * n1x;
    len = n2x * n2x + n2y * n2y + n2z * n2z;
    if (len > 0) {
      len = 1 / Math.sqrt(len);
      n2x *= len;
      n2y *= len;
      n2z *= len;
    }

    // r is radius of circle of intersection
    const rInt = Math.sqrt(rA * rA - dmp * dmp);

    n1x *= rInt;
    n1y *= rInt;
    n1z *= rInt;
    n2x *= rInt;
    n2y *= rInt;
    n2z *= rInt;
    // a->b scaled to the circle's centre
    ax *= dmp;
    ay *= dmp;
    az *= dmp;

    const mx = ax + px[a],
      my = ay + py[a],
      mz = az + pz[a];

    // Probe positions lie on the circle (centre `mid`, radius `rInt`), so only neighbours with
    // |c - mid| < rInt + r_c can obscure one. The margin covers rounding, mostly the Float32 cos/sin
    // tables (|p - mid| <= rInt * (1 + 6e-8) + coordinate rounding). A NaN `rInt` (no intersection) keeps none, as before
    // nothing is obscured and no voxel is written.
    let tc = 0;
    for (let j = 0; j < nbCount; ++j) {
      const cx = nbX[j],
        cy = nbY[j],
        cz = nbZ[j];
      const cdx = cx - mx;
      const cdy = cy - my;
      const cdz = cz - mz;
      const reach = rInt * (1 + 1e-6) + nbR[j] + 1e-4;
      // written always, kept by advancing the count
      tX[tc] = cx;
      tY[tc] = cy;
      tZ[tc] = cz;
      tRSq[tc] = nbRSq[j];
      tc += cdx * cdx + cdy * cdy + cdz * cdz < reach * reach && nbIndex[j] !== b ? 1 : 0;
    }
    lastClip = -1;

    for (let i = 0; i < probePositions; ++i) {
      const cost = cosTable[i];
      const sint = sinTable[i];

      const px = mx + cost * n1x + sint * n2x;
      const py = my + cost * n1y + sint * n2y;
      const pz = mz + cost * n1z + sint * n2z;

      if (!obscured(px, py, pz, tX, tY, tZ, tRSq, tc)) {
        const iax = Math.floor(scaleFactor * (px - minX));
        const iay = Math.floor(scaleFactor * (py - minY));
        const iaz = Math.floor(scaleFactor * (pz - minZ));

        const begX = Math.max(0, iax - ngTorus);
        const begY = Math.max(0, iay - ngTorus);
        const begZ = Math.max(0, iaz - ngTorus);

        const endX = Math.min(dimX, iax + ngTorus + 2);
        const endY = Math.min(dimY, iay + ngTorus + 2);
        const endZ = Math.min(dimZ, iaz + ngTorus + 2);

        for (let xi = begX; xi < endX; ++xi) {
          const dx = px - gridx[xi];
          const xIdx = xi * iuv;

          for (let yi = begY; yi < endY; ++yi) {
            const dy = py - gridy[yi];
            const dxySq = dx * dx + dy * dy;
            const xyIdx = yi * iu + xIdx;

            for (let zi = begZ; zi < endZ; ++zi) {
              const dz = pz - gridz[zi];
              const dSq = dxySq + dz * dz;

              const idx = zi + xyIdx;
              const current = data[idx];

              if (current > 0.0 && dSq < current * current) {
                data[idx] = Math.sqrt(dSq);
                // Is this grid point closer to a or b?
                // Take dot product of atob and gridpoint->p (dx, dy, dz)
                const dp = dx * ax + dy * ay + dz * az;
                idData[idx] = dp < 0.0 ? idB : idA;
              }
            }
          }
        }
      }
    }
  }

  function projectToriiRange(begI: number, endI: number) {
    for (let i = begI; i < endI; ++i) {
      const k = OrderedSet.getAt(indices, i);
      lookup3d.find(px[k], py[k], pz[k], radius[k]);
      loadNeighbours(k);
      const idK = id[i];
      for (let j = 0; j < nbCount; ++j) {
        const l = nbIndex[j];
        if (k < l) projectTorus(k, l, idK, id[nbOrdinal[j]]);
      }
    }
  }

  async function projectTorii() {
    for (let i = 0; i < n; i += updateChunk) {
      projectToriiRange(i, Math.min(i + updateChunk, n));

      if (ctx.shouldUpdate) {
        await ctx.update({ message: 'projecting torii', current: i, max: n });
      }
    }
  }

  // console.time('MolecularSurface')
  // console.time('MolecularSurface createState')
  const { resolution, probeRadius, probePositions } = props;
  const scaleFactor = 1 / resolution;
  const ngTorus = 2 + Math.floor(probeRadius * scaleFactor);
  const ngPoints = probeRadius < resolution * 2 ? 1 : 0;

  const cellSize = Vec3.create(maxRadius, maxRadius, maxRadius);
  Vec3.scale(cellSize, cellSize, 2);
  const lookup3d = GridLookup3D(position, boundary, cellSize);
  const neighbours = lookup3d.result;
  if (box === null) box = lookup3d.boundary.box;

  const { indices, x: px, y: py, z: pz, id, radius } = position;
  const n = OrderedSet.size(indices);

  const pad = maxRadius + resolution;
  const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
  const [minX, minY, minZ] = expandedBox.min;
  const scaledBox = Box3D.scale(Box3D(), expandedBox, scaleFactor);
  const dim = Box3D.size(Vec3(), scaledBox);
  Vec3.ceil(dim, dim);

  const [dimX, dimY, dimZ] = dim;
  const iu = dimZ,
    iv = dimY,
    iuv = iu * iv;

  const { cosTable, sinTable } = getAngleTables(probePositions);

  const space = Tensor.Space(dim, [0, 1, 2], Float32Array);
  const data = space.create();
  const idData = space.create();

  data.fill(-1001.0);
  idData.fill(-1);

  const gridx = fillGridDim(dimX, minX, resolution);
  const gridy = fillGridDim(dimY, minY, resolution);
  const gridz = fillGridDim(dimZ, minZ, resolution);

  // atoms per yield by the work an atom costs: the voxels the points loop visits plus a fixed cost (`find`, the
  // neighbour load, the tori) in voxel-visit units, which dominates at coarse resolution. `1e5 / (maxRadius⁹ · scale)`
  // gave whole structures for small radii (perf-audit-2 02/B2); voxels alone gave seconds per chunk at res 8
  // (perf-audit-3 02/B6)
  const ngAtom = Math.ceil((maxRadius + probeRadius) * scaleFactor) + ngPoints;
  const updateChunk = Math.max(1, Math.ceil(2e6 / (Math.pow(2 * ngAtom + 2, 3) + 2000)));
  // console.timeEnd('MolecularSurface createState')

  // console.time('MolecularSurface projectPoints')
  await projectPoints();
  // console.timeEnd('MolecularSurface projectPoints')

  // console.time('MolecularSurface projectTorii')
  await projectTorii();
  // console.timeEnd('MolecularSurface projectTorii')
  // console.timeEnd('MolecularSurface')

  const field = Tensor.create(space, data);
  const idField = Tensor.create(space, idData);

  const transform = Mat4.identity();
  Mat4.fromScaling(transform, Vec3.create(resolution, resolution, resolution));
  Mat4.setTranslation(transform, expandedBox.min);
  // console.log({ field, idField, transform, updateChunk })
  return { field, idField, transform, resolution, maxRadius };
}
