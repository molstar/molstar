/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ShrakeRupleyContext, VdWLookup } from './common.js';
import type { RuntimeContext } from '@molstar/core/task';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Unit } from '@molstar/model/model/structure';

// TODO
// - iterate over units and elements
// - avoid using serial-element index whenever possible
// - calculate atomRadiusType only for invariant units
// - factor serialResidueIndex out

const updateChunk = 5000;
export async function computeArea(runtime: RuntimeContext, ctx: ShrakeRupleyContext) {
  const { atomRadiusType: atomRadius } = ctx;
  const spherePoints = flattenSpherePoints(ctx.spherePoints);
  const neighbors = createNeighbors(256);
  for (let i = 0; i < atomRadius.length; i += updateChunk) {
    if (runtime.shouldUpdate) {
      await runtime.update({
        message: 'Computing per residue surface accessibility...',
        current: i,
        max: atomRadius.length,
      });
    }

    computeRange(ctx, spherePoints, neighbors, i, Math.min(i + updateChunk, atomRadius.length));
  }
}

function flattenSpherePoints(points: Vec3[]) {
  const flat = new Float64Array(3 * points.length);
  for (let i = 0; i < points.length; ++i) {
    flat[3 * i] = points[i][0];
    flat[3 * i + 1] = points[i][1];
    flat[3 * i + 2] = points[i][2];
  }
  return flat;
}

/** Neighbours of the current atom: the projection threshold and the offset from the atom; reused for every atom. */
interface Neighbors {
  sqRadius: Float64Array;
  x: Float64Array;
  y: Float64Array;
  z: Float64Array;
}

function createNeighbors(capacity: number): Neighbors {
  return {
    sqRadius: new Float64Array(capacity),
    x: new Float64Array(capacity),
    y: new Float64Array(capacity),
    z: new Float64Array(capacity),
  };
}

function ensureNeighborCapacity(neighbors: Neighbors, count: number) {
  if (count <= neighbors.sqRadius.length) return;
  const size = Math.max(count, 2 * neighbors.sqRadius.length);
  neighbors.sqRadius = new Float64Array(size);
  neighbors.x = new Float64Array(size);
  neighbors.y = new Float64Array(size);
  neighbors.z = new Float64Array(size);
}

function computeRange(
  ctx: ShrakeRupleyContext,
  spherePoints: Float64Array,
  neighbors: Neighbors,
  begin: number,
  end: number,
) {
  const { structure, atomRadiusType, serialResidueIndex, area, scalingConstant, maxLookupRadius, probeSize } = ctx;
  const { lookup3d, serialMapping, unitIndexMap, units } = structure;
  const { cumulativeUnitElementCount, elementIndices, unitIndices } = serialMapping;
  const missing = VdWLookup[0];
  const spherePointCount = spherePoints.length / 3;

  for (let aI = begin; aI < end; ++aI) {
    const vdw1 = VdWLookup[atomRadiusType[aI]];
    if (vdw1 === missing) continue;

    const aUnit = units[unitIndices[aI]];
    const aElementIndex = elementIndices[aI];
    const aX = aUnit.conformation.x(aElementIndex);
    const aY = aUnit.conformation.y(aElementIndex);
    const aZ = aUnit.conformation.z(aElementIndex);

    // pre-filter by lookup3d (provides >10x speed-up compared to naive evaluation)
    const { count, units: lUnits, indices, squaredDistances } = lookup3d.find(aX, aY, aZ, maxLookupRadius);

    ensureNeighborCapacity(neighbors, count);
    const { sqRadius, x: nX, y: nY, z: nZ } = neighbors;

    // see optimizations proposed in Eisenhaber et al., 1995 (https://doi.org/10.1002/jcc.540160303)
    // collect neighbors for each atom
    const radius1 = probeSize + vdw1;
    const cutoff1 = probeSize + radius1;
    let neighborCount = 0;
    // the lookup returns the hits of one unit together
    let bUnit: Unit | undefined;
    let bUnitOffset = 0;
    for (let iI = 0; iI < count; ++iI) {
      if (lUnits[iI] !== bUnit) {
        bUnit = lUnits[iI];
        bUnitOffset = cumulativeUnitElementCount[unitIndexMap.get(bUnit.id)];
      }
      const bI = bUnitOffset + indices[iI];

      const vdw2 = VdWLookup[atomRadiusType[bI]];
      if (bI === aI || vdw2 === missing) continue;

      const radius2 = probeSize + vdw2;
      if (squaredDistances[iI] < (cutoff1 + vdw2) * (cutoff1 + vdw2)) {
        const bElementIndex = elementIndices[bI];
        // while here: compute values for later lookup
        sqRadius[neighborCount] = (squaredDistances[iI] + radius1 * radius1 - radius2 * radius2) / (2 * radius1);
        nX[neighborCount] = bUnit.conformation.x(bElementIndex) - aX;
        nY[neighborCount] = bUnit.conformation.y(bElementIndex) - aY;
        nZ[neighborCount] = bUnit.conformation.z(bElementIndex) - aZ;
        ++neighborCount;
      }
    }

    // A point is buried if any neighbor occludes it, so the order of the neighbors does not change the count. Instead
    // of sorting them by distance, test the last occluder first: neighboring sphere points tend to share it
    // (Eisenhaber et al., 1995).
    let accessiblePointCount = 0;
    if (neighborCount === 0) {
      accessiblePointCount = spherePointCount;
    } else {
      let lastOccluder = 0;
      sl: for (let sI = 0; sI < spherePointCount; ++sI) {
        const sX = spherePoints[3 * sI];
        const sY = spherePoints[3 * sI + 1];
        const sZ = spherePoints[3 * sI + 2];
        if (sX * nX[lastOccluder] + sY * nY[lastOccluder] + sZ * nZ[lastOccluder] > sqRadius[lastOccluder]) continue;
        for (let nI = 0; nI < neighborCount; ++nI) {
          const dot = sX * nX[nI] + sY * nY[nI] + sZ * nZ[nI];
          if (dot > sqRadius[nI]) {
            lastOccluder = nI;
            continue sl;
          }
        }
        ++accessiblePointCount;
      }
    }

    area[serialResidueIndex[aI]] += scalingConstant * accessiblePointCount * radius1 * radius1;
  }
}
