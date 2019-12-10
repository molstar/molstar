/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ShrakeRupleyContext, VdWLookup } from './common';
import { Vec3 } from '../../../../mol-math/linear-algebra';

export async function computePerResidue(ctx: ShrakeRupleyContext) {
    const { rtctx, updateChunk, atomRadius } = ctx;
    for (let i = 0; i < atomRadius.length; i += updateChunk) {
        computeRange(ctx, i, Math.min(i + updateChunk, atomRadius.length));

        if (rtctx.shouldUpdate) {
            rtctx.update({ message: 'Computing per residue surface accessibility...', current: i, max: atomRadius.length });
        }
    }
    // async yields 4x speed-up
}

const aPos = Vec3.zero();
const bPos = Vec3.zero();
let testPoint = Vec3.zero();
function computeRange(ctx: ShrakeRupleyContext, begin: number, end: number) {
    const { structure, atomRadius, accessibleSurfaceArea, spherePoints, cons, maxLookupRadius, probeSize } = ctx;
        const { model } = structure;
        const { x, y, z } = model.atomicConformation;
        const { residueAtomSegments } = model.atomicHierarchy;
        const { lookup3d } = structure;

    const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i]);

    // console.log(`computing ASA chunk ${ begin } to ${ end }`);
    for (let aI = begin; aI < end; aI++) {
        const radius1 = VdWLookup[atomRadius[aI]];
        if (radius1 === VdWLookup[0]) continue;

        // pre-filter by lookup3d
        // lookup provides >10x speed-up compared to naive evaluation
        const { count, units, indices, squaredDistances } = lookup3d.find(x[aI], y[aI], z[aI], maxLookupRadius);
        // we could keep track of all found neighbors of each atom, however slow and memory intensive

        // collect neighbors for each atom
        const cutoff1 = probeSize + probeSize + radius1;
        const neighbors = [];
        for (let iI = 0; iI < count; ++iI) {
            const bI = units[iI].elements[indices[iI]];

            const radius2 = VdWLookup[atomRadius[bI]];
            if (aI === bI || radius2 === VdWLookup[0]) continue;

            const cutoff2 = (cutoff1 + radius2) * (cutoff1 + radius2);
            if (squaredDistances[iI] < cutoff2) {
                neighbors[neighbors.length] = bI;
            }
        }

        // for all neighbors: test all sphere points
        position(aI, aPos);
        const scalar = probeSize + radius1;
        let accessiblePointCount = 0;
        for (let sI = 0; sI < spherePoints.length; ++sI) {
            const spherePoint = spherePoints[sI];
            testPoint = [spherePoint[0] * scalar + aPos[0],
                spherePoint[1] * scalar + aPos[1],
                spherePoint[2] * scalar + aPos[2]] as Vec3;
            let accessible = true;

            for (let _nI = 0; _nI < neighbors.length; ++_nI) {
                const nI = neighbors[_nI];
                position(nI, bPos);
                const radius3 = VdWLookup[atomRadius[nI]];
                const cutoff3 = (radius3 + probeSize) * (radius3 + probeSize);
                if (Vec3.squaredDistance(testPoint, bPos) < cutoff3) {
                    accessible = false;
                    break;
                }
            }

            if (accessible) ++accessiblePointCount;
        }

        accessibleSurfaceArea[residueAtomSegments.index[aI]] += cons * accessiblePointCount * scalar * scalar;
    }
}