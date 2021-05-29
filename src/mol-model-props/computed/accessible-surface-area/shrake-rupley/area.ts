/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShrakeRupleyContext, VdWLookup } from './common';
import { RuntimeContext } from '../../../../mol-task';

// TODO
// - iterate over units and elements
// - avoid using serial-element index whenever possible
// - calculate atomRadiusType only for invariant units
// - factor serialResidueIndex out

const updateChunk = 5000;
export async function computeArea(runtime: RuntimeContext, ctx: ShrakeRupleyContext) {
    const { atomRadiusType: atomRadius } = ctx;
    for (let i = 0; i < atomRadius.length; i += updateChunk) {
        if (runtime.shouldUpdate) {
            await runtime.update({ message: 'Computing per residue surface accessibility...', current: i, max: atomRadius.length });
        }

        computeRange(ctx, i, Math.min(i + updateChunk, atomRadius.length));
    }
}

function computeRange(ctx: ShrakeRupleyContext, begin: number, end: number) {
    const { structure, atomRadiusType, serialResidueIndex, area, spherePoints, scalingConstant, maxLookupRadius, probeSize } = ctx;
    const { lookup3d, serialMapping, unitIndexMap, units } = structure;
    const { cumulativeUnitElementCount, elementIndices, unitIndices } = serialMapping;

    for (let aI = begin; aI < end; ++aI) {
        const vdw1 = VdWLookup[atomRadiusType[aI]];
        if (vdw1 === VdWLookup[0]) continue;

        const aUnit = units[unitIndices[aI]];
        const aElementIndex = elementIndices[aI];
        const aX = aUnit.conformation.x(aElementIndex);
        const aY = aUnit.conformation.y(aElementIndex);
        const aZ = aUnit.conformation.z(aElementIndex);

        // pre-filter by lookup3d (provides >10x speed-up compared to naive evaluation)
        const { count, units: lUnits, indices, squaredDistances } = lookup3d.find(aX, aY, aZ, maxLookupRadius);

        // see optimizations proposed in Eisenhaber et al., 1995 (https://doi.org/10.1002/jcc.540160303)
        // collect neighbors for each atom
        const radius1 = probeSize + vdw1;
        const cutoff1 = probeSize + radius1;
        const neighbors = []; // TODO reuse
        for (let iI = 0; iI < count; ++iI) {
            const bUnit = lUnits[iI];
            const bI = cumulativeUnitElementCount[unitIndexMap.get(bUnit.id)] + indices[iI];
            const bElementIndex = elementIndices[bI];

            const vdw2 = VdWLookup[atomRadiusType[bI]];
            if ((aUnit === bUnit && aElementIndex === bElementIndex) || vdw2 === VdWLookup[0]) continue;

            const radius2 = probeSize + vdw2;
            if (squaredDistances[iI] < (cutoff1 + vdw2) * (cutoff1 + vdw2)) {
                const bElementIndex = elementIndices[bI];
                // while here: compute values for later lookup
                neighbors[neighbors.length] = [squaredDistances[iI],
                    (squaredDistances[iI] + radius1 * radius1 - radius2 * radius2) / (2 * radius1),
                    bUnit.conformation.x(bElementIndex) - aX,
                    bUnit.conformation.y(bElementIndex) - aY,
                    bUnit.conformation.z(bElementIndex) - aZ];
            }
        }

        // sort ascendingly by distance for improved downstream performance
        neighbors.sort((a, b) => a[0] - b[0]);

        let accessiblePointCount = 0;
        sl: for (let sI = 0; sI < spherePoints.length; ++sI) {
            const [sX, sY, sZ] = spherePoints[sI];
            for (let nI = 0; nI < neighbors.length; ++nI) {
                const [, sqRadius, nX, nY, nZ] = neighbors[nI];
                const dot = sX * nX + sY * nY + sZ * nZ;
                if (dot > sqRadius) {
                    continue sl;
                }
            }
            ++accessiblePointCount;
        }

        area[serialResidueIndex[aI]] += scalingConstant * accessiblePointCount * radius1 * radius1;
    }
}