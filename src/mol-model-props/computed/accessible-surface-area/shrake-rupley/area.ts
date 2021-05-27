/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShrakeRupleyContext, VdWLookup } from './common';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { RuntimeContext } from '../../../../mol-task';
import { StructureProperties, StructureElement, Structure } from '../../../../mol-model/structure/structure';

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

const aPos = Vec3();

function setLocation(l: StructureElement.Location, structure: Structure, serialIndex: number) {
    l.structure = structure;
    l.unit = structure.units[structure.serialMapping.unitIndices[serialIndex]];
    l.element = structure.serialMapping.elementIndices[serialIndex];
    return l;
}

function computeRange(ctx: ShrakeRupleyContext, begin: number, end: number) {
    const { structure, atomRadiusType, serialResidueIndex, area, spherePoints, scalingConstant, maxLookupRadius, probeSize } = ctx;
    const aLoc = StructureElement.Location.create(structure);
    const bLoc = StructureElement.Location.create(structure);
    const { x, y, z } = StructureProperties.atom;
    const { lookup3d, serialMapping, unitIndexMap } = structure;
    const { cumulativeUnitElementCount } = serialMapping;

    for (let aI = begin; aI < end; ++aI) {
        const vdw1 = VdWLookup[atomRadiusType[aI]];
        if (vdw1 === VdWLookup[0]) continue;

        setLocation(aLoc, structure, aI);
        const aX = x(aLoc);
        const aY = y(aLoc);
        const aZ = z(aLoc);
        Vec3.set(aPos, x(aLoc), y(aLoc), z(aLoc));

        // pre-filter by lookup3d (provides >10x speed-up compared to naive evaluation)
        const { count, units, indices, squaredDistances } = lookup3d.find(aX, aY, aZ, maxLookupRadius);

        // see optimizations proposed in Eisenhaber et al., 1995 (https://doi.org/10.1002/jcc.540160303)
        // collect neighbors for each atom
        const radius1 = probeSize + vdw1;
        const cutoff1 = probeSize + radius1;
        let neighbors = []; // TODO reuse
        for (let iI = 0; iI < count; ++iI) {
            const bUnit = units[iI];
            StructureElement.Location.set(bLoc, ctx.structure, bUnit, bUnit.elements[indices[iI]]);
            const bI = cumulativeUnitElementCount[unitIndexMap.get(bUnit.id)] + indices[iI];

            const vdw2 = VdWLookup[atomRadiusType[bI]];
            if (StructureElement.Location.areEqual(aLoc, bLoc) || vdw2 === VdWLookup[0]) continue;

            const cutoff2 = (cutoff1 + vdw2) * (cutoff1 + vdw2);
            const radius2 = probeSize + vdw2;
            if (squaredDistances[iI] < cutoff2) {
                setLocation(bLoc, structure, bI);
                // while here: compute values for later lookup
                neighbors[neighbors.length] = [squaredDistances[iI],
                    (squaredDistances[iI] + radius1 * radius1 - radius2 * radius2) / (2 * radius1),
                    x(bLoc) - aPos[0],
                    y(bLoc) - aPos[1],
                    z(bLoc) - aPos[2]];
            }
        }

        // sort ascendingly by distance for improved downstream performance
        neighbors = neighbors.sort((a, b) => a[0] - b[0]);

        let accessiblePointCount = 0;
        sl: for (let sI = 0; sI < spherePoints.length; ++sI) {
            const [sx, sy, sz] = spherePoints[sI];
            for (let nI = 0; nI < neighbors.length; ++nI) {
                const [, sqRadius, nx, ny, nz] = neighbors[nI];
                const dot = sx * nx + sy * ny + sz * nz;
                if (dot > sqRadius) {
                    continue sl;
                }
            }
            ++accessiblePointCount;
        }

        area[serialResidueIndex[aI]] += scalingConstant * accessiblePointCount * radius1 * radius1;
    }
}