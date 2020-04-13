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
const bPos = Vec3();
const testPoint = Vec3();
const aLoc = StructureElement.Location.create(void 0 as any);
const bLoc = StructureElement.Location.create(void 0 as any);

function setLocation(l: StructureElement.Location, structure: Structure, serialIndex: number) {
    l.structure = structure;
    l.unit = structure.units[structure.serialMapping.unitIndices[serialIndex]];
    l.element = structure.serialMapping.elementIndices[serialIndex];
    return l;
}

function computeRange(ctx: ShrakeRupleyContext, begin: number, end: number) {
    const { structure, atomRadiusType, serialResidueIndex, area, spherePoints, scalingConstant, maxLookupRadius, probeSize } = ctx;
    const { x, y, z } = StructureProperties.atom;
    const { lookup3d, serialMapping, unitIndexMap } = structure;
    const { cumulativeUnitElementCount } = serialMapping;

    for (let aI = begin; aI < end; ++aI) {
        const radius1 = VdWLookup[atomRadiusType[aI]];
        if (radius1 === VdWLookup[0]) continue;

        setLocation(aLoc, structure, aI);
        const aX = x(aLoc);
        const aY = y(aLoc);
        const aZ = z(aLoc);

        // pre-filter by lookup3d (provides >10x speed-up compared to naive evaluation)
        const { count, units, indices, squaredDistances } = lookup3d.find(aX, aY, aZ, maxLookupRadius);

        // collect neighbors for each atom
        const cutoff1 = probeSize + probeSize + radius1;
        const neighbors = []; // TODO reuse
        for (let iI = 0; iI < count; ++iI) {
            const bUnit = units[iI];
            StructureElement.Location.set(bLoc, ctx.structure, bUnit, bUnit.elements[indices[iI]]);
            const bI = cumulativeUnitElementCount[unitIndexMap.get(bUnit.id)] + indices[iI];

            const radius2 = VdWLookup[atomRadiusType[bI]];
            if (StructureElement.Location.areEqual(aLoc, bLoc) || radius2 === VdWLookup[0]) continue;

            const cutoff2 = (cutoff1 + radius2) * (cutoff1 + radius2);
            if (squaredDistances[iI] < cutoff2) {
                neighbors[neighbors.length] = bI;
            }
        }

        // for all neighbors: test all sphere points
        Vec3.set(aPos, aX, aY, aZ);
        const scale = probeSize + radius1;
        let accessiblePointCount = 0;
        for (let sI = 0; sI < spherePoints.length; ++sI) {
            Vec3.scaleAndAdd(testPoint, aPos, spherePoints[sI], scale);
            let accessible = true;

            for (let _nI = 0; _nI < neighbors.length; ++_nI) {
                const nI = neighbors[_nI];
                setLocation(bLoc, structure, nI);
                Vec3.set(bPos, x(bLoc), y(bLoc), z(bLoc));
                const radius3 = VdWLookup[atomRadiusType[nI]];
                const cutoff3 = (radius3 + probeSize) * (radius3 + probeSize);
                if (Vec3.squaredDistance(testPoint, bPos) < cutoff3) {
                    accessible = false;
                    break;
                }
            }

            if (accessible) ++accessiblePointCount;
        }

        area[serialResidueIndex[aI]] += scalingConstant * accessiblePointCount * scale * scale;
    }
}