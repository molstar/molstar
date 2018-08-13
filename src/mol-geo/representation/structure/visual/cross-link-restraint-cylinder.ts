/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link, Structure, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, ComplexVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { LinkCylinderProps, DefaultLinkCylinderProps, createLinkCylinderMesh } from './util/link';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { ComplexMeshVisual } from '../complex-visual';
import { LocationIterator } from './util/location-iterator';
import { Interval } from 'mol-data/int';
import { SizeThemeProps } from 'mol-view/theme/size';

async function createCrossLinkRestraintCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {

    const crossLinks = structure.crossLinkRestraints
    if (!crossLinks.count) return Mesh.createEmpty(mesh)

    const builderProps = {
        linkCount: crossLinks.count,
        referencePosition: (edgeIndex: number) => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = crossLinks.pairs[edgeIndex]
            // console.log(b)
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
            // console.log(posA, posB)
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => 0
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultCrossLinkRestraintProps = {
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeThemeProps,
    flipSided: false,
    flatShaded: false,
}
export type CrossLinkRestraintProps = typeof DefaultCrossLinkRestraintProps

export function CrossLinkRestraintVisual(): ComplexVisual<CrossLinkRestraintProps> {
    return ComplexMeshVisual<CrossLinkRestraintProps>({
        defaultProps: DefaultCrossLinkRestraintProps,
        createMesh: createCrossLinkRestraintCylinderMesh,
        createLocationIterator: CrossLinkRestraintIterator,
        getLoci: getLinkLoci,
        mark: markLink
    })
}

function CrossLinkRestraintIterator(structure: Structure): LocationIterator {
    const { pairs } = structure.crossLinkRestraints
    const elementCount = pairs.length
    const instanceCount = 1
    const location = Link.Location()
    const getLocation = (elementIndex: number, instanceIndex: number) => {
        const pair = pairs[elementIndex]
        location.aUnit = pair.unitA
        location.aIndex = pair.indexA
        location.bUnit = pair.unitB
        location.bIndex = pair.indexB
        return location
    }
    return LocationIterator(elementCount, instanceCount, getLocation)
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const pair = structure.crossLinkRestraints.pairs[elementId]
        if (pair) {
            return Link.Loci([
                Link.Location(
                    pair.unitA, pair.indexA as StructureElement.UnitIndex,
                    pair.unitB, pair.indexB as StructureElement.UnitIndex
                )
            ])
        }
    }
    return EmptyLoci
}

function markLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const crossLinks = structure.crossLinkRestraints

    let changed = false
    if (Link.isLoci(loci)) {
        for (const b of loci.links) {
            const indices = crossLinks.getPairIndices(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
            if (indices) {
                for (let i = 0, il = indices.length; i < il; ++i) {
                    if (apply(Interval.ofSingleton(indices[i]))) changed = true
                }
            }
        }
    }
    return changed
}