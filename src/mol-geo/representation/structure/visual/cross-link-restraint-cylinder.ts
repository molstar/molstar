/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { Link, Structure, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, ComplexVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { LinkCylinderProps, DefaultLinkCylinderProps, createLinkCylinderMesh } from './util/link';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, MarkerData } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { ComplexMeshVisual } from '../complex-visual';
import { LocationIterator } from './util/location-iterator';

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
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeTheme,
    flipSided: false,
    flatShaded: false,
}
export type CrossLinkRestraintProps = typeof DefaultCrossLinkRestraintProps

// TODO update & create in-place
// if (currentProps.radialSegments !== newProps.radialSegments) return false

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

function markLink(loci: Loci, action: MarkerAction, structure: Structure, values: MarkerData) {
    const tMarker = values.tMarker

    const crossLinks = structure.crossLinkRestraints
    const elementCount = crossLinks.count
    const instanceCount = 1

    let changed = false
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        applyMarkerAction(array, 0, elementCount * instanceCount, action)
        changed = true
    } else if (Link.isLoci(loci)) {
        for (const b of loci.links) {
            const indices = crossLinks.getPairIndices(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
            if (indices) {
                for (let i = 0, il = indices.length; i < il; ++i) {
                    const idx = indices[i]
                    if (applyMarkerAction(array, idx, idx + 1, action) && !changed) {
                        changed = true
                    }
                }
            }
        }
    } else {
        return
    }
    if (changed) {
        ValueCell.update(tMarker, tMarker.ref.value)
    }
}