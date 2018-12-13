/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link, Structure, StructureElement } from 'mol-model/structure';
import { ComplexVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { createLinkCylinderMesh, LinkCylinderParams } from './util/link';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { ComplexMeshVisual, ComplexMeshParams } from '../complex-visual';
import { Interval } from 'mol-data/int';
import { BitFlags } from 'mol-util';
import { LinkType } from 'mol-model/structure/model/types';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';

function createCrossLinkRestraintCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CrossLinkRestraintParams>, mesh?: Mesh) {

    const crossLinks = structure.crossLinkRestraints
    if (!crossLinks.count) return Mesh.createEmpty(mesh)
    const { sizeFactor } = props

    const location = StructureElement.create()

    const builderProps = {
        linkCount: crossLinks.count,
        referencePosition: () => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = crossLinks.pairs[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
        },
        order: () => 1,
        flags: () => BitFlags.create(LinkType.Flag.None),
        radius: (edgeIndex: number) => {
            const b = crossLinks.pairs[edgeIndex]
            location.unit = b.unitA
            location.element = b.unitA.elements[b.indexA]
            return theme.size.size(location) * sizeFactor
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const CrossLinkRestraintParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
}
export type CrossLinkRestraintParams = typeof CrossLinkRestraintParams

export function CrossLinkRestraintVisual(): ComplexVisual<CrossLinkRestraintParams> {
    return ComplexMeshVisual<CrossLinkRestraintParams>({
        defaultProps: PD.getDefaultValues(CrossLinkRestraintParams),
        createGeometry: createCrossLinkRestraintCylinderMesh,
        createLocationIterator: CrossLinkRestraintIterator,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CrossLinkRestraintParams>, currentProps: PD.Values<CrossLinkRestraintParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            )
        }
    })
}

function CrossLinkRestraintIterator(structure: Structure): LocationIterator {
    const { pairs } = structure.crossLinkRestraints
    const groupCount = pairs.length
    const instanceCount = 1
    const location = Link.Location()
    const getLocation = (groupIndex: number) => {
        const pair = pairs[groupIndex]
        location.aUnit = pair.unitA
        location.aIndex = pair.indexA
        location.bUnit = pair.unitB
        location.bIndex = pair.indexB
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true)
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const pair = structure.crossLinkRestraints.pairs[groupId]
        if (pair) {
            return Link.Loci(structure, [
                Link.Location(pair.unitA, pair.indexA, pair.unitB, pair.indexB),
                Link.Location(pair.unitB, pair.indexB, pair.unitA, pair.indexA)
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