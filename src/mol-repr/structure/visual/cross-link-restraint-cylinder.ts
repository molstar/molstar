/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link, Structure, StructureElement } from 'mol-model/structure';
import { ComplexVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { LinkCylinderProps, createLinkCylinderMesh, LinkCylinderParams } from './util/link';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { ComplexMeshVisual, ComplexMeshParams } from '../complex-visual';
import { Interval } from 'mol-data/int';
import { SizeThemeOptions, SizeThemeName } from 'mol-theme/size';
import { BitFlags } from 'mol-util';
import { LinkType } from 'mol-model/structure/model/types';
import { SelectParam, NumberParam, paramDefaultValues } from 'mol-util/parameter';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

async function createCrossLinkRestraintCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: LinkCylinderProps, mesh?: Mesh) {

    const crossLinks = structure.crossLinkRestraints
    if (!crossLinks.count) return Mesh.createEmpty(mesh)

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
            return theme.size.size(location)
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const CrossLinkRestraintParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 1, 0, 20, 0.1),
}
export const DefaultCrossLinkRestraintProps = paramDefaultValues(CrossLinkRestraintParams)
export type CrossLinkRestraintProps = typeof DefaultCrossLinkRestraintProps

export function CrossLinkRestraintVisual(): ComplexVisual<CrossLinkRestraintProps> {
    return ComplexMeshVisual<CrossLinkRestraintProps>({
        defaultProps: DefaultCrossLinkRestraintProps,
        createGeometry: createCrossLinkRestraintCylinderMesh,
        createLocationIterator: CrossLinkRestraintIterator,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: (state: VisualUpdateState, newProps: CrossLinkRestraintProps, currentProps: CrossLinkRestraintProps) => {
            state.createGeometry = newProps.radialSegments !== currentProps.radialSegments
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
            return Link.Loci([ Link.Location(pair.unitA, pair.indexA, pair.unitB, pair.indexB) ])
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