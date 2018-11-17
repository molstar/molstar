/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link, Structure, StructureElement } from 'mol-model/structure';
import { ComplexVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { createLinkCylinderMesh, LinkIterator, LinkCylinderParams } from './util/link';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { ComplexMeshVisual, ComplexMeshParams } from '../complex-visual';
import { Interval } from 'mol-data/int';
import { BitFlags } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualContext } from 'mol-repr/representation';
import { Theme } from 'mol-theme/theme';

async function createInterUnitLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitLinkParams>, mesh?: Mesh) {
    const links = structure.links
    const { bondCount, bonds } = links
    const { sizeFactor } = props

    if (!bondCount) return Mesh.createEmpty(mesh)

    const location = StructureElement.create()

    const builderProps = {
        linkCount: bondCount,
        referencePosition: (edgeIndex: number) => null, // TODO
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = bonds[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
        },
        order: (edgeIndex: number) => bonds[edgeIndex].order,
        flags: (edgeIndex: number) => BitFlags.create(bonds[edgeIndex].flag),
        radius: (edgeIndex: number) => {
            const b = bonds[edgeIndex]
            location.unit = b.unitA
            location.element = b.unitA.elements[b.indexA]
            return theme.size.size(location) * sizeFactor
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const InterUnitLinkParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
}
export type InterUnitLinkParams = typeof InterUnitLinkParams

export function InterUnitLinkVisual(): ComplexVisual<InterUnitLinkParams> {
    return ComplexMeshVisual<InterUnitLinkParams>({
        defaultProps: PD.getDefaultValues(InterUnitLinkParams),
        createGeometry: createInterUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromStructure,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitLinkParams>, currentProps: PD.Values<InterUnitLinkParams>) => {
            if (newProps.linkScale !== currentProps.linkScale) state.createGeometry = true
            if (newProps.linkSpacing !== currentProps.linkSpacing) state.createGeometry = true
            if (newProps.radialSegments !== currentProps.radialSegments) state.createGeometry = true
        }
    })
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const bond = structure.links.bonds[groupId]
        return Link.Loci(structure, [
            Link.Location(
                bond.unitA, bond.indexA as StructureElement.UnitIndex,
                bond.unitB, bond.indexB as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function markLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!Link.isLoci(loci)) return false
    for (const b of loci.links) {
        const idx = structure.links.getBondIndex(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
        if (idx !== -1) {
            if (apply(Interval.ofSingleton(idx))) changed = true
        }
    }
    return changed
}