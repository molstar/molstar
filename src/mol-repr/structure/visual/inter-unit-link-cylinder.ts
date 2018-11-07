/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Link, Structure, StructureElement } from 'mol-model/structure';
import { ComplexVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { LinkCylinderProps, createLinkCylinderMesh, LinkIterator, LinkCylinderParams } from './util/link';
import { Vec3 } from 'mol-math/linear-algebra';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { ComplexMeshVisual, ComplexMeshParams } from '../complex-visual';
import { Interval } from 'mol-data/int';
import { BitFlags } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

async function createInterUnitLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: LinkCylinderProps, mesh?: Mesh) {
    const links = structure.links
    const { bondCount, bonds } = links

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
            return theme.size.size(location)
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const InterUnitLinkParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
}
export const DefaultInterUnitLinkProps = PD.getDefaultValues(InterUnitLinkParams)
export type InterUnitLinkProps = typeof DefaultInterUnitLinkProps

export function InterUnitLinkVisual(): ComplexVisual<InterUnitLinkProps> {
    return ComplexMeshVisual<InterUnitLinkProps>({
        defaultProps: DefaultInterUnitLinkProps,
        createGeometry: createInterUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromStructure,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: (state: VisualUpdateState, newProps: InterUnitLinkProps, currentProps: InterUnitLinkProps) => {
            state.createGeometry = newProps.radialSegments !== currentProps.radialSegments
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