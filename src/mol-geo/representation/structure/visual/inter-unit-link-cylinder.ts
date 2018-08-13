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
import { SizeTheme } from '../../../theme';
import { LinkIterator } from './util/location-iterator';
import { ComplexMeshVisual } from '../complex-visual';
import { Interval } from 'mol-data/int';

async function createInterUnitLinkCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {
    const links = structure.links
    const { bondCount, bonds } = links

    if (!bondCount) return Mesh.createEmpty(mesh)

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
        flags: (edgeIndex: number) => bonds[edgeIndex].flag
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultInterUnitLinkProps = {
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeTheme,
}
export type InterUnitLinkProps = typeof DefaultInterUnitLinkProps

export function InterUnitLinkVisual(): ComplexVisual<InterUnitLinkProps> {
    return ComplexMeshVisual<InterUnitLinkProps>({
        defaultProps: DefaultInterUnitLinkProps,
        createMesh: createInterUnitLinkCylinderMesh,
        createLocationIterator: LinkIterator.fromStructure,
        getLoci: getLinkLoci,
        mark: markLink
    })
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const bond = structure.links.bonds[elementId]
        return Link.Loci([
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
    if (Link.isLoci(loci)) {
        for (const b of loci.links) {
            const idx = structure.links.getBondIndex(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    }
    return changed
}