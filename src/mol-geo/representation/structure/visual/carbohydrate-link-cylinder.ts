/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, Link, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, ComplexVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { DefaultMeshProps } from '../../util';
import { Vec3 } from 'mol-math/linear-algebra';
import { LocationIterator } from './util/location-iterator';
import { createLinkCylinderMesh, DefaultLinkCylinderProps, LinkCylinderProps } from './util/link';
import { OrderedSet, Interval } from 'mol-data/int';
import { ComplexMeshVisual } from '../complex-visual';

// TODO create seperate visual
// for (let i = 0, il = carbohydrates.terminalLinks.length; i < il; ++i) {
//     const tl = carbohydrates.terminalLinks[i]
//     const center = carbohydrates.elements[tl.carbohydrateIndex].geometry.center
//     tl.elementUnit.conformation.position(tl.elementUnit.elements[tl.elementIndex], p)
//     if (tl.fromCarbohydrate) {
//         builder.addCylinder(center, p, 0.5, linkParams)
//     } else {
//         builder.addCylinder(p, center, 0.5, linkParams)
//     }
// }

async function createCarbohydrateLinkCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {
    const { links, elements } = structure.carbohydrates

    const builderProps = {
        linkCount: links.length,
        referencePosition: (edgeIndex: number) => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const l = links[edgeIndex]
            Vec3.copy(posA, elements[l.carbohydrateIndexA].geometry.center)
            Vec3.copy(posB, elements[l.carbohydrateIndexB].geometry.center)
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => 0
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultCarbohydrateLinkProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type CarbohydrateLinkProps = typeof DefaultCarbohydrateLinkProps

export function CarbohydrateLinkVisual(): ComplexVisual<CarbohydrateLinkProps> {
    return ComplexMeshVisual<CarbohydrateLinkProps>({
        defaultProps: DefaultCarbohydrateLinkProps,
        createMesh: createCarbohydrateLinkCylinderMesh,
        createLocationIterator: CarbohydrateLinkIterator,
        getLoci: getLinkLoci,
        mark: markLink
    })
}

function CarbohydrateLinkIterator(structure: Structure): LocationIterator {
    const { elements, links } = structure.carbohydrates
    const elementCount = links.length
    const instanceCount = 1
    const location = Link.Location()
    const getLocation = (elementIndex: number, instanceIndex: number) => {
        const link = links[elementIndex]
        const carbA = elements[link.carbohydrateIndexA]
        const carbB = elements[link.carbohydrateIndexB]
        const indexA = OrderedSet.findPredecessorIndex(carbA.unit.elements, carbA.anomericCarbon)
        const indexB = OrderedSet.findPredecessorIndex(carbB.unit.elements, carbB.anomericCarbon)
        location.aUnit = carbA.unit
        location.aIndex = indexA as StructureElement.UnitIndex
        location.bUnit = carbB.unit
        location.bIndex = indexB as StructureElement.UnitIndex
        return location
    }
    return LocationIterator(elementCount, instanceCount, getLocation)
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const { links, elements } = structure.carbohydrates
        const l = links[elementId]
        const carbA = elements[l.carbohydrateIndexA]
        const carbB = elements[l.carbohydrateIndexB]
        const indexA = OrderedSet.findPredecessorIndex(carbA.unit.elements, carbA.anomericCarbon)
        const indexB = OrderedSet.findPredecessorIndex(carbB.unit.elements, carbB.anomericCarbon)
        return Link.Loci([
            Link.Location(
                carbA.unit, indexA as StructureElement.UnitIndex,
                carbB.unit, indexB as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function markLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getLinkIndex } = structure.carbohydrates

    let changed = false
    if (Link.isLoci(loci)) {
        for (const l of loci.links) {
            const idx = getLinkIndex(l.aUnit, l.aUnit.elements[l.aIndex], l.bUnit, l.bUnit.elements[l.bIndex])
            if (idx !== undefined) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    }
    return changed
}