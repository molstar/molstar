/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Link, StructureElement } from 'mol-model/structure';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { Vec3 } from 'mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderProps, LinkCylinderParams } from './util/link';
import { OrderedSet, Interval } from 'mol-data/int';
import { ComplexMeshVisual, ComplexVisual } from '../complex-visual';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { LinkType } from 'mol-model/structure/model/types';
import { BitFlags } from 'mol-util';
import { UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualUpdateState } from '../../util';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

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

const radiusFactor = 0.3

async function createCarbohydrateLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: LinkCylinderProps, mesh?: Mesh) {
    const { links, elements } = structure.carbohydrates
    const location = StructureElement.create()

    const builderProps = {
        linkCount: links.length,
        referencePosition: (edgeIndex: number) => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const l = links[edgeIndex]
            Vec3.copy(posA, elements[l.carbohydrateIndexA].geometry.center)
            Vec3.copy(posB, elements[l.carbohydrateIndexB].geometry.center)
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => BitFlags.create(LinkType.Flag.None),
        radius: (edgeIndex: number) => {
            const l = links[edgeIndex]
            location.unit = elements[l.carbohydrateIndexA].unit
            location.element = elements[l.carbohydrateIndexA].anomericCarbon
            return theme.size.size(location) * radiusFactor
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const CarbohydrateLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    sizeTheme: PD.Select<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
    sizeValue: PD.Numeric('Size Value', '', 1, 0, 20, 0.1),
    detail: PD.Numeric('Sphere Detail', '', 0, 0, 3, 1),
}
export const DefaultCarbohydrateLinkProps = PD.getDefaultValues(CarbohydrateLinkParams)
export type CarbohydrateLinkProps = typeof DefaultCarbohydrateLinkProps

export function CarbohydrateLinkVisual(): ComplexVisual<CarbohydrateLinkProps> {
    return ComplexMeshVisual<CarbohydrateLinkProps>({
        defaultProps: DefaultCarbohydrateLinkProps,
        createGeometry: createCarbohydrateLinkCylinderMesh,
        createLocationIterator: CarbohydrateLinkIterator,
        getLoci: getLinkLoci,
        mark: markLink,
        setUpdateState: (state: VisualUpdateState, newProps: CarbohydrateLinkProps, currentProps: CarbohydrateLinkProps) => {
            state.createGeometry = newProps.radialSegments !== currentProps.radialSegments
        }
    })
}

function CarbohydrateLinkIterator(structure: Structure): LocationIterator {
    const { elements, links } = structure.carbohydrates
    const groupCount = links.length
    const instanceCount = 1
    const location = Link.Location()
    const getLocation = (groupIndex: number) => {
        const link = links[groupIndex]
        const carbA = elements[link.carbohydrateIndexA]
        const carbB = elements[link.carbohydrateIndexB]
        const indexA = OrderedSet.indexOf(carbA.unit.elements, carbA.anomericCarbon)
        const indexB = OrderedSet.indexOf(carbB.unit.elements, carbB.anomericCarbon)
        location.aUnit = carbA.unit
        location.aIndex = indexA as StructureElement.UnitIndex
        location.bUnit = carbB.unit
        location.bIndex = indexB as StructureElement.UnitIndex
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true)
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const { links, elements } = structure.carbohydrates
        const l = links[groupId]
        const carbA = elements[l.carbohydrateIndexA]
        const carbB = elements[l.carbohydrateIndexB]
        const indexA = OrderedSet.indexOf(carbA.unit.elements, carbA.anomericCarbon)
        const indexB = OrderedSet.indexOf(carbB.unit.elements, carbB.anomericCarbon)
        if (indexA !== -1 && indexB !== -1) {
            return Link.Loci(structure, [
                Link.Location(
                    carbA.unit, indexA as StructureElement.UnitIndex,
                    carbB.unit, indexB as StructureElement.UnitIndex
                )
            ])
        }
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