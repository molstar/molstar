/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Link, StructureElement } from 'mol-model/structure';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { Vec3 } from 'mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderParams } from './util/link';
import { OrderedSet, Interval } from 'mol-data/int';
import { ComplexMeshVisual, ComplexVisual } from '../complex-visual';
import { LinkType } from 'mol-model/structure/model/types';
import { BitFlags } from 'mol-util';
import { UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { VisualUpdateState } from '../../util';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';

function createCarbohydrateLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateLinkParams>, mesh?: Mesh) {
    const { links, elements } = structure.carbohydrates
    const { linkSizeFactor } = props

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
            return theme.size.size(location) * linkSizeFactor
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const CarbohydrateLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    linkSizeFactor: PD.Numeric(0.3, { min: 0, max: 3, step: 0.01 }),
}
export type CarbohydrateLinkParams = typeof CarbohydrateLinkParams

export function CarbohydrateLinkVisual(): ComplexVisual<CarbohydrateLinkParams> {
    return ComplexMeshVisual<CarbohydrateLinkParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateLinkParams),
        createGeometry: createCarbohydrateLinkCylinderMesh,
        createLocationIterator: CarbohydrateLinkIterator,
        getLoci: getLinkLoci,
        eachLocation: eachCarbohydrateLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateLinkParams>, currentProps: PD.Values<CarbohydrateLinkParams>) => {
            state.createGeometry = (
                newProps.linkSizeFactor !== currentProps.linkSizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            )
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
                ),
                Link.Location(
                    carbB.unit, indexB as StructureElement.UnitIndex,
                    carbA.unit, indexA as StructureElement.UnitIndex
                )
            ])
        }
    }
    return EmptyLoci
}

function eachCarbohydrateLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Link.isLoci(loci)) {
        if (loci.structure !== structure) return false
        const { getLinkIndex } = structure.carbohydrates
        for (const l of loci.links) {
            const idx = getLinkIndex(l.aUnit, l.aUnit.elements[l.aIndex], l.bUnit, l.bUnit.elements[l.bIndex])
            if (idx !== undefined) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    } else if (StructureElement.isLoci(loci)) {
        if (loci.structure !== structure) return false
        // TODO mark link only when both of the link elements are in a StructureElement.Loci
        const { getElementIndex, getLinkIndices, elements } = structure.carbohydrates
        for (const e of loci.elements) {
            OrderedSet.forEach(e.indices, v => {
                const carbI = getElementIndex(e.unit, e.unit.elements[v])
                if (carbI !== undefined) {
                    const carb = elements[carbI]
                    const indices = getLinkIndices(carb.unit, carb.anomericCarbon)
                    for (let i = 0, il = indices.length; i < il; ++i) {
                        if (apply(Interval.ofSingleton(indices[i]))) changed = true
                    }
                }
            })
        }
    }
    return changed
}