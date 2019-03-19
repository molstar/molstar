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

function createCarbohydrateTerminalLinkCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateTerminalLinkParams>, mesh?: Mesh) {
    const { terminalLinks, elements } = structure.carbohydrates
    const { linkSizeFactor } = props

    const location = StructureElement.create()

    const builderProps = {
        linkCount: terminalLinks.length,
        referencePosition: (edgeIndex: number) => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const l = terminalLinks[edgeIndex]
            if (l.fromCarbohydrate) {
                Vec3.copy(posA, elements[l.carbohydrateIndex].geometry.center)
                l.elementUnit.conformation.position(l.elementUnit.elements[l.elementIndex], posB)
            } else {
                l.elementUnit.conformation.position(l.elementUnit.elements[l.elementIndex], posA)
                Vec3.copy(posB, elements[l.carbohydrateIndex].geometry.center)
            }
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => BitFlags.create(LinkType.Flag.None),
        radius: (edgeIndex: number) => {
            const l = terminalLinks[edgeIndex]
            if (l.fromCarbohydrate) {
                location.unit = elements[l.carbohydrateIndex].unit
                location.element = elements[l.carbohydrateIndex].anomericCarbon
            } else {
                location.unit = l.elementUnit
                location.element = l.elementUnit.elements[l.elementIndex]
            }
            return theme.size.size(location) * linkSizeFactor
        }
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const CarbohydrateTerminalLinkParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    linkSizeFactor: PD.Numeric(0.3, { min: 0, max: 3, step: 0.01 }),
}
export type CarbohydrateTerminalLinkParams = typeof CarbohydrateTerminalLinkParams

export function CarbohydrateTerminalLinkVisual(materialId: number): ComplexVisual<CarbohydrateTerminalLinkParams> {
    return ComplexMeshVisual<CarbohydrateTerminalLinkParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateTerminalLinkParams),
        createGeometry: createCarbohydrateTerminalLinkCylinderMesh,
        createLocationIterator: CarbohydrateTerminalLinkIterator,
        getLoci: getTerminalLinkLoci,
        eachLocation: eachTerminalLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateTerminalLinkParams>, currentProps: PD.Values<CarbohydrateTerminalLinkParams>) => {
            state.createGeometry = (
                newProps.linkSizeFactor !== currentProps.linkSizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            )
        }
    }, materialId)
}

function CarbohydrateTerminalLinkIterator(structure: Structure): LocationIterator {
    const { elements, terminalLinks } = structure.carbohydrates
    const groupCount = terminalLinks.length
    const instanceCount = 1
    const location = Link.Location()
    const getLocation = (groupIndex: number) => {
        const terminalLink = terminalLinks[groupIndex]
        const carb = elements[terminalLink.carbohydrateIndex]
        const indexCarb = OrderedSet.indexOf(carb.unit.elements, carb.anomericCarbon)
        if (terminalLink.fromCarbohydrate) {
            location.aUnit = carb.unit
            location.aIndex = indexCarb as StructureElement.UnitIndex
            location.bUnit = terminalLink.elementUnit
            location.bIndex = terminalLink.elementIndex
        } else {
            location.aUnit = terminalLink.elementUnit
            location.aIndex = terminalLink.elementIndex
            location.bUnit = carb.unit
            location.bIndex = indexCarb as StructureElement.UnitIndex
        }
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true)
}

function getTerminalLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const { terminalLinks, elements } = structure.carbohydrates
        const l = terminalLinks[groupId]
        const carb = elements[l.carbohydrateIndex]
        const carbIndex = OrderedSet.indexOf(carb.unit.elements, carb.anomericCarbon)

        return Link.Loci(structure, [
            Link.Location(
                carb.unit, carbIndex as StructureElement.UnitIndex,
                l.elementUnit, l.elementIndex
            ),
            Link.Location(
                l.elementUnit, l.elementIndex,
                carb.unit, carbIndex as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

// TODO for each link when both of the link elements are in a StructureElement.Loci
function eachTerminalLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getTerminalLinkIndex } = structure.carbohydrates
    let changed = false
    if (Link.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        for (const l of loci.links) {
            const idx = getTerminalLinkIndex(l.aUnit, l.aUnit.elements[l.aIndex], l.bUnit, l.bUnit.elements[l.bIndex])
            if (idx !== undefined) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    } else if (StructureElement.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        // TODO mark link only when both of the link elements are in a StructureElement.Loci
        const { getElementIndex, getTerminalLinkIndices, elements } = structure.carbohydrates
        for (const e of loci.elements) {
            OrderedSet.forEach(e.indices, v => {
                const carbI = getElementIndex(e.unit, e.unit.elements[v])
                if (carbI !== undefined) {
                    const carb = elements[carbI]
                    const indices = getTerminalLinkIndices(carb.unit, carb.anomericCarbon)
                    for (let i = 0, il = indices.length; i < il; ++i) {
                        if (apply(Interval.ofSingleton(indices[i]))) changed = true
                    }
                }
            })
        }
    }
    return changed
}