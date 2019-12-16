/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Bond, Structure, StructureElement } from '../../../mol-model/structure';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval, SortedArray, OrderedSet } from '../../../mol-data/int';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { VisualContext } from '../../visual';
import { Theme } from '../../../mol-theme/theme';
import { BondType } from '../../../mol-model/structure/model/types';
import { InteractionsProvider } from '../../../mol-model-props/computed/interactions';
import { createBondCylinderMesh, BondCylinderParams } from './util/bond';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, StructureGroup } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';

async function createIntraUnitInteractionsCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<InteractionsIntraUnitParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const interactions = InteractionsProvider.getValue(structure).value!
    const { features, links } = interactions.get(unit.id)!

    const { x, y, z, offsets, members } = features
    const { edgeCount, a, b } = links
    const { sizeFactor } = props
    const { elements } = unit
    const rootElements = structure.root.unitMap.get(unit.id).elements

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const builderProps = {
        bondCount: edgeCount * 2,
        referencePosition: () => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            Vec3.set(posA, x[a[edgeIndex]], y[a[edgeIndex]], z[a[edgeIndex]])
            Vec3.set(posB, x[b[edgeIndex]], y[b[edgeIndex]], z[b[edgeIndex]])
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => BondType.Flag.MetallicCoordination, // TODO
        radius: (edgeIndex: number) => sizeFactor,
        ignore: elements !== rootElements ? (edgeIndex: number) => {
            for (let i = offsets[a[edgeIndex]], il = offsets[a[edgeIndex] + 1]; i < il; ++i) {
                if (!SortedArray.has(elements, rootElements[members[i]])) return true
            }
            for (let i = offsets[b[edgeIndex]], il = offsets[b[edgeIndex] + 1]; i < il; ++i) {
                if (!SortedArray.has(elements, rootElements[members[i]])) return true
            }
            return false
        } : () => false
    }

    return createBondCylinderMesh(ctx, builderProps, props, mesh)
}

export const InteractionsIntraUnitParams = {
    ...UnitsMeshParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
}
export type InteractionsIntraUnitParams = typeof InteractionsIntraUnitParams

export function InteractionsIntraUnitVisual(materialId: number): UnitsVisual<InteractionsIntraUnitParams> {
    return UnitsMeshVisual<InteractionsIntraUnitParams>({
        defaultProps: PD.getDefaultValues(InteractionsIntraUnitParams),
        createGeometry: createIntraUnitInteractionsCylinderMesh,
        createLocationIterator: createInteractionsIterator,
        getLoci: getLinkLoci,
        eachLocation: eachInteraction,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InteractionsIntraUnitParams>, currentProps: PD.Values<InteractionsIntraUnitParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.bondScale !== currentProps.bondScale ||
                newProps.bondSpacing !== currentProps.bondSpacing
            )
        }
    }, materialId)
}

function getLinkLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        // use corresponding unit from root structure
        const unit = structure.root.unitMap.get(group.units[instanceId].id)
        if (Unit.isAtomic(unit)) {
            const interactions = InteractionsProvider.getValue(structure).value!
            const { features, links } = interactions.get(unit.id)!
            const { members, offsets } = features
            // TODO this uses the first member elements of the features of an interaction as a representative
            return Bond.Loci(structure.root, [
                Bond.Location(
                    unit, members[offsets[links.a[groupId]]],
                    unit, members[offsets[links.b[groupId]]]
                ),
                Bond.Location(
                    unit, members[offsets[links.b[groupId]]],
                    unit, members[offsets[links.a[groupId]]]
                )
            ])
        }
    }
    return EmptyLoci
}

function eachInteraction(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Bond.isLoci(loci)) {
        const { structure, group } = structureGroup
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        const unit = group.units[0]
        if (!Unit.isAtomic(unit)) return false
        const interactions = InteractionsProvider.getValue(structure).value!
        const { links, getLinkIndex } = interactions.get(unit.id)!
        const groupCount = links.edgeCount * 2
        for (const b of loci.bonds) {
            const unitIdx = group.unitIndexMap.get(b.aUnit.id)
            if (unitIdx !== undefined) {
                const idx = getLinkIndex(b.aIndex, b.bIndex)
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        const { structure, group } = structureGroup
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        const unit = group.units[0]
        if (!Unit.isAtomic(unit)) return false
        const interactions = InteractionsProvider.getValue(structure).value!
        const { links, elementsIndex } = interactions.get(unit.id)!
        const groupCount = links.edgeCount * 2
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id)
            if (unitIdx !== undefined) {
                const { offset } = links
                const { indices, offsets } = elementsIndex
                const rootElements = structure.root.unitMap.get(unit.id).elements
                OrderedSet.forEach(e.indices, _v => {
                    const v = SortedArray.indexOf(rootElements, e.unit.elements[_v])
                    for (let i = offsets[v], il = offsets[v + 1]; i < il; ++i) {
                        const f = indices[i]
                        for (let t = offset[f], _t = offset[f + 1]; t < _t; t++) {
                            if (apply(Interval.ofSingleton(unitIdx * groupCount + t))) changed = true
                        }
                    }
                })
            }
        }
    }
    return changed
}

function createInteractionsIterator(structureGroup: StructureGroup): LocationIterator {
    const { structure, group } = structureGroup
    const unit = group.units[0]
    const interactions = InteractionsProvider.getValue(structure).value!
    const { links, features } = interactions.get(unit.id)!
    const { members, offsets } = features
    const groupCount = links.edgeCount * 2
    const instanceCount = group.units.length
    const location = Bond.Location()
    const getLocation = (groupIndex: number, instanceIndex: number) => {
        const fA = links.a[groupIndex]
        const fB = links.b[groupIndex]
        const instanceUnit = group.units[instanceIndex]
        location.aUnit = instanceUnit
        location.aIndex = members[offsets[fA]]
        location.bUnit = instanceUnit
        location.bIndex = members[offsets[fB]]
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation)
}