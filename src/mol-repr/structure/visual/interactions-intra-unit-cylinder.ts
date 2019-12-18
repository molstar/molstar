/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../../mol-model/structure';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval } from '../../../mol-data/int';
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
import { Interactions } from '../../../mol-model-props/computed/interactions/interactions';

async function createIntraUnitInteractionsCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<InteractionsIntraUnitParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const interactions = InteractionsProvider.getValue(structure).value!
    const features = interactions.unitsFeatures.get(unit.id)
    const links = interactions.unitsLinks.get(unit.id)

    const { x, y, z } = features
    const { edgeCount, a, b } = links
    const { sizeFactor } = props

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
        ignore: () => false
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
        getLoci: getInteractionLoci,
        eachLocation: eachInteraction,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InteractionsIntraUnitParams>, currentProps: PD.Values<InteractionsIntraUnitParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            )
        }
    }, materialId)
}

function getInteractionLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = structure.unitMap.get(group.units[instanceId].id)
        const interactions = InteractionsProvider.getValue(structure).value!
        const links = interactions.unitsLinks.get(unit.id)
        return Interactions.Loci(structure, interactions, [
            { unitA: unit, indexA: links.a[groupId], unitB: unit, indexB: links.b[groupId] },
            { unitA: unit, indexA: links.b[groupId], unitB: unit, indexB: links.a[groupId] },
        ])
    }
    return EmptyLoci
}

function eachInteraction(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Interactions.isLoci(loci)) {
        const { structure, group } = structureGroup
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        const interactions = InteractionsProvider.getValue(structure).value!
        if (loci.interactions !== interactions) return false
        const unit = group.units[0]
        const links = interactions.unitsLinks.get(unit.id)
        const groupCount = links.edgeCount * 2
        for (const l of loci.links) {
            const unitIdx = group.unitIndexMap.get(l.unitA.id)
            if (unitIdx !== undefined) {
                const idx = links.getDirectedEdgeIndex(l.indexA, l.indexB)
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true
                }
            }
        }
    }
    return changed
}

function createInteractionsIterator(structureGroup: StructureGroup): LocationIterator {
    const { structure, group } = structureGroup
    const unit = group.units[0]
    const interactions = InteractionsProvider.getValue(structure).value!
    const links = interactions.unitsLinks.get(unit.id)
    const groupCount = links.edgeCount * 2
    const instanceCount = group.units.length
    const location = Interactions.Location(interactions)
    const getLocation = (groupIndex: number, instanceIndex: number) => {
        const instanceUnit = group.units[instanceIndex]
        location.unitA = instanceUnit
        location.indexA = links.a[groupIndex]
        location.unitB = instanceUnit
        location.indexB = links.b[groupIndex]
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation)
}