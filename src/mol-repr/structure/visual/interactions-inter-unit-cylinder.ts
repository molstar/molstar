/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { createBondCylinderMesh, BondCylinderParams } from './util/bond';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval } from '../../../mol-data/int';
import { BondType } from '../../../mol-model/structure/model/types';
import { Interactions } from '../../../mol-model-props/computed/interactions/interactions';
import { InteractionsProvider } from '../../../mol-model-props/computed/interactions';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';

function createInterUnitInteractionCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InteractionsInterUnitParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh)

    const interactions = InteractionsProvider.getValue(structure).value!
    const { links, unitsFeatures } = interactions

    const { edgeCount, edges } = links
    const { sizeFactor } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const builderProps = {
        bondCount: edgeCount,
        referencePosition: () => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const { unitA, indexA, unitB, indexB } = edges[edgeIndex]
            const fA = unitsFeatures.get(unitA.id)
            const fB = unitsFeatures.get(unitB.id)
            Vec3.set(posA, fA.x[indexA], fA.y[indexA], fA.z[indexA])
            Vec3.transformMat4(posA, posA, unitA.conformation.operator.matrix)
            Vec3.set(posB, fB.x[indexB], fB.y[indexB], fB.z[indexB])
            Vec3.transformMat4(posB, posB, unitB.conformation.operator.matrix)
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => BondType.Flag.MetallicCoordination, // TODO
        radius: (edgeIndex: number) => sizeFactor,
        ignore: () => false
    }

    return createBondCylinderMesh(ctx, builderProps, props, mesh)
}

export const InteractionsInterUnitParams = {
    ...ComplexMeshParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
}
export type InteractionsInterUnitParams = typeof InteractionsInterUnitParams

export function InteractionsInterUnitVisual(materialId: number): ComplexVisual<InteractionsInterUnitParams> {
    return ComplexMeshVisual<InteractionsInterUnitParams>({
        defaultProps: PD.getDefaultValues(InteractionsInterUnitParams),
        createGeometry: createInterUnitInteractionCylinderMesh,
        createLocationIterator: createInteractionsIterator,
        getLoci: getInteractionLoci,
        eachLocation: eachInteraction,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InteractionsInterUnitParams>, currentProps: PD.Values<InteractionsInterUnitParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            )

            const interactionsHash = InteractionsProvider.getValue(newStructure).version
            if ((state.info.interactionsHash as number) !== interactionsHash) {
                state.createGeometry = true
                state.updateTransform = true
                state.info.interactionsHash = interactionsHash
            }
        }
    }, materialId)
}

function getInteractionLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        const interactions = InteractionsProvider.getValue(structure).value!
        const l = interactions.links.edges[groupId]
        return Interactions.Loci(structure, interactions, [
            { unitA: l.unitA, indexA: l.indexA, unitB: l.unitB, indexB: l.indexB },
            { unitA: l.unitB, indexA: l.indexB, unitB: l.unitA, indexB: l.indexA },
        ])
    }
    return EmptyLoci
}

function eachInteraction(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Interactions.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false
        const interactions = InteractionsProvider.getValue(structure).value!
        if (loci.interactions !== interactions) return false
        const { links } = interactions

        for (const l of loci.links) {
            const idx = links.getEdgeIndex(l.indexA, l.unitA, l.indexB, l.unitB)
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true
            }
        }
    }
    return changed
}

function createInteractionsIterator(structure: Structure): LocationIterator {
    const interactions = InteractionsProvider.getValue(structure).value!
    const { links } = interactions
    const groupCount = links.edgeCount
    const instanceCount = 1
    const location = Interactions.Location(interactions)
    const getLocation = (groupIndex: number) => {
        const link = links.edges[groupIndex]
        location.unitA = link.unitA
        location.indexA = link.indexA
        location.unitB = link.unitB
        location.indexB = link.indexB
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true)
}