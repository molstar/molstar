/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../../mol-repr/visual';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderParams, LinkCylinderStyle } from '../../../mol-repr/structure/visual/util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../../../mol-repr/structure/complex-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval } from '../../../mol-data/int';
import { Interactions } from '../interactions/interactions';
import { InteractionsProvider } from '../interactions';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { InteractionFlag } from '../interactions/common';

const tmpLoc = StructureElement.Location.create(void 0);

function createInterUnitInteractionCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InteractionsInterUnitParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh);

    const interactions = InteractionsProvider.get(structure).value!;
    const { contacts, unitsFeatures } = interactions;

    const { edgeCount, edges } = contacts;
    const { sizeFactor } = props;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const builderProps = {
        linkCount: edgeCount,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const { unitA, indexA, unitB, indexB } = edges[edgeIndex];
            const fA = unitsFeatures.get(unitA.id);
            const fB = unitsFeatures.get(unitB.id);

            Vec3.set(posA, fA.x[indexA], fA.y[indexA], fA.z[indexA]);
            Vec3.transformMat4(posA, posA, unitA.conformation.operator.matrix);

            Vec3.set(posB, fB.x[indexB], fB.y[indexB], fB.z[indexB]);
            Vec3.transformMat4(posB, posB, unitB.conformation.operator.matrix);
        },
        style: (edgeIndex: number) => LinkCylinderStyle.Dashed,
        radius: (edgeIndex: number) => {
            const b = edges[edgeIndex];
            const fA = unitsFeatures.get(b.unitA.id);
            tmpLoc.structure = structure;
            tmpLoc.unit = b.unitA;
            tmpLoc.element = b.unitA.elements[fA.members[fA.offsets[b.indexA]]];
            const sizeA = theme.size.size(tmpLoc);
            const fB = unitsFeatures.get(b.unitB.id);
            tmpLoc.unit = b.unitB;
            tmpLoc.element = b.unitB.elements[fB.members[fB.offsets[b.indexB]]];
            const sizeB = theme.size.size(tmpLoc);
            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: (edgeIndex: number) => edges[edgeIndex].props.flag === InteractionFlag.Filtered
    };

    return createLinkCylinderMesh(ctx, builderProps, props, mesh);
}

export const InteractionsInterUnitParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
};
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
            );

            const interactionsHash = InteractionsProvider.get(newStructure).version;
            if ((state.info.interactionsHash as number) !== interactionsHash) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.info.interactionsHash = interactionsHash;
            }
        }
    }, materialId);
}

function getInteractionLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const interactions = InteractionsProvider.get(structure).value!;
        const c = interactions.contacts.edges[groupId];
        return Interactions.Loci(structure, interactions, [
            { unitA: c.unitA, indexA: c.indexA, unitB: c.unitB, indexB: c.indexB },
            { unitA: c.unitB, indexA: c.indexB, unitB: c.unitA, indexB: c.indexA },
        ]);
    }
    return EmptyLoci;
}

function eachInteraction(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Interactions.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
        const interactions = InteractionsProvider.get(structure).value!;
        if (loci.data.interactions !== interactions) return false;
        const { contacts } = interactions;

        for (const c of loci.elements) {
            const idx = contacts.getEdgeIndex(c.indexA, c.unitA, c.indexB, c.unitB);
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true;
            }
        }
    }
    return changed;
}

function createInteractionsIterator(structure: Structure): LocationIterator {
    const interactions = InteractionsProvider.get(structure).value!;
    const { contacts } = interactions;
    const groupCount = contacts.edgeCount;
    const instanceCount = 1;
    const location = Interactions.Location(interactions, structure);
    const { element } = location;
    const getLocation = (groupIndex: number) => {
        const c = contacts.edges[groupIndex];
        element.unitA = c.unitA;
        element.indexA = c.indexA;
        element.unitB = c.unitB;
        element.indexB = c.indexB;
        return location;
    };
    return LocationIterator(groupCount, instanceCount, getLocation, true);
}