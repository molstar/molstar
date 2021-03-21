/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../../mol-repr/visual';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderParams, LinkStyle } from '../../../mol-repr/structure/visual/util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../../../mol-repr/structure/complex-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval, OrderedSet, SortedArray } from '../../../mol-data/int';
import { Interactions } from '../interactions/interactions';
import { InteractionsProvider } from '../interactions';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { InteractionFlag } from '../interactions/common';
import { Unit } from '../../../mol-model/structure/structure';
import { Sphere3D } from '../../../mol-math/geometry';

function createInterUnitInteractionCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InteractionsInterUnitParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh);

    const l = StructureElement.Location.create(structure);
    const interactions = InteractionsProvider.get(structure).value!;
    const { contacts, unitsFeatures } = interactions;

    const { edgeCount, edges } = contacts;
    const { sizeFactor } = props;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const { child } = structure;

    const builderProps = {
        linkCount: edgeCount,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const { unitA, indexA, unitB, indexB } = edges[edgeIndex];
            const fA = unitsFeatures.get(unitA);
            const fB = unitsFeatures.get(unitB);
            const uA = structure.unitMap.get(unitA);
            const uB = structure.unitMap.get(unitB);

            Vec3.set(posA, fA.x[indexA], fA.y[indexA], fA.z[indexA]);
            Vec3.transformMat4(posA, posA, uA.conformation.operator.matrix);

            Vec3.set(posB, fB.x[indexB], fB.y[indexB], fB.z[indexB]);
            Vec3.transformMat4(posB, posB, uB.conformation.operator.matrix);
        },
        style: (edgeIndex: number) => LinkStyle.Dashed,
        radius: (edgeIndex: number) => {
            const b = edges[edgeIndex];
            const fA = unitsFeatures.get(b.unitA);
            l.unit = structure.unitMap.get(b.unitA);
            l.element = l.unit.elements[fA.members[fA.offsets[b.indexA]]];
            const sizeA = theme.size.size(l);
            const fB = unitsFeatures.get(b.unitB);
            l.unit = structure.unitMap.get(b.unitB);
            l.element = l.unit.elements[fB.members[fB.offsets[b.indexB]]];
            const sizeB = theme.size.size(l);
            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: (edgeIndex: number) => {
            if (edges[edgeIndex].props.flag === InteractionFlag.Filtered) return true;

            if (child) {
                const b = edges[edgeIndex];
                const childUnitA = child.unitMap.get(b.unitA);
                if (!childUnitA) return true;

                const unitA = structure.unitMap.get(b.unitA);
                const fA = unitsFeatures.get(b.unitA);
                // TODO: check all members
                const eA = unitA.elements[fA.members[fA.offsets[b.indexA]]];
                if (!SortedArray.has(childUnitA.elements, eA)) return true;
            }

            return false;
        }
    };

    const m = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    const sphere = Sphere3D.expand(Sphere3D(), (child ?? structure).boundary.sphere, 1 * sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export const InteractionsInterUnitParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    dashCount: PD.Numeric(6, { min: 2, max: 10, step: 2 }),
    dashScale: PD.Numeric(0.4, { min: 0, max: 2, step: 0.1 }),
    includeParent: PD.Boolean(false),
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
                newProps.dashCount !== currentProps.dashCount ||
                newProps.dashScale !== currentProps.dashScale ||
                newProps.dashCap !== currentProps.dashCap ||
                newProps.radialSegments !== currentProps.radialSegments
            );

            const interactionsHash = InteractionsProvider.get(newStructure).version;
            if ((state.info.interactionsHash as number) !== interactionsHash) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
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
        const unitA = structure.unitMap.get(c.unitA);
        const unitB = structure.unitMap.get(c.unitB);
        return Interactions.Loci(structure, interactions, [
            { unitA: unitA, indexA: c.indexA, unitB: unitB, indexB: c.indexB },
            { unitA: unitB, indexA: c.indexB, unitB: unitA, indexB: c.indexA },
        ]);
    }
    return EmptyLoci;
}

function eachInteraction(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, isMarking: boolean) {
    let changed = false;
    if (Interactions.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
        const interactions = InteractionsProvider.get(structure).value!;
        if (loci.data.interactions !== interactions) return false;
        const { contacts } = interactions;

        for (const c of loci.elements) {
            const idx = contacts.getEdgeIndex(c.indexA, c.unitA.id, c.indexB, c.unitB.id);
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true;
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        if (isMarking && loci.elements.length === 1) return false; // only a single unit

        const contacts = InteractionsProvider.get(structure).value?.contacts;
        if (!contacts) return false;

        // TODO when isMarking, all elements of contact features need to be in the loci
        for (const e of loci.elements) {
            const { unit } = e;
            if (!Unit.isAtomic(unit)) continue;
            if (isMarking && OrderedSet.size(e.indices) === 1) continue;

            OrderedSet.forEach(e.indices, v => {
                for (const idx of contacts.getContactIndicesForElement(v, unit)) {
                    if (apply(Interval.ofSingleton(idx))) changed = true;
                }
            });
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
        element.unitA = structure.unitMap.get(c.unitA);
        element.indexA = c.indexA;
        element.unitB = structure.unitMap.get(c.unitB);
        element.indexB = c.indexB;
        return location;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}