/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval } from '../../../mol-data/int';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme } from '../../../mol-theme/theme';
import { InteractionsProvider } from '../interactions';
import { createLinkCylinderMesh, LinkCylinderParams, LinkCylinderStyle } from '../../../mol-repr/structure/visual/util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, StructureGroup } from '../../../mol-repr/structure/units-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Interactions } from '../interactions/interactions';
import { InteractionFlag } from '../interactions/common';

async function createIntraUnitInteractionsCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<InteractionsIntraUnitParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const location = StructureElement.Location.create(structure, unit);

    const interactions = InteractionsProvider.get(structure).value!;
    const features = interactions.unitsFeatures.get(unit.id);
    const contacts = interactions.unitsContacts.get(unit.id);

    const { x, y, z, members, offsets } = features;
    const { edgeCount, a, b, edgeProps: { flag } } = contacts;
    const { sizeFactor } = props;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const builderProps = {
        linkCount: edgeCount * 2,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            Vec3.set(posA, x[a[edgeIndex]], y[a[edgeIndex]], z[a[edgeIndex]]);
            Vec3.set(posB, x[b[edgeIndex]], y[b[edgeIndex]], z[b[edgeIndex]]);
        },
        style: (edgeIndex: number) => LinkCylinderStyle.Dashed,
        radius: (edgeIndex: number) => {
            location.element = unit.elements[members[offsets[a[edgeIndex]]]];
            const sizeA = theme.size.size(location);
            location.element = unit.elements[members[offsets[b[edgeIndex]]]];
            const sizeB = theme.size.size(location);
            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: (edgeIndex: number) => flag[edgeIndex] === InteractionFlag.Filtered
    };

    return createLinkCylinderMesh(ctx, builderProps, props, mesh);
}

export const InteractionsIntraUnitParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
};
export type InteractionsIntraUnitParams = typeof InteractionsIntraUnitParams

export function InteractionsIntraUnitVisual(materialId: number): UnitsVisual<InteractionsIntraUnitParams> {
    return UnitsMeshVisual<InteractionsIntraUnitParams>({
        defaultProps: PD.getDefaultValues(InteractionsIntraUnitParams),
        createGeometry: createIntraUnitInteractionsCylinderMesh,
        createLocationIterator: createInteractionsIterator,
        getLoci: getInteractionLoci,
        eachLocation: eachInteraction,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InteractionsIntraUnitParams>, currentProps: PD.Values<InteractionsIntraUnitParams>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            );

            const interactionsHash = InteractionsProvider.get(newStructureGroup.structure).version;
            if ((state.info.interactionsHash as number) !== interactionsHash) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.info.interactionsHash = interactionsHash;
            }
        }
    }, materialId);
}

function getInteractionLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = structure.unitMap.get(group.units[instanceId].id);
        const interactions = InteractionsProvider.get(structure).value!;
        const { a, b } = interactions.unitsContacts.get(unit.id);
        return Interactions.Loci(structure, interactions, [
            { unitA: unit, indexA: a[groupId], unitB: unit, indexB: b[groupId] },
            { unitA: unit, indexA: b[groupId], unitB: unit, indexB: a[groupId] },
        ]);
    }
    return EmptyLoci;
}

function eachInteraction(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Interactions.isLoci(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
        const interactions = InteractionsProvider.get(structure).value!;
        if (loci.data.interactions !== interactions) return false;
        const unit = group.units[0];
        const contacts = interactions.unitsContacts.get(unit.id);
        const groupCount = contacts.edgeCount * 2;
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unitA.id);
            if (unitIdx !== undefined) {
                const idx = contacts.getDirectedEdgeIndex(e.indexA, e.indexB);
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true;
                }
            }
        }
    }
    return changed;
}

function createInteractionsIterator(structureGroup: StructureGroup): LocationIterator {
    const { structure, group } = structureGroup;
    const unit = group.units[0];
    const interactions = InteractionsProvider.get(structure).value!;
    const contacts = interactions.unitsContacts.get(unit.id);
    const groupCount = contacts.edgeCount * 2;
    const instanceCount = group.units.length;
    const location = Interactions.Location(interactions, structure);
    const { element } = location;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
        const instanceUnit = group.units[instanceIndex];
        element.unitA = instanceUnit;
        element.indexA = contacts.a[groupIndex];
        element.unitB = instanceUnit;
        element.indexB = contacts.b[groupIndex];
        return location;
    };
    return LocationIterator(groupCount, instanceCount, getLocation);
}