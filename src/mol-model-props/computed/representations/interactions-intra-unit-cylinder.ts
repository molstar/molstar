/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval, OrderedSet, SortedArray } from '../../../mol-data/int';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme } from '../../../mol-theme/theme';
import { InteractionsProvider } from '../interactions';
import { createLinkCylinderMesh, LinkCylinderParams, LinkStyle } from '../../../mol-repr/structure/visual/util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../../../mol-repr/structure/units-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Interactions } from '../interactions/interactions';
import { FeatureType, InteractionFlag } from '../interactions/common';
import { Sphere3D } from '../../../mol-math/geometry';
import { StructureGroup, isHydrogen } from '../../../mol-repr/structure/visual/util/common';
import { assertUnreachable } from '../../../mol-util/type-helpers';
import { InteractionsSharedParams } from './shared';
import { InteractionType } from '../interactions/common';
import { eachIntraBondedAtom } from '../chemistry/util';

async function createIntraUnitInteractionsCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<InteractionsIntraUnitParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Mesh.createEmpty(mesh);

    const location = StructureElement.Location.create(structure, unit);

    const interactions = InteractionsProvider.get(structure).value!;
    const features = interactions.unitsFeatures.get(unit.id);
    const contacts = interactions.unitsContacts.get(unit.id);

    const { x, y, z, members, offsets, types } = features;
    const { edgeCount, a, b, edgeProps: { flag, type } } = contacts;
    const { sizeFactor, ignoreHydrogens, ignoreHydrogensVariant, parentDisplay } = props;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const { elements, conformation: c } = unit;
    const p = Vec3();
    const pA = Vec3();
    const pB = Vec3();

    const builderProps = {
        linkCount: edgeCount * 2,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const t = type[edgeIndex];
            if ((!ignoreHydrogens || ignoreHydrogensVariant !== 'all') && (
                t === InteractionType.HydrogenBond || (t === InteractionType.WeakHydrogenBond && ignoreHydrogensVariant !== 'non-polar'))
            ) {
                const idxA = members[offsets[a[edgeIndex]]];
                const idxB = members[offsets[b[edgeIndex]]];
                c.invariantPosition(elements[idxA], pA);
                c.invariantPosition(elements[idxB], pB);
                let minDistA = Vec3.distance(pA, pB);
                let minDistB = minDistA;
                Vec3.copy(posA, pA);
                Vec3.copy(posB, pB);
                const donorType = t === InteractionType.HydrogenBond ? FeatureType.HydrogenDonor : FeatureType.WeakHydrogenDonor;
                const isHydrogenDonorA = types[offsets[a[edgeIndex]]] === donorType;

                if (isHydrogenDonorA) {
                    eachIntraBondedAtom(unit, idxA, (_, idx) => {
                        if (isHydrogen(structure, unit, elements[idx], 'all')) {
                            c.invariantPosition(elements[idx], p);
                            const dist = Vec3.distance(p, pB);
                            if (dist < minDistA) {
                                minDistA = dist;
                                Vec3.copy(posA, p);
                            }
                        }
                    });
                } else {
                    eachIntraBondedAtom(unit, idxB, (_, idx) => {
                        if (isHydrogen(structure, unit, elements[idx], 'all')) {
                            c.invariantPosition(elements[idx], p);
                            const dist = Vec3.distance(p, pA);
                            if (dist < minDistB) {
                                minDistB = dist;
                                Vec3.copy(posB, p);
                            }
                        }
                    });
                }
            } else {
                Vec3.set(posA, x[a[edgeIndex]], y[a[edgeIndex]], z[a[edgeIndex]]);
                Vec3.set(posB, x[b[edgeIndex]], y[b[edgeIndex]], z[b[edgeIndex]]);
            }
        },
        style: (edgeIndex: number) => LinkStyle.Dashed,
        radius: (edgeIndex: number) => {
            location.element = elements[members[offsets[a[edgeIndex]]]];
            const sizeA = theme.size.size(location);
            location.element = elements[members[offsets[b[edgeIndex]]]];
            const sizeB = theme.size.size(location);
            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: (edgeIndex: number) => {
            if (flag[edgeIndex] === InteractionFlag.Filtered) return true;

            if (childUnit) {
                if (parentDisplay === 'stub') {
                    const f = a[edgeIndex];
                    for (let i = offsets[f], il = offsets[f + 1]; i < il; ++i) {
                        const e = elements[members[offsets[i]]];
                        if (!SortedArray.has(childUnit.elements, e)) return true;
                    }
                } else if (parentDisplay === 'full' || parentDisplay === 'between') {
                    let flagA = false;
                    let flagB = false;

                    const fA = a[edgeIndex];
                    for (let i = offsets[fA], il = offsets[fA + 1]; i < il; ++i) {
                        const eA = elements[members[offsets[i]]];
                        if (!SortedArray.has(childUnit.elements, eA)) flagA = true;
                    }

                    const fB = b[edgeIndex];
                    for (let i = offsets[fB], il = offsets[fB + 1]; i < il; ++i) {
                        const eB = elements[members[offsets[i]]];
                        if (!SortedArray.has(childUnit.elements, eB)) flagB = true;
                    }

                    return parentDisplay === 'full' ? flagA && flagB : flagA === flagB;
                } else {
                    assertUnreachable(parentDisplay);
                }
            }

            return false;
        }
    };

    const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    if (boundingSphere) {
        m.setBoundingSphere(boundingSphere);
    } else if (m.triangleCount > 0) {
        const sphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, 1 * sizeFactor);
        m.setBoundingSphere(sphere);
    }

    return m;
}

export const InteractionsIntraUnitParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    ...InteractionsSharedParams,
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
                newProps.dashCount !== currentProps.dashCount ||
                newProps.dashScale !== currentProps.dashScale ||
                newProps.dashCap !== currentProps.dashCap ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.parentDisplay !== currentProps.parentDisplay
            );

            const interactionsHash = InteractionsProvider.get(newStructureGroup.structure).version;
            if ((state.info.interactionsHash as number) !== interactionsHash) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
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

const __contactIndicesSet = new Set<number>();

function eachInteraction(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean, isMarking: boolean) {
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
            if (e.unitA !== e.unitB) continue;
            const unitIdx = group.unitIndexMap.get(e.unitA.id);
            if (unitIdx !== undefined) {
                const idx = contacts.getDirectedEdgeIndex(e.indexA, e.indexB);
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true;
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areEquivalent(loci.structure, structure)) return false;

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;
        const unit = group.units[0];
        const contacts = interactions.unitsContacts.get(unit.id);
        const features = interactions.unitsFeatures.get(unit.id);
        const groupCount = contacts.edgeCount * 2;

        const { offset } = contacts;
        const { offsets: fOffsets, indices: fIndices } = features.elementsIndex;
        const { members, offsets } = features;

        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id);
            if (unitIdx === undefined) continue;

            OrderedSet.forEach(e.indices, v => {
                for (let i = fOffsets[v], il = fOffsets[v + 1]; i < il; ++i) {
                    const fI = fIndices[i];
                    for (let j = offset[fI], jl = offset[fI + 1]; j < jl; ++j) {
                        __contactIndicesSet.add(j);
                    }
                }
            });

            __contactIndicesSet.forEach(i => {
                if (isMarking) {
                    const fA = contacts.a[i];
                    for (let j = offsets[fA], jl = offsets[fA + 1]; j < jl; ++j) {
                        if (!OrderedSet.has(e.indices, members[j])) return;
                    }
                    const fB = contacts.b[i];
                    for (let j = offsets[fB], jl = offsets[fB + 1]; j < jl; ++j) {
                        if (!OrderedSet.has(e.indices, members[j])) return;
                    }
                }

                if (apply(Interval.ofSingleton(unitIdx * groupCount + i))) changed = true;
            });

            __contactIndicesSet.clear();
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
    return LocationIterator(groupCount, instanceCount, 1, getLocation);
}