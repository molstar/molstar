/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paul.pillot@tandemai.com>
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
import { FeatureType, InteractionFlag, InteractionType } from '../interactions/common';
import { Unit } from '../../../mol-model/structure/structure';
import { Sphere3D } from '../../../mol-math/geometry';
import { assertUnreachable } from '../../../mol-util/type-helpers';
import { InteractionsSharedParams } from './shared';
import { eachBondedAtom } from '../chemistry/util';
import { isHydrogen } from '../../../mol-repr/structure/visual/util/common';

function createInterUnitInteractionCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InteractionsInterUnitParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh);

    const l = StructureElement.Location.create(structure);
    const interactions = InteractionsProvider.get(structure).value!;
    const { contacts, unitsFeatures } = interactions;

    const { edgeCount, edges } = contacts;
    const { sizeFactor, ignoreHydrogens, ignoreHydrogensVariant, parentDisplay } = props;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const { child } = structure;
    const p = Vec3();
    const pA = Vec3();
    const pB = Vec3();

    const builderProps = {
        linkCount: edgeCount,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const { unitA, indexA, unitB, indexB, props: { type: t } } = edges[edgeIndex];
            const fA = unitsFeatures.get(unitA);
            const fB = unitsFeatures.get(unitB);
            const uA = structure.unitMap.get(unitA) as Unit.Atomic;
            const uB = structure.unitMap.get(unitB) as Unit.Atomic;

            if ((!ignoreHydrogens || ignoreHydrogensVariant !== 'all') && (
                t === InteractionType.HydrogenBond || (t === InteractionType.WeakHydrogenBond && ignoreHydrogensVariant !== 'non-polar'))
            ) {
                const idxA = fA.members[fA.offsets[indexA]];
                const idxB = fB.members[fB.offsets[indexB]];
                uA.conformation.position(uA.elements[idxA], pA);
                uB.conformation.position(uB.elements[idxB], pB);
                let minDistA = Vec3.distance(pA, pB);
                let minDistB = minDistA;
                Vec3.copy(posA, pA);
                Vec3.copy(posB, pB);
                const donorType = t === InteractionType.HydrogenBond ? FeatureType.HydrogenDonor : FeatureType.WeakHydrogenDonor;
                const isHydrogenDonorA = fA.types[fA.offsets[indexA]] === donorType;

                if (isHydrogenDonorA) {
                    eachBondedAtom(structure, uA, idxA, (u, idx) => {
                        const eI = u.elements[idx];
                        if (isHydrogen(structure, u, eI, 'all')) {
                            u.conformation.position(eI, p);
                            const dist = Vec3.distance(p, pB);
                            if (dist < minDistA) {
                                minDistA = dist;
                                Vec3.copy(posA, p);
                            }
                        }
                    });
                } else {
                    eachBondedAtom(structure, uB, idxB, (u, idx) => {
                        const eI = u.elements[idx];
                        if (isHydrogen(structure, u, eI, 'all')) {
                            u.conformation.position(eI, p);
                            const dist = Vec3.distance(p, pA);
                            if (dist < minDistB) {
                                minDistB = dist;
                                Vec3.copy(posB, p);
                            }
                        }
                    });
                }
            } else {
                Vec3.set(posA, fA.x[indexA], fA.y[indexA], fA.z[indexA]);
                Vec3.transformMat4(posA, posA, uA.conformation.operator.matrix);

                Vec3.set(posB, fB.x[indexB], fB.y[indexB], fB.z[indexB]);
                Vec3.transformMat4(posB, posB, uB.conformation.operator.matrix);
            }
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

                if (parentDisplay === 'stub') {
                    const childUnitA = child.unitMap.get(b.unitA);
                    if (!childUnitA) return true;

                    const unitA = structure.unitMap.get(b.unitA);
                    const { offsets, members } = unitsFeatures.get(b.unitA);
                    for (let i = offsets[b.indexA], il = offsets[b.indexA + 1]; i < il; ++i) {
                        const eA = unitA.elements[members[i]];
                        if (!SortedArray.has(childUnitA.elements, eA)) return true;
                    }
                } else if (parentDisplay === 'full' || parentDisplay === 'between') {
                    let flagA = false;
                    let flagB = false;

                    const childUnitA = child.unitMap.get(b.unitA);
                    if (!childUnitA) {
                        flagA = true;
                    } else {
                        const unitA = structure.unitMap.get(b.unitA);
                        const { offsets, members } = unitsFeatures.get(b.unitA);
                        for (let i = offsets[b.indexA], il = offsets[b.indexA + 1]; i < il; ++i) {
                            const eA = unitA.elements[members[i]];
                            if (!SortedArray.has(childUnitA.elements, eA)) flagA = true;
                        }
                    }

                    const childUnitB = child.unitMap.get(b.unitB);
                    if (!childUnitB) {
                        flagB = true;
                    } else {
                        const unitB = structure.unitMap.get(b.unitB);
                        const { offsets, members } = unitsFeatures.get(b.unitB);
                        for (let i = offsets[b.indexB], il = offsets[b.indexB + 1]; i < il; ++i) {
                            const eB = unitB.elements[members[i]];
                            if (!SortedArray.has(childUnitB.elements, eB)) flagB = true;
                        }
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
        const { child } = structure;
        const sphere = Sphere3D.expand(Sphere3D(), (child ?? structure).boundary.sphere, 1 * sizeFactor);
        m.setBoundingSphere(sphere);
    }

    return m;
}

export const InteractionsInterUnitParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    ...InteractionsSharedParams,
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
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.parentDisplay !== currentProps.parentDisplay
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

const __unitMap = new Map<number, OrderedSet<StructureElement.UnitIndex>>();
const __contactIndicesSet = new Set<number>();

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

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;

        const { contacts, unitsFeatures } = interactions;

        for (const e of loci.elements) __unitMap.set(e.unit.id, e.indices);

        for (const e of loci.elements) {
            const { unit } = e;
            if (!Unit.isAtomic(unit)) continue;

            OrderedSet.forEach(e.indices, v => {
                for (const idx of contacts.getContactIndicesForElement(v, unit)) {
                    __contactIndicesSet.add(idx);
                }
            });
        }

        __contactIndicesSet.forEach(i => {
            if (isMarking) {
                const { indexA, unitA, indexB, unitB } = contacts.edges[i];

                const indicesA = __unitMap.get(unitA);
                const indicesB = __unitMap.get(unitB);
                if (!indicesA || !indicesB) return;

                const { offsets: offsetsA, members: membersA } = unitsFeatures.get(unitA);
                for (let j = offsetsA[indexA], jl = offsetsA[indexA + 1]; j < jl; ++j) {
                    if (!OrderedSet.has(indicesA, membersA[j])) return;
                }

                const { offsets: offsetsB, members: membersB } = unitsFeatures.get(unitB);
                for (let j = offsetsB[indexB], jl = offsetsB[indexB + 1]; j < jl; ++j) {
                    if (!OrderedSet.has(indicesB, membersB[j])) return;
                }
            }

            if (apply(Interval.ofSingleton(i))) changed = true;
        });

        __unitMap.clear();
        __contactIndicesSet.clear();
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