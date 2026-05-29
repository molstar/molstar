/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../../mol-repr/visual';
import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { createLinkCylinderMesh, LinkCylinderParams, LinkStyle } from '../../../mol-repr/structure/visual/util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../../../mol-repr/structure/complex-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { NullLocation } from '../../../mol-model/location';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { InteractionsProvider } from '../interactions';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { BridgeContacts, Bridges } from '../interactions/interactions';
import { Sphere3D } from '../../../mol-math/geometry';
import { InteractionsSharedParams } from './shared';
import { Features } from '../interactions/features';

type CanonicalLegIndices = {
    endpointA: Int32Array
    endpointB: Int32Array
};

const CanonicalLegIndicesCache = new WeakMap<BridgeContacts, CanonicalLegIndices>();

function getCanonicalLegIndices(bridges: BridgeContacts): CanonicalLegIndices {
    const cached = CanonicalLegIndicesCache.get(bridges);
    if (cached) return cached;

    const n = bridges.length;
    const endpointA = new Int32Array(n);
    const endpointB = new Int32Array(n);

    const legA = new Map<string, number>();
    const legB = new Map<string, number>();

    for (let i = 0; i < n; i++) {
        const b = bridges[i];

        const kA = `${b.unitA}|${b.indexA}|${b.unitM}|${b.indexMA}`;
        const kB = `${b.unitM}|${b.indexMB}|${b.unitB}|${b.indexB}`;

        let ai = legA.get(kA);
        if (ai === undefined) {
            ai = i;
            legA.set(kA, i);
        }
        endpointA[i] = ai;

        let bi = legB.get(kB);
        if (bi === undefined) {
            bi = i;
            legB.set(kB, i);
        }
        endpointB[i] = bi;
    }

    const indices = { endpointA, endpointB };
    CanonicalLegIndicesCache.set(bridges, indices);
    return indices;
}

function getFeatureMember(features: Features, featureIndex: Features.FeatureIndex): StructureElement.UnitIndex {
    return features.members[features.offsets[featureIndex]] as StructureElement.UnitIndex;
}

function atomPosition(unit: Unit.Atomic, features: Features, featureIndex: Features.FeatureIndex, out: Vec3) {
    const atomLocalIdx = getFeatureMember(features, featureIndex);
    unit.conformation.position(unit.elements[atomLocalIdx], out);
}

function setFeatureLocation(
    structure: Structure,
    location: StructureElement.Location,
    unitId: number,
    features: Features,
    featureIndex: Features.FeatureIndex
) {
    const unit = structure.unitMap.get(unitId) as Unit.Atomic;
    const atomLocalIdx = getFeatureMember(features, featureIndex);

    location.unit = unit;
    location.element = unit.elements[atomLocalIdx];
}

function applyLegA(
    bridgeIndex: number,
    bridgeCount: number,
    canonical: CanonicalLegIndices,
    apply: (interval: Interval) => boolean
) {
    let changed = false;
    const i = canonical.endpointA[bridgeIndex];

    if (apply(Interval.ofSingleton(i))) changed = true;
    if (apply(Interval.ofSingleton(i + bridgeCount))) changed = true;

    return changed;
}

function applyLegB(
    bridgeIndex: number,
    bridgeCount: number,
    canonical: CanonicalLegIndices,
    apply: (interval: Interval) => boolean
) {
    let changed = false;
    const i = canonical.endpointB[bridgeIndex];

    if (apply(Interval.ofSingleton(i + 2 * bridgeCount))) changed = true;
    if (apply(Interval.ofSingleton(i + 3 * bridgeCount))) changed = true;

    return changed;
}

function createBridgeCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<BridgeParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh);

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return Mesh.createEmpty(mesh);

    const { bridges, unitsFeatures } = interactions;

    const n = bridges.length;
    if (!n) return Mesh.createEmpty(mesh);

    const l = StructureElement.Location.create(structure);
    const { sizeFactor } = props;
    const canonical = getCanonicalLegIndices(bridges);

    const builderProps = {
        // Four half-cylinders per bridge; createLinkCylinderMesh draws the A-side half per call:
        //   [0,   n): A→mediator, forward  (A side)
        //   [n,  2n): A→mediator, backward (mediator side)
        //   [2n, 3n): mediator→B, forward  (mediator side)
        //   [3n, 4n): mediator→B, backward (B side)
        //
        // When multiple bridges share the same physical leg, only the first
        // occurrence is drawn; later ones map back to the canonical edge index.
        linkCount: 4 * n,

        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = bridges[edgeIndex % n];
            const uM = structure.unitMap.get(b.unitM) as Unit.Atomic;
            const fM = unitsFeatures.get(b.unitM);
            const leg = Math.floor(edgeIndex / n);

            if (leg === 0) {
                const uA = structure.unitMap.get(b.unitA) as Unit.Atomic;
                const fA = unitsFeatures.get(b.unitA);
                atomPosition(uA, fA, b.indexA, posA);
                atomPosition(uM, fM, b.indexMA, posB);
            } else if (leg === 1) {
                const uA = structure.unitMap.get(b.unitA) as Unit.Atomic;
                const fA = unitsFeatures.get(b.unitA);
                atomPosition(uM, fM, b.indexMA, posA);
                atomPosition(uA, fA, b.indexA, posB);
            } else if (leg === 2) {
                const uB = structure.unitMap.get(b.unitB) as Unit.Atomic;
                const fB = unitsFeatures.get(b.unitB);
                atomPosition(uM, fM, b.indexMB, posA);
                atomPosition(uB, fB, b.indexB, posB);
            } else {
                const uB = structure.unitMap.get(b.unitB) as Unit.Atomic;
                const fB = unitsFeatures.get(b.unitB);
                atomPosition(uB, fB, b.indexB, posA);
                atomPosition(uM, fM, b.indexMB, posB);
            }
        },

        ignore: (edgeIndex: number) => {
            const bi = edgeIndex % n;
            const leg = Math.floor(edgeIndex / n);

            return leg <= 1
                ? canonical.endpointA[bi] !== bi
                : canonical.endpointB[bi] !== bi;
        },

        style: (_edgeIndex: number) => LinkStyle.Dashed,

        radius: (edgeIndex: number) => {
            const b = bridges[edgeIndex % n];
            const leg = Math.floor(edgeIndex / n);
            const isLegA = leg <= 1;

            if (isLegA) {
                const fA = unitsFeatures.get(b.unitA);
                const fM = unitsFeatures.get(b.unitM);

                setFeatureLocation(structure, l, b.unitA, fA, b.indexA);
                const sizeA = theme.size.size(l);

                setFeatureLocation(structure, l, b.unitM, fM, b.indexMA);
                const sizeM = theme.size.size(l);

                return Math.min(sizeA, sizeM) * sizeFactor;
            } else {
                const fM = unitsFeatures.get(b.unitM);
                const fB = unitsFeatures.get(b.unitB);

                setFeatureLocation(structure, l, b.unitM, fM, b.indexMB);
                const sizeM = theme.size.size(l);

                setFeatureLocation(structure, l, b.unitB, fB, b.indexB);
                const sizeB = theme.size.size(l);

                return Math.min(sizeM, sizeB) * sizeFactor;
            }
        },
    };

    const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    if (boundingSphere) {
        m.setBoundingSphere(boundingSphere);
    } else if (m.triangleCount > 0) {
        const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, sizeFactor);
        m.setBoundingSphere(sphere);
    }

    return m;
}

export const BridgeParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    ...InteractionsSharedParams,
};
export type BridgeParams = typeof BridgeParams

export function BridgeVisual(materialId: number): ComplexVisual<BridgeParams> {
    return ComplexMeshVisual<BridgeParams>({
        defaultProps: PD.getDefaultValues(BridgeParams),
        createGeometry: createBridgeCylinderMesh,
        createLocationIterator: createBridgeIterator,
        getLoci: getBridgeLoci,
        eachLocation: eachBridgeInteraction,

        setUpdateState: (
            state: VisualUpdateState,
            newProps: PD.Values<BridgeParams>,
            currentProps: PD.Values<BridgeParams>,
            newTheme: Theme,
            currentTheme: Theme,
            newStructure: Structure,
            _currentStructure: Structure
        ) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.dashScale !== currentProps.dashScale ||
                newProps.dashCap !== currentProps.dashCap ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newTheme.size !== currentTheme.size
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

function getBridgeLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id !== objectId) return EmptyLoci;

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return EmptyLoci;

    const { bridges, unitsFeatures } = interactions;
    const n = bridges.length;

    if (!n || groupId < 0 || groupId >= 4 * n) return EmptyLoci;

    const bridgeIndex = groupId % n;

    return Bridges.Loci({ structure, bridges, unitsFeatures }, [{ bridgeIndex }]);
}

const __unitMap = new Map<number, OrderedSet<StructureElement.UnitIndex>>();

function eachBridgeInteraction(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, _isMarking: boolean) {
    let changed = false;

    if (Bridges.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;

        const { bridges } = interactions;
        const n = bridges.length;
        if (!n) return false;

        const canonical = getCanonicalLegIndices(bridges);

        for (const e of loci.elements) {
            if (e.bridgeIndex < 0 || e.bridgeIndex >= n) continue;

            if (applyLegA(e.bridgeIndex, n, canonical, apply)) changed = true;
            if (applyLegB(e.bridgeIndex, n, canonical, apply)) changed = true;
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;

        const { bridges, unitsFeatures } = interactions;
        const n = bridges.length;
        if (!n) return false;

        const canonical = getCanonicalLegIndices(bridges);

        __unitMap.clear();
        for (const e of loci.elements) {
            __unitMap.set(e.unit.id, e.indices);
        }

        for (let i = 0; i < n; i++) {
            const b = bridges[i];

            const indicesA = __unitMap.get(b.unitA);
            const indicesM = __unitMap.get(b.unitM);
            const indicesB = __unitMap.get(b.unitB);

            if (!indicesA && !indicesM && !indicesB) continue;

            let hitA = false;
            if (indicesA) {
                const fA = unitsFeatures.get(b.unitA);
                const mi = getFeatureMember(fA, b.indexA);
                hitA = OrderedSet.has(indicesA, mi);
            }

            let hitM = false;
            if (indicesM) {
                const fM = unitsFeatures.get(b.unitM);
                const miA = getFeatureMember(fM, b.indexMA);
                const miB = getFeatureMember(fM, b.indexMB);
                hitM = OrderedSet.has(indicesM, miA) || OrderedSet.has(indicesM, miB);
            }

            let hitB = false;
            if (indicesB) {
                const fB = unitsFeatures.get(b.unitB);
                const mi = getFeatureMember(fB, b.indexB);
                hitB = OrderedSet.has(indicesB, mi);
            }

            if (hitA || hitM) {
                if (applyLegA(i, n, canonical, apply)) changed = true;
            }

            if (hitB || hitM) {
                if (applyLegB(i, n, canonical, apply)) changed = true;
            }
        }

        __unitMap.clear();
    }

    return changed;
}

function createBridgeIterator(structure: Structure): LocationIterator {
    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return LocationIterator(0, 1, 1, () => NullLocation, true);

    const { bridges, unitsFeatures } = interactions;

    const n = bridges.length;
    const groupCount = 4 * n;
    const instanceCount = 1;

    const data: Bridges.Data = { structure, bridges, unitsFeatures };
    const location = Bridges.Location(data);
    const { element } = location;

    const getLocation = (groupIndex: number) => {
        element.bridgeIndex = n === 0 ? 0 : groupIndex % n;
        return location;
    };

    return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}
