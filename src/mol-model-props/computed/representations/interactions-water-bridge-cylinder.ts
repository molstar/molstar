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
import { WaterBridges } from '../interactions/water-bridges';
import { Sphere3D } from '../../../mol-math/geometry';
import { InteractionsSharedParams } from './shared';
import { Features } from '../interactions/features';

type WaterBridgeContacts = WaterBridges.Data['waterBridges'];

type CanonicalLegIndices = {
    donor: Int32Array
    acceptor: Int32Array
};

const CanonicalLegIndicesCache = new WeakMap<WaterBridgeContacts, CanonicalLegIndices>();

function getCanonicalLegIndices(waterBridges: WaterBridgeContacts): CanonicalLegIndices {
    const cached = CanonicalLegIndicesCache.get(waterBridges);
    if (cached) return cached;

    const n = waterBridges.length;
    const donor = new Int32Array(n);
    const acceptor = new Int32Array(n);

    const donorLegs = new Map<string, number>();
    const acceptorLegs = new Map<string, number>();

    for (let i = 0; i < n; i++) {
        const wb = waterBridges[i];

        const dk = `${wb.unitA}|${wb.indexA}|${wb.unitW}|${wb.indexWA}`;
        const ak = `${wb.unitW}|${wb.indexWD}|${wb.unitB}|${wb.indexB}`;

        let di = donorLegs.get(dk);
        if (di === undefined) {
            di = i;
            donorLegs.set(dk, i);
        }
        donor[i] = di;

        let ai = acceptorLegs.get(ak);
        if (ai === undefined) {
            ai = i;
            acceptorLegs.set(ak, i);
        }
        acceptor[i] = ai;
    }

    const indices = { donor, acceptor };
    CanonicalLegIndicesCache.set(waterBridges, indices);
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

function applyDonorLeg(
    bridgeIndex: number,
    bridgeCount: number,
    canonical: CanonicalLegIndices,
    apply: (interval: Interval) => boolean
) {
    let changed = false;
    const i = canonical.donor[bridgeIndex];

    if (apply(Interval.ofSingleton(i))) changed = true;
    if (apply(Interval.ofSingleton(i + bridgeCount))) changed = true;

    return changed;
}

function applyAcceptorLeg(
    bridgeIndex: number,
    bridgeCount: number,
    canonical: CanonicalLegIndices,
    apply: (interval: Interval) => boolean
) {
    let changed = false;
    const i = canonical.acceptor[bridgeIndex];

    if (apply(Interval.ofSingleton(i + 2 * bridgeCount))) changed = true;
    if (apply(Interval.ofSingleton(i + 3 * bridgeCount))) changed = true;

    return changed;
}

function createWaterBridgeCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<WaterBridgeInterUnitParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh);

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return Mesh.createEmpty(mesh);

    const { waterBridges, unitsFeatures } = interactions;

    const n = waterBridges.length;
    if (!n) return Mesh.createEmpty(mesh);

    const l = StructureElement.Location.create(structure);
    const { sizeFactor } = props;
    const canonical = getCanonicalLegIndices(waterBridges);

    const builderProps = {
        // Four half-cylinders per bridge; createLinkCylinderMesh draws the A-side half per call:
        //   [0,   n): donor→water,    forward  (donor side)
        //   [n,  2n): donor→water,    backward (water side)
        //   [2n, 3n): water→acceptor, forward  (water side)
        //   [3n, 4n): water→acceptor, backward (acceptor side)
        //
        // When multiple bridges share the same physical leg, only the first
        // occurrence is drawn. Marking later maps duplicate legs back to the
        // canonical drawn edge index.
        linkCount: 4 * n,

        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const wb = waterBridges[edgeIndex % n];
            const uW = structure.unitMap.get(wb.unitW) as Unit.Atomic;
            const fW = unitsFeatures.get(wb.unitW);
            const leg = Math.floor(edgeIndex / n);

            if (leg === 0) {
                // donor→water, A-side: draw donor→mid
                const uA = structure.unitMap.get(wb.unitA) as Unit.Atomic;
                const fA = unitsFeatures.get(wb.unitA);
                atomPosition(uA, fA, wb.indexA, posA);
                atomPosition(uW, fW, wb.indexWA, posB);
            } else if (leg === 1) {
                // donor→water, B-side: draw water→mid
                const uA = structure.unitMap.get(wb.unitA) as Unit.Atomic;
                const fA = unitsFeatures.get(wb.unitA);
                atomPosition(uW, fW, wb.indexWA, posA);
                atomPosition(uA, fA, wb.indexA, posB);
            } else if (leg === 2) {
                // water→acceptor, A-side: draw water→mid
                const uB = structure.unitMap.get(wb.unitB) as Unit.Atomic;
                const fB = unitsFeatures.get(wb.unitB);
                atomPosition(uW, fW, wb.indexWD, posA);
                atomPosition(uB, fB, wb.indexB, posB);
            } else {
                // water→acceptor, B-side: draw acceptor→mid
                const uB = structure.unitMap.get(wb.unitB) as Unit.Atomic;
                const fB = unitsFeatures.get(wb.unitB);
                atomPosition(uB, fB, wb.indexB, posA);
                atomPosition(uW, fW, wb.indexWD, posB);
            }
        },

        ignore: (edgeIndex: number) => {
            const bi = edgeIndex % n;
            const leg = Math.floor(edgeIndex / n);

            return leg <= 1
                ? canonical.donor[bi] !== bi
                : canonical.acceptor[bi] !== bi;
        },

        style: (_edgeIndex: number) => LinkStyle.Dashed,

        radius: (edgeIndex: number) => {
            const wb = waterBridges[edgeIndex % n];
            const leg = Math.floor(edgeIndex / n);
            const isDonorWaterLeg = leg <= 1;

            if (isDonorWaterLeg) {
                const fA = unitsFeatures.get(wb.unitA);
                const fW = unitsFeatures.get(wb.unitW);

                setFeatureLocation(structure, l, wb.unitA, fA, wb.indexA);
                const sizeA = theme.size.size(l);

                setFeatureLocation(structure, l, wb.unitW, fW, wb.indexWA);
                const sizeW = theme.size.size(l);

                return Math.min(sizeA, sizeW) * sizeFactor;
            } else {
                const fW = unitsFeatures.get(wb.unitW);
                const fB = unitsFeatures.get(wb.unitB);

                setFeatureLocation(structure, l, wb.unitW, fW, wb.indexWD);
                const sizeW = theme.size.size(l);

                setFeatureLocation(structure, l, wb.unitB, fB, wb.indexB);
                const sizeB = theme.size.size(l);

                return Math.min(sizeW, sizeB) * sizeFactor;
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

export const WaterBridgeInterUnitParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    ...InteractionsSharedParams,
};
export type WaterBridgeInterUnitParams = typeof WaterBridgeInterUnitParams

export function WaterBridgeInterUnitVisual(materialId: number): ComplexVisual<WaterBridgeInterUnitParams> {
    return ComplexMeshVisual<WaterBridgeInterUnitParams>({
        defaultProps: PD.getDefaultValues(WaterBridgeInterUnitParams),
        createGeometry: createWaterBridgeCylinderMesh,
        createLocationIterator: createWaterBridgeIterator,
        getLoci: getWaterBridgeLoci,
        eachLocation: eachWaterBridgeInteraction,

        setUpdateState: (
            state: VisualUpdateState,
            newProps: PD.Values<WaterBridgeInterUnitParams>,
            currentProps: PD.Values<WaterBridgeInterUnitParams>,
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

function getWaterBridgeLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id !== objectId) return EmptyLoci;

    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return EmptyLoci;

    const { waterBridges, unitsFeatures } = interactions;
    const n = waterBridges.length;

    if (!n || groupId < 0 || groupId >= 4 * n) return EmptyLoci;

    const bridgeIndex = groupId % n;

    return WaterBridges.Loci({ structure, waterBridges, unitsFeatures }, [{ bridgeIndex }]);
}

const __unitMap = new Map<number, OrderedSet<StructureElement.UnitIndex>>();

function eachWaterBridgeInteraction(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, _isMarking: boolean) {
    let changed = false;

    if (WaterBridges.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;

        const n = interactions.waterBridges.length;
        if (!n) return false;

        const canonical = getCanonicalLegIndices(interactions.waterBridges);

        for (const e of loci.elements) {
            if (e.bridgeIndex < 0 || e.bridgeIndex >= n) continue;

            if (applyDonorLeg(e.bridgeIndex, n, canonical, apply)) changed = true;
            if (applyAcceptorLeg(e.bridgeIndex, n, canonical, apply)) changed = true;
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;

        const { waterBridges, unitsFeatures } = interactions;
        const n = waterBridges.length;
        if (!n) return false;

        const canonical = getCanonicalLegIndices(waterBridges);

        __unitMap.clear();
        for (const e of loci.elements) {
            __unitMap.set(e.unit.id, e.indices);
        }

        for (let i = 0; i < n; i++) {
            const wb = waterBridges[i];

            const indicesA = __unitMap.get(wb.unitA);
            const indicesW = __unitMap.get(wb.unitW);
            const indicesB = __unitMap.get(wb.unitB);

            if (!indicesA && !indicesW && !indicesB) continue;

            let hitA = false;
            if (indicesA) {
                const fA = unitsFeatures.get(wb.unitA);
                const mi = getFeatureMember(fA, wb.indexA);
                hitA = OrderedSet.has(indicesA, mi);
            }

            let hitW = false;
            if (indicesW) {
                const fW = unitsFeatures.get(wb.unitW);
                const miA = getFeatureMember(fW, wb.indexWA);
                const miD = getFeatureMember(fW, wb.indexWD);
                hitW = OrderedSet.has(indicesW, miA) || OrderedSet.has(indicesW, miD);
            }

            let hitB = false;
            if (indicesB) {
                const fB = unitsFeatures.get(wb.unitB);
                const mi = getFeatureMember(fB, wb.indexB);
                hitB = OrderedSet.has(indicesB, mi);
            }

            if (hitA || hitW) {
                if (applyDonorLeg(i, n, canonical, apply)) changed = true;
            }

            if (hitB || hitW) {
                if (applyAcceptorLeg(i, n, canonical, apply)) changed = true;
            }
        }

        __unitMap.clear();
    }

    return changed;
}

function createWaterBridgeIterator(structure: Structure): LocationIterator {
    const interactions = InteractionsProvider.get(structure).value;
    if (!interactions) return LocationIterator(0, 1, 1, () => NullLocation, true);

    const { waterBridges, unitsFeatures } = interactions;

    const n = waterBridges.length;
    const groupCount = 4 * n;
    const instanceCount = 1;

    const data: WaterBridges.Data = { structure, waterBridges, unitsFeatures };
    const location = WaterBridges.Location(data);
    const { element } = location;

    const getLocation = (groupIndex: number) => {
        element.bridgeIndex = n === 0 ? 0 : groupIndex % n;
        return location;
    };

    return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}