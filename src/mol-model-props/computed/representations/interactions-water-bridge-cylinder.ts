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
import { Interval, OrderedSet } from '../../../mol-data/int';
import { InteractionsProvider } from '../interactions';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { WaterBridges } from '../interactions/water-bridges';
import { Sphere3D } from '../../../mol-math/geometry';
import { InteractionsSharedParams } from './shared';
import { Features } from '../interactions/features';

function createWaterBridgeCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<WaterBridgeInterUnitParams>, mesh?: Mesh) {
    if (!structure.hasAtomic) return Mesh.createEmpty(mesh);

    const interactions = InteractionsProvider.get(structure).value!;
    const { waterBridges, unitsFeatures } = interactions;

    const n = waterBridges.length;
    if (!n) return Mesh.createEmpty(mesh);

    const l = StructureElement.Location.create(structure);
    const { sizeFactor } = props;

    const builderProps = {
        // Four half-cylinders per bridge (createLinkCylinderMesh draws only the A-side half per call):
        //   [0,   n): donor→water,    forward  (donor side)
        //   [n,  2n): donor→water,    backward (water side)
        //   [2n, 3n): water→acceptor, forward  (water side)
        //   [3n, 4n): water→acceptor, backward (acceptor side)
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
                // donor→water, B-side: draw water→mid (posA/posB swapped)
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
                // water→acceptor, B-side: draw acceptor→mid (posA/posB swapped)
                const uB = structure.unitMap.get(wb.unitB) as Unit.Atomic;
                const fB = unitsFeatures.get(wb.unitB);
                atomPosition(uB, fB, wb.indexB, posA);
                atomPosition(uW, fW, wb.indexWD, posB);
            }
        },
        style: (_edgeIndex: number) => LinkStyle.Dashed,
        radius: (edgeIndex: number) => {
            const wb = waterBridges[edgeIndex % n];
            const leg = Math.floor(edgeIndex / n);
            const isDonorWaterLeg = leg <= 1;

            const fA = unitsFeatures.get(isDonorWaterLeg ? wb.unitA : wb.unitW);
            const unitIdA = isDonorWaterLeg ? wb.unitA : wb.unitW;
            const indexA = isDonorWaterLeg ? wb.indexA : wb.indexWD;
            l.unit = structure.unitMap.get(unitIdA);
            l.element = l.unit.elements[fA.members[fA.offsets[indexA]]];
            const sizeA = theme.size.size(l);

            const fB = unitsFeatures.get(isDonorWaterLeg ? wb.unitW : wb.unitB);
            const unitIdB = isDonorWaterLeg ? wb.unitW : wb.unitB;
            const indexB = isDonorWaterLeg ? wb.indexWA : wb.indexB;
            l.unit = structure.unitMap.get(unitIdB);
            l.element = l.unit.elements[fB.members[fB.offsets[indexB]]];
            const sizeB = theme.size.size(l);

            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: (_edgeIndex: number) => false,
    };

    const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    if (boundingSphere) {
        m.setBoundingSphere(boundingSphere);
    } else if (m.triangleCount > 0) {
        const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, 1 * sizeFactor);
        m.setBoundingSphere(sphere);
    }

    return m;
}

function atomPosition(unit: Unit.Atomic, features: Features, featureIndex: Features.FeatureIndex, out: Vec3) {
    const atomLocalIdx = features.members[features.offsets[featureIndex]];
    unit.conformation.position(unit.elements[atomLocalIdx], out);
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
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<WaterBridgeInterUnitParams>, currentProps: PD.Values<WaterBridgeInterUnitParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
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

function getWaterBridgeLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id !== objectId) return EmptyLoci;

    const interactions = InteractionsProvider.get(structure).value!;
    const { waterBridges, unitsFeatures } = interactions;
    const bridgeIndex = groupId % waterBridges.length;

    return WaterBridges.Loci({ structure, waterBridges, unitsFeatures }, [{ bridgeIndex }]);
}

const __unitMap = new Map<number, OrderedSet<StructureElement.UnitIndex>>();

function eachWaterBridgeInteraction(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, _isMarking: boolean) {
    let changed = false;

    if (WaterBridges.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
        const interactions = InteractionsProvider.get(structure).value!;
        const n = interactions.waterBridges.length;

        for (const e of loci.elements) {
            // Apply all four half-cylinders for this bridge.
            if (apply(Interval.ofSingleton(e.bridgeIndex))) changed = true;
            if (apply(Interval.ofSingleton(e.bridgeIndex + n))) changed = true;
            if (apply(Interval.ofSingleton(e.bridgeIndex + 2 * n))) changed = true;
            if (apply(Interval.ofSingleton(e.bridgeIndex + 3 * n))) changed = true;
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;

        const interactions = InteractionsProvider.get(structure).value;
        if (!interactions) return false;

        const { waterBridges, unitsFeatures } = interactions;
        const n = waterBridges.length;

        for (const e of loci.elements) __unitMap.set(e.unit.id, e.indices);

        for (let i = 0; i < n; i++) {
            const wb = waterBridges[i];
            const indicesA = __unitMap.get(wb.unitA);
            const indicesB = __unitMap.get(wb.unitB);
            if (!indicesA && !indicesB) continue;

            let hitA = false;
            if (indicesA) {
                const fA = unitsFeatures.get(wb.unitA);
                const mi = fA.members[fA.offsets[wb.indexA]] as StructureElement.UnitIndex;
                hitA = OrderedSet.has(indicesA, mi);
            }

            let hitB = false;
            if (indicesB) {
                const fB = unitsFeatures.get(wb.unitB);
                const mi = fB.members[fB.offsets[wb.indexB]] as StructureElement.UnitIndex;
                hitB = OrderedSet.has(indicesB, mi);
            }

            if (hitA) {
                if (apply(Interval.ofSingleton(i))) changed = true;
                if (apply(Interval.ofSingleton(i + n))) changed = true;
            }
            if (hitB) {
                if (apply(Interval.ofSingleton(i + 2 * n))) changed = true;
                if (apply(Interval.ofSingleton(i + 3 * n))) changed = true;
            }
        }

        __unitMap.clear();
    }

    return changed;
}

function createWaterBridgeIterator(structure: Structure): LocationIterator {
    const interactions = InteractionsProvider.get(structure).value!;
    const { waterBridges, unitsFeatures } = interactions;
    const n = waterBridges.length;
    const groupCount = 4 * n;
    const instanceCount = 1;
    const data: WaterBridges.Data = { structure, waterBridges, unitsFeatures };
    const location = WaterBridges.Location(data);
    const { element } = location;
    const getLocation = (groupIndex: number) => {
        element.bridgeIndex = groupIndex % n;
        return location;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}
