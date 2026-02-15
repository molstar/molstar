/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../mol-math/linear-algebra';
import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { ComplexMeshParams, ComplexMeshVisual } from '../complex-visual';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme } from '../../../mol-theme/theme';
import { convexHull } from '../../../mol-math/geometry/convex-hull';
import { SortedArray } from '../../../mol-data/int/sorted-array';

export const PolyhedronMeshParams = {
    ...ComplexMeshParams,
    includeParent: PD.Boolean(false),
    minCoordination: PD.Numeric(4, { min: 4, max: 12, step: 1 }, { description: 'Minimum number of coordinating atoms to draw a polyhedron' }),
    maxCoordination: PD.Numeric(12, { min: 4, max: 24, step: 1 }, { description: 'Maximum coordination number' }),
};
export type PolyhedronMeshParams = typeof PolyhedronMeshParams

function createPolyhedronMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<PolyhedronMeshParams>, mesh?: Mesh) {
    const { minCoordination, maxCoordination } = props;
    const { child, coordination: { sites } } = structure;

    const count = sites.length * 4;
    const builderState = MeshBuilder.createState(count, count / 2, mesh);

    for (let i = 0; i < sites.length; i++) {
        const site = sites[i];
        const coord = site.ligandPositions.length;
        if (coord < minCoordination || coord > maxCoordination) continue;

        if (child) {
            const childUnit = child.unitMap.get(site.unit.id);
            if (!childUnit || !SortedArray.has(childUnit.elements, site.element)) continue;
        }

        const hull = convexHull(site.ligandPositions as Vec3[]);
        if (hull) {
            builderState.currentGroup = i;
            for (let i = 0; i < hull.indices.length; i += 3) {
                const a = site.ligandPositions[hull.indices[i]];
                const b = site.ligandPositions[hull.indices[i + 1]];
                const c = site.ligandPositions[hull.indices[i + 2]];
                MeshBuilder.addTriangle(builderState, a, b, c);
            }
        }
    }

    return MeshBuilder.getMesh(builderState);
}

function PolyhedronIterator(structure: Structure, props: PD.Values<PolyhedronMeshParams>): LocationIterator {
    const { sites } = structure.coordination;

    const groupCount = sites.length;
    const instanceCount = 1;
    const location = StructureElement.Location.create(structure);

    function getLocation(groupIndex: number) {
        if (groupIndex < sites.length) {
            const site = sites[groupIndex];
            location.unit = site.unit;
            location.element = site.element;
        }
        return location;
    }

    return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}

function getPolyhedronLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        if (groupId === PickingId.Null) {
            return Structure.Loci(structure);
        }
        const { sites } = structure.coordination;
        if (groupId < sites.length) {
            const site = sites[groupId];
            const unitIndex = SortedArray.indexOf(site.unit.elements, site.element);
            if (unitIndex >= 0) {
                return StructureElement.Loci(structure, [{
                    unit: site.unit,
                    indices: OrderedSet.ofSingleton(unitIndex as StructureElement.UnitIndex)
                }]);
            }
        }
        return Structure.Loci(structure);
    }
    return EmptyLoci;
}

function eachPolyhedron(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    const { getSiteIndices } = structure.coordination;

    for (const { unit, indices } of loci.elements) {
        if (!Unit.isAtomic(unit)) continue;
        OrderedSet.forEach(indices, v => {
            const element = unit.elements[v];
            for (const groupIndex of getSiteIndices(unit, element)) {
                if (apply(Interval.ofSingleton(groupIndex))) changed = true;
            }
        });
    }

    return changed;
}

export function PolyhedronMeshVisual(materialId: number): ComplexVisual<PolyhedronMeshParams> {
    return ComplexMeshVisual<PolyhedronMeshParams>({
        defaultProps: PD.getDefaultValues(PolyhedronMeshParams),
        createGeometry: createPolyhedronMesh,
        createLocationIterator: PolyhedronIterator,
        getLoci: getPolyhedronLoci,
        eachLocation: eachPolyhedron,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolyhedronMeshParams>, currentProps: PD.Values<PolyhedronMeshParams>) => {
            state.createGeometry = (
                newProps.minCoordination !== currentProps.minCoordination ||
                newProps.maxCoordination !== currentProps.maxCoordination
            );
        }
    }, materialId);
}
