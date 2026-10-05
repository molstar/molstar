/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { ComplexMeshParams, ComplexMeshVisual } from '../complex-visual.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ComplexVisual } from '../representation.js';
import type { VisualUpdateState } from '../../util.js';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Interval, OrderedSet } from '@molstar/core/data/int';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import type { VisualContext } from '../../visual.js';
import type { Theme } from '@molstar/graphics/theme/theme';
import { convexHull } from '@molstar/core/math/geometry/convex-hull';
import { SortedArray } from '@molstar/core/data/int/sorted-array';
import type { CoordinationIndex } from '@molstar/model/model/structure/structure/coordination/data';

export const CoordinationPolyhedronMeshParams = {
    ...ComplexMeshParams,
    includeParent: PD.Boolean(false),
    minCoordination: PD.Numeric(4, { min: 4, max: 12, step: 1 }, { description: 'Minimum number of coordinating atoms to draw a polyhedron' }),
    maxCoordination: PD.Numeric(12, { min: 4, max: 24, step: 1 }, { description: 'Maximum coordination number' }),
};
export type CoordinationPolyhedronMeshParams = typeof CoordinationPolyhedronMeshParams

function createCoordinationPolyhedronMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CoordinationPolyhedronMeshParams>, mesh?: Mesh) {
    const { minCoordination, maxCoordination } = props;
    const { child, coordination: { sites, eachLigand } } = structure;

    const count = sites.count * 4;
    const builderState = MeshBuilder.createState(count, count / 2, mesh);

    for (let i = 0; i < sites.count; i++) {
        const number = sites.numbers[i];
        if (number < minCoordination || number > maxCoordination) continue;

        if (child) {
            const unit = structure.unitMap.get(sites.unitIds[i]);
            const childUnit = child.unitMap.get(unit.id);
            const element = unit.elements[sites.indices[i]];
            if (!childUnit || !SortedArray.has(childUnit.elements, element)) continue;
        }

        const positions: Vec3[] = [];
        eachLigand(i as CoordinationIndex, l => {
            const p = l.unit.conformation.position(l.element, Vec3());
            positions.push(p);
        });

        const hull = convexHull(positions);
        if (hull) {
            builderState.currentGroup = i;
            for (let i = 0; i < hull.indices.length; i += 3) {
                const a = positions[hull.indices[i]];
                const b = positions[hull.indices[i + 1]];
                const c = positions[hull.indices[i + 2]];
                MeshBuilder.addTriangle(builderState, a, b, c);
            }
        }
    }

    return MeshBuilder.getMesh(builderState);
}

function CoordinationPolyhedronIterator(structure: Structure, props: PD.Values<CoordinationPolyhedronMeshParams>): LocationIterator {
    const { sites } = structure.coordination;

    const groupCount = sites.count;
    const instanceCount = 1;
    const location = StructureElement.Location.create(structure);

    function getLocation(groupIndex: number) {
        if (groupIndex < sites.count) {
            const u = structure.unitMap.get(sites.unitIds[groupIndex]);
            location.unit = u;
            location.element = u.elements[sites.indices[groupIndex]];
        }
        return location;
    }

    return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
}

function getCoordinationPolyhedronLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        if (groupId === PickingId.Null) {
            return Structure.Loci(structure);
        }
        const { sites } = structure.coordination;
        if (groupId < sites.count) {
            return StructureElement.Loci(structure, [{
                unit: structure.unitMap.get(sites.unitIds[groupId]),
                indices: OrderedSet.ofSingleton(sites.indices[groupId])
            }]);
        }
        return Structure.Loci(structure);
    }
    return EmptyLoci;
}

function eachCoordinationPolyhedron(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    const { getSiteIndex } = structure.coordination;
    for (const { unit, indices } of loci.elements) {
        if (!Unit.isAtomic(unit)) continue;
        OrderedSet.forEach(indices, v => {
            const groupIndex = getSiteIndex(unit, unit.elements[v]);
            if (groupIndex >= 0) {
                if (apply(Interval.ofSingleton(groupIndex))) changed = true;
            }
        });
    }
    return changed;
}

export function CoordinationPolyhedronMeshVisual(materialId: number): ComplexVisual<CoordinationPolyhedronMeshParams> {
    return ComplexMeshVisual<CoordinationPolyhedronMeshParams>({
        defaultProps: PD.getDefaultValues(CoordinationPolyhedronMeshParams),
        createGeometry: createCoordinationPolyhedronMesh,
        createLocationIterator: CoordinationPolyhedronIterator,
        getLoci: getCoordinationPolyhedronLoci,
        eachLocation: eachCoordinationPolyhedron,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CoordinationPolyhedronMeshParams>, currentProps: PD.Values<CoordinationPolyhedronMeshParams>) => {
            state.createGeometry = (
                newProps.minCoordination !== currentProps.minCoordination ||
                newProps.maxCoordination !== currentProps.maxCoordination
            );
        }
    }, materialId);
}
