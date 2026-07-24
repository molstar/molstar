/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../../../mol-data/int/ordered-set';
import { SortedArray } from '../../../../mol-data/int/sorted-array';
import { PositionData } from '../../../../mol-math/geometry/common';
import { GridLookup3D } from '../../../../mol-math/geometry/lookup3d/grid';
import { Box3D } from '../../../../mol-math/geometry/primitives/box3d';
import { Sphere3D } from '../../../../mol-math/geometry/primitives/sphere3d';
import { Vec3 } from '../../../../mol-math/linear-algebra/3d/vec3';
import { ElementIndex } from '../../model/indexing';
import { Unit } from '../unit';

/**
 * Merge the elements of units from the same model that are known to share
 * the same symmetry operator into a single unit.
 */
export function mergeUnitsWithSameOperator(units: readonly Unit[], id?: number): Unit {
    const u = units[0];
    if (id === undefined) id = u.id;

    let start = -1 as ElementIndex, end = -1 as ElementIndex;
    let elements = SortedArray.Empty as SortedArray<ElementIndex>;

    for (let i = 0, il = units.length; i < il; ++i) {
        const e = units[i].elements;
        if (SortedArray.isRange(e)) {
            if (end !== -1 && e[0] === end + 1) {
                // extend range
                end = e[e.length - 1];
            } else {
                if (end !== -1) {
                    // pending range
                    elements = SortedArray.union(elements, SortedArray.ofRange(start, end));
                }
                // new range
                start = e[0];
                end = e[e.length - 1];
            }
        } else {
            if (end !== -1) {
                // pending range
                elements = SortedArray.union(elements, SortedArray.ofRange(start, end));
                start = -1 as ElementIndex, end = -1 as ElementIndex;
            }
            elements = SortedArray.union(elements, e);
        }
    }

    if (end !== -1) {
        // pending range
        elements = SortedArray.union(elements, SortedArray.ofRange(start, end));
    }

    return Unit.create(u.id, u.id, u.chainGroupId, u.traits | Unit.Trait.MultiChain, u.kind, u.model, u.conformation.operator, elements);
}

/**
 * Partition the untransformed units of the same model into a grid and merge
 * the units in each grid cell that share the same symmetry operator.
 */
export function partitionUntransformedUnits(units: readonly Unit[], cellSize: number) {
    const unitCount = units.length;
    const mergedUnits: Unit[] = [];

    const box = Box3D.setEmpty(Box3D());
    const x = new Float32Array(unitCount);
    const y = new Float32Array(unitCount);
    const z = new Float32Array(unitCount);
    const indices = OrderedSet.ofBounds(0, unitCount);

    for (let i = 0, il = unitCount; i < il; ++i) {
        const v = units[i].boundary.sphere.center;
        x[i] = v[0];
        y[i] = v[1];
        z[i] = v[2];
        Box3D.add(box, v);
    }
    Box3D.expand(box, box, Vec3.create(1, 1, 1));

    const positionData: PositionData = { x, y, z, indices };
    const boundary = { box, sphere: Sphere3D.fromBox3D(Sphere3D(), box) };
    const lookup = GridLookup3D(positionData, boundary, Vec3.create(cellSize, cellSize, cellSize));

    const { array, offset, count } = lookup.buckets;

    for (let i = 0, il = offset.length; i < il; ++i) {
        const start = offset[i];
        const size = count[i];
        const cellUnits: Unit[] = [];
        for (let j = start, jl = start + size; j < jl; ++j) {
            cellUnits.push(units[array[j]]);
        }
        mergedUnits.push(mergeUnitsWithSameOperator(cellUnits, i));
    }

    return mergedUnits;
}
