/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from '../structure'
import Element from '../element'
import { Lookup3D, GridLookup3D, Result, Box3D, Sphere3D } from 'mol-math/geometry';
import { ElementSet, Unit } from '../../structure';
import { Vec3 } from 'mol-math/linear-algebra';
import { OrderedSet } from 'mol-data/int';
import { computeStructureBoundary } from './boundary';

interface ElementSetLookup3D extends Lookup3D<Element> {}

namespace ElementSetLookup3D {
    class Impl implements ElementSetLookup3D {
        private unitLookup: Lookup3D;
        private result = Result.create<Element>();
        private pivot = Vec3.zero();

        find(x: number, y: number, z: number, radius: number): Result<Element> {
            Result.reset(this.result);
            const { units, elements } = this.structure;
            const closeUnits = this.unitLookup.find(x, y, z, radius);
            if (closeUnits.count === 0) return this.result;

            for (let t = 0, _t = closeUnits.count; t < _t; t++) {
                const i = closeUnits.indices[t];
                const unitId = ElementSet.groupUnitIndex(elements, i);
                const group = ElementSet.groupAt(elements, i);
                const unit = units[unitId];
                Vec3.set(this.pivot, x, y, z);
                if (!unit.operator.isIdentity) {
                    Vec3.transformMat4(this.pivot, this.pivot, unit.operator.inverse);
                }
                const groupLookup = Unit.getLookup3d(unit, group);
                const groupResult = groupLookup.find(this.pivot[0], this.pivot[1], this.pivot[2], radius);
                for (let j = 0, _j = groupResult.count; j < _j; j++) {
                    Result.add(this.result, Element.create(unitId, groupResult.indices[j]), groupResult.squaredDistances[j]);
                }
            }

            return this.result;
        }

        check(x: number, y: number, z: number, radius: number): boolean {
            const { units, elements } = this.structure;
            const closeUnits = this.unitLookup.find(x, y, z, radius);
            if (closeUnits.count === 0) return false;

            for (let t = 0, _t = closeUnits.count; t < _t; t++) {
                const i = closeUnits.indices[t];
                const unitId = ElementSet.groupUnitIndex(elements, i);
                const group = ElementSet.groupAt(elements, i);
                const unit = units[unitId];
                Vec3.set(this.pivot, x, y, z);
                if (!unit.operator.isIdentity) {
                    Vec3.transformMat4(this.pivot, this.pivot, unit.operator.inverse);
                }
                const groupLookup = Unit.getLookup3d(unit, group);
                if (groupLookup.check(this.pivot[0], this.pivot[1], this.pivot[2], radius)) return true;
            }

            return false;
        }

        boundary: { box: Box3D; sphere: Sphere3D; };

        constructor(private structure: Structure) {
            const { units, elements } = structure;
            const unitCount = ElementSet.groupCount(elements);
            const xs = new Float32Array(unitCount);
            const ys = new Float32Array(unitCount);
            const zs = new Float32Array(unitCount);
            const radius = new Float32Array(unitCount);

            const center = Vec3.zero();
            for (let i = 0; i < unitCount; i++) {
                const group = ElementSet.groupAt(elements, i);
                const unit = units[ElementSet.groupUnitIndex(elements, i)];
                const lookup = Unit.getLookup3d(unit, group);
                const s = lookup.boundary.sphere;

                Vec3.transformMat4(center, s.center, unit.operator.matrix);

                xs[i] = center[0];
                ys[i] = center[1];
                zs[i] = center[2];
                radius[i] = s.radius;
            }

            this.unitLookup = GridLookup3D({ x: xs, y: ys, z: zs, radius, indices: OrderedSet.ofBounds(0, unitCount) });
            this.boundary = computeStructureBoundary(structure);
        }
    }

    export function create(s: Structure): ElementSetLookup3D {
        return new Impl(s);
    }
}

export { ElementSetLookup3D }