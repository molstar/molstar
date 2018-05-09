/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from '../structure'
import Element from '../element'
import { Lookup3D, GridLookup3D, Result, Box3D, Sphere3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { computeStructureBoundary } from './boundary';
import { OrderedSet } from 'mol-data/int';

interface StructureLookup3D extends Lookup3D<Element> {}

namespace StructureLookup3D {
    class Impl implements StructureLookup3D {
        private unitLookup: Lookup3D;
        private result = Result.create<Element>();
        private pivot = Vec3.zero();

        find(x: number, y: number, z: number, radius: number): Result<Element> {
            Result.reset(this.result);
            const { units } = this.structure;
            const closeUnits = this.unitLookup.find(x, y, z, radius);
            if (closeUnits.count === 0) return this.result;

            for (let t = 0, _t = closeUnits.count; t < _t; t++) {
                const unit = units[closeUnits.indices[t]];
                Vec3.set(this.pivot, x, y, z);
                if (!unit.conformation.operator.isIdentity) {
                    Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
                }
                const unitLookup = unit.lookup3d;
                const groupResult = unitLookup.find(this.pivot[0], this.pivot[1], this.pivot[2], radius);
                for (let j = 0, _j = groupResult.count; j < _j; j++) {
                    Result.add(this.result, Element.create(unit.id, groupResult.indices[j]), groupResult.squaredDistances[j]);
                }
            }

            return this.result;
        }

        check(x: number, y: number, z: number, radius: number): boolean {
            const { units } = this.structure;
            const closeUnits = this.unitLookup.find(x, y, z, radius);
            if (closeUnits.count === 0) return false;

            for (let t = 0, _t = closeUnits.count; t < _t; t++) {
                const unit = units[closeUnits.indices[t]];
                Vec3.set(this.pivot, x, y, z);
                if (!unit.conformation.operator.isIdentity) {
                    Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
                }
                const groupLookup = unit.lookup3d;
                if (groupLookup.check(this.pivot[0], this.pivot[1], this.pivot[2], radius)) return true;
            }

            return false;
        }

        boundary: { box: Box3D; sphere: Sphere3D; };

        constructor(private structure: Structure) {
            const { units } = structure;
            const unitCount = units.length;
            const xs = new Float32Array(unitCount);
            const ys = new Float32Array(unitCount);
            const zs = new Float32Array(unitCount);
            const radius = new Float32Array(unitCount);

            const center = Vec3.zero();
            for (let i = 0; i < unitCount; i++) {
                const unit = units[i];
                const lookup = unit.lookup3d;
                const s = lookup.boundary.sphere;

                Vec3.transformMat4(center, s.center, unit.conformation.operator.matrix);

                xs[i] = center[0];
                ys[i] = center[1];
                zs[i] = center[2];
                radius[i] = s.radius;
            }

            this.unitLookup = GridLookup3D({ x: xs, y: ys, z: zs, radius, indices: OrderedSet.ofBounds(0, unitCount) });
            this.boundary = computeStructureBoundary(structure);
        }
    }

    export function create(s: Structure): StructureLookup3D {
        return new Impl(s);
    }
}

export { StructureLookup3D }