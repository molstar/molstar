/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Structure from '../structure';
import { Lookup3D, GridLookup3D, Box3D, Sphere3D, Result } from '../../../../mol-math/geometry';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { computeStructureBoundary } from './boundary';
import { OrderedSet } from '../../../../mol-data/int';
import { StructureUniqueSubsetBuilder } from './unique-subset-builder';
import StructureElement from '../element';
import Unit from '../unit';
import { getBoundary } from '../../../../mol-math/geometry/boundary';

export interface StructureResult extends Result<StructureElement.UnitIndex> {
    units: Unit[]
}

export namespace StructureResult {
    export function add(result: StructureResult, unit: Unit, index: StructureElement.UnitIndex, distSq: number) {
        result.indices[result.count] = index;
        result.units[result.count] = unit;
        result.squaredDistances[result.count] = distSq;
        result.count++;
    }

    export function create(): StructureResult {
        return { count: 0, indices: [], units: [], squaredDistances: [] };
    }

    export function copy(out: StructureResult, result: StructureResult) {
        for (let i = 0; i < result.count; ++i) {
            out.indices[i] = result.indices[i];
            out.units[i] = result.units[i];
            out.squaredDistances[i] = result.squaredDistances[i];
        }
        out.count = result.count;
        return out;
    }
}

export class StructureLookup3D {
    private unitLookup: Lookup3D;
    private pivot = Vec3.zero();

    findUnitIndices(x: number, y: number, z: number, radius: number): Result<number> {
        return this.unitLookup.find(x, y, z, radius);
    }

    private result = StructureResult.create();
    find(x: number, y: number, z: number, radius: number): StructureResult {
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
                StructureResult.add(this.result, unit, groupResult.indices[j], groupResult.squaredDistances[j]);
            }
        }
        return this.result;
    }

    findIntoBuilder(x: number, y: number, z: number, radius: number, builder: StructureUniqueSubsetBuilder) {
        const { units } = this.structure;
        const closeUnits = this.unitLookup.find(x, y, z, radius);
        if (closeUnits.count === 0) return;

        for (let t = 0, _t = closeUnits.count; t < _t; t++) {
            const unit = units[closeUnits.indices[t]];
            Vec3.set(this.pivot, x, y, z);
            if (!unit.conformation.operator.isIdentity) {
                Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
            }
            const unitLookup = unit.lookup3d;
            const groupResult = unitLookup.find(this.pivot[0], this.pivot[1], this.pivot[2], radius);
            if (groupResult.count === 0) continue;

            const elements = unit.elements;
            builder.beginUnit(unit.id);
            for (let j = 0, _j = groupResult.count; j < _j; j++) {
                builder.addElement(elements[groupResult.indices[j]]);
            }
            builder.commitUnit();
        }
    }

    findIntoBuilderWithRadius(x: number, y: number, z: number, pivotR: number, maxRadius: number, radius: number, eRadius: StructureElement.Property<number>, builder: StructureUniqueSubsetBuilder) {
        const { units } = this.structure;
        const closeUnits = this.unitLookup.find(x, y, z, radius);
        if (closeUnits.count === 0) return;

        const se = StructureElement.Location.create(this.structure);
        const queryRadius = pivotR + maxRadius + radius;

        for (let t = 0, _t = closeUnits.count; t < _t; t++) {
            const unit = units[closeUnits.indices[t]];
            Vec3.set(this.pivot, x, y, z);
            if (!unit.conformation.operator.isIdentity) {
                Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
            }
            const unitLookup = unit.lookup3d;
            const groupResult = unitLookup.find(this.pivot[0], this.pivot[1], this.pivot[2], queryRadius);
            if (groupResult.count === 0) continue;

            const elements = unit.elements;
            se.unit = unit;
            builder.beginUnit(unit.id);
            for (let j = 0, _j = groupResult.count; j < _j; j++) {
                se.element = elements[groupResult.indices[j]];
                const rr = eRadius(se);
                if (Math.sqrt(groupResult.squaredDistances[j]) - pivotR - rr > radius) continue;
                builder.addElement(elements[groupResult.indices[j]]);
            }
            builder.commitUnit();
        }
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

    _boundary: { box: Box3D; sphere: Sphere3D; } | undefined = void 0;

    get boundary() {
        if (this._boundary) return this._boundary!;
        this._boundary = computeStructureBoundary(this.structure);
        return this._boundary!;
    }

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

        const position = { x: xs, y: ys, z: zs, radius, indices: OrderedSet.ofBounds(0, unitCount) };
        this.unitLookup = GridLookup3D(position, getBoundary(position));
    }
}