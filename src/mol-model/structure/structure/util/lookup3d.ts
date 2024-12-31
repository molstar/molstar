/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Structure } from '../structure';
import { Lookup3D, GridLookup3D, Result } from '../../../../mol-math/geometry';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { OrderedSet } from '../../../../mol-data/int';
import { StructureUniqueSubsetBuilder } from './unique-subset-builder';
import { StructureElement } from '../element';
import { Unit } from '../unit';
import { UnitIndex } from '../element/util';
import { FibonacciHeap } from '../../../../mol-util/fibonacci-heap';

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

export interface StructureLookup3DResultContext {
    result: StructureResult,
    closeUnitsResult: Result<number>,
    unitGroupResult: Result<UnitIndex>,
}

export function StructureLookup3DResultContext(): StructureLookup3DResultContext {
    return { result: StructureResult.create(), closeUnitsResult: Result.create(), unitGroupResult: Result.create() };
}

export class StructureLookup3D {
    private unitLookup: Lookup3D;
    private pivot = Vec3();
    private heap = new FibonacciHeap();

    findUnitIndices(x: number, y: number, z: number, radius: number): Result<number> {
        return this.unitLookup.find(x, y, z, radius);
    }

    private findContext = StructureLookup3DResultContext();

    find(x: number, y: number, z: number, radius: number, ctx?: StructureLookup3DResultContext): StructureResult {
        return this._find(x, y, z, radius, ctx ?? this.findContext);
    }

    private _find(x: number, y: number, z: number, radius: number, ctx: StructureLookup3DResultContext): StructureResult {
        Result.reset(ctx.result);
        const { units } = this.structure;
        const closeUnits = this.unitLookup.find(x, y, z, radius, ctx.closeUnitsResult);
        if (closeUnits.count === 0) return ctx.result;

        for (let t = 0, _t = closeUnits.count; t < _t; t++) {
            const unit = units[closeUnits.indices[t]];
            Vec3.set(this.pivot, x, y, z);
            if (!unit.conformation.operator.isIdentity) {
                Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
            }
            const unitLookup = unit.lookup3d;
            const groupResult = unitLookup.find(this.pivot[0], this.pivot[1], this.pivot[2], radius, ctx.unitGroupResult);
            for (let j = 0, _j = groupResult.count; j < _j; j++) {
                StructureResult.add(ctx.result, unit, groupResult.indices[j], groupResult.squaredDistances[j]);
            }
        }
        return ctx.result;
    }

    nearest(x: number, y: number, z: number, k: number = 1, ctx?: StructureLookup3DResultContext): StructureResult {
        return this._nearest(x, y, z, k, ctx ?? this.findContext);
    }

    _nearest(x: number, y: number, z: number, k: number, ctx: StructureLookup3DResultContext): StructureResult {
        const result = ctx.result, heap = this.heap;
        Result.reset(result);
        heap.clear();
        const { units } = this.structure;
        let elementsCount = 0;
        const closeUnits = this.unitLookup.nearest(x, y, z, units.length, (uid: number) => (elementsCount += units[uid].elements.length) >= k, ctx.closeUnitsResult); // sort units based on distance to the point
        if (closeUnits.count === 0) return result;
        let totalCount = 0, maxDistResult = -Number.MAX_VALUE;
        for (let t = 0, _t = closeUnits.count; t < _t; t++) {
            const unitSqDist = closeUnits.squaredDistances[t];
            if (totalCount >= k && maxDistResult < unitSqDist) break;
            Vec3.set(this.pivot, x, y, z);
            const unit = units[closeUnits.indices[t]];
            if (!unit.conformation.operator.isIdentity) {
                Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
            }
            const unitLookup = unit.lookup3d;
            const groupResult = unitLookup.nearest(this.pivot[0], this.pivot[1], this.pivot[2], k, void 0, ctx.unitGroupResult);
            if (groupResult.count === 0) continue;
            totalCount += groupResult.count;
            maxDistResult = Math.max(maxDistResult, groupResult.squaredDistances[groupResult.count - 1]);
            for (let j = 0, _j = groupResult.count; j < _j; j++) {
                heap.insert(groupResult.squaredDistances[j], { index: groupResult.indices[j], unit: unit });
            }
        }
        if (k === 1) {
            const node = heap.findMinimum();
            if (node) {
                const { key: squaredDistance } = node;
                const { unit, index } = node.value as { index: UnitIndex, unit: Unit };
                StructureResult.add(result, unit as Unit, index as UnitIndex, squaredDistance as number);
            }
        } else {
            while (!heap.isEmpty() && result.count < k) {
                const node = heap.extractMinimum();
                const { key: squaredDistance } = node!;
                const { unit, index } = node!.value as { index: UnitIndex, unit: Unit };
                StructureResult.add(result, unit as Unit, index as UnitIndex, squaredDistance as number);
            }
        }
        return result;
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

    findIntoBuilderIf(x: number, y: number, z: number, radius: number, builder: StructureUniqueSubsetBuilder, test: (l: StructureElement.Location) => boolean) {
        const { units } = this.structure;
        const closeUnits = this.unitLookup.find(x, y, z, radius);
        if (closeUnits.count === 0) return;

        const loc = StructureElement.Location.create(this.structure);

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
            loc.unit = unit;
            builder.beginUnit(unit.id);
            for (let j = 0, _j = groupResult.count; j < _j; j++) {
                loc.element = elements[groupResult.indices[j]];
                if (test(loc)) {
                    builder.addElement(loc.element);
                }
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

    approxNearest(x: number, y: number, z: number, radius: number, ctx?: StructureLookup3DResultContext): StructureResult {
        return this._approxNearest(x, y, z, radius, ctx ?? this.findContext);
    }

    _approxNearest(x: number, y: number, z: number, radius: number, ctx: StructureLookup3DResultContext): StructureResult {
        Result.reset(ctx.result);
        const { units } = this.structure;
        const closeUnits = this.unitLookup.find(x, y, z, radius, ctx.closeUnitsResult);
        if (closeUnits.count === 0) return ctx.result;

        let minDistSq = Number.MAX_VALUE;
        for (let t = 0, _t = closeUnits.count; t < _t; t++) {
            const unit = units[closeUnits.indices[t]];
            Vec3.set(this.pivot, x, y, z);
            if (!unit.conformation.operator.isIdentity) {
                Vec3.transformMat4(this.pivot, this.pivot, unit.conformation.operator.inverse);
            }
            const unitLookup = unit.lookup3d;
            const groupResult = unitLookup.approxNearest(this.pivot[0], this.pivot[1], this.pivot[2], radius, ctx.unitGroupResult);
            for (let j = 0, _j = groupResult.count; j < _j; j++) {
                if (groupResult.squaredDistances[j] < minDistSq) {
                    StructureResult.add(ctx.result, unit, groupResult.indices[j], groupResult.squaredDistances[j]);
                    minDistSq = groupResult.squaredDistances[j];
                }
            }
        }

        return ctx.result;
    }

    get boundary() {
        return this.structure.boundary;
    }

    constructor(private structure: Structure) {
        const { units, boundary } = structure;
        const unitCount = units.length;
        const xs = new Float32Array(unitCount);
        const ys = new Float32Array(unitCount);
        const zs = new Float32Array(unitCount);
        const radius = new Float32Array(unitCount);

        const center = Vec3();
        for (let i = 0; i < unitCount; i++) {
            const unit = units[i];
            const s = unit.boundary.sphere;

            Vec3.transformMat4(center, s.center, unit.conformation.operator.matrix);

            xs[i] = center[0];
            ys[i] = center[1];
            zs[i] = center[2];
            radius[i] = s.radius;
        }

        const position = { x: xs, y: ys, z: zs, radius, indices: OrderedSet.ofBounds(0, unitCount) };
        this.unitLookup = GridLookup3D(position, boundary);
    }
}
