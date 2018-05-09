/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IntMap, SortedArray, Iterator } from 'mol-data/int'
import { UniqueArray } from 'mol-data/generic'
import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import { Model, Format } from '../model'
import { sortArray, sort, arraySwap, hash1 } from 'mol-data/util';
import Element from './element'
import Unit from './unit'
import { StructureLookup3D } from './util/lookup3d';
import StructureSymmetry from './symmetry';

class Structure {
    readonly unitMap: IntMap<Unit>;
    readonly units: ReadonlyArray<Unit>;
    readonly elementCount: number;

    private _hashCode = 0;

    subsetBuilder(isSorted: boolean) {
        return new Structure.SubsetBuilder(this, isSorted);
    }

    get hashCode() {
        if (this._hashCode !== 0) return this._hashCode;
        return this.computeHash();
    }

    private computeHash() {
        let hash = 23;
        for (let i = 0, _i = this.units.length; i < _i; i++) {
            const u = this.units[i];
            hash = (31 * hash + u.id) | 0;
            hash = (31 * hash + SortedArray.hashCode(u.elements)) | 0;
        }
        hash = (31 * hash + this.elementCount) | 0;
        hash = hash1(hash);
        this._hashCode = hash;
        return hash;
    }

    elementLocations(): Iterator<Element.Location> {
        return new Structure.ElementLocationIterator(this);
    }

    get boundary() {
        return this.lookup3d.boundary;
    }

    private _lookup3d?: StructureLookup3D = void 0;
    get lookup3d() {
        if (this._lookup3d) return this._lookup3d;
        this._lookup3d = StructureLookup3D.create(this);
        return this._lookup3d;
    }

    constructor(units: ArrayLike<Unit>) {
        const map = IntMap.Mutable<Unit>();
        let elementCount = 0;
        let isSorted = true;
        let lastId = units.length > 0 ? units[0].id : 0;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const u = units[i];
            map.set(u.id, u);
            elementCount += u.elements.length;
            if (u.id < lastId) isSorted = false;
            lastId = u.id;
        }
        if (!isSorted) sort(units, 0, units.length, cmpUnits, arraySwap)
        this.unitMap = map;
        this.units = units as ReadonlyArray<Unit>;
        this.elementCount = elementCount;
    }
}

function cmpUnits(units: ArrayLike<Unit>, i: number, j: number) { return units[i].id - units[j].id; }

namespace Structure {
    export const Empty = new Structure([]);

    export function create(units: ReadonlyArray<Unit>): Structure { return new Structure(units); }

    export function ofData(format: Format) {
        const models = Model.create(format);
        return models.map(ofModel);
    }

    export function ofModel(model: Model): Structure {
        const chains = model.hierarchy.chainSegments;
        const builder = new StructureBuilder();

        for (let c = 0; c < chains.count; c++) {
            const elements = SortedArray.ofBounds(chains.segments[c], chains.segments[c + 1]);
            builder.addUnit(Unit.Kind.Atomic, model, SymmetryOperator.Default, elements);
        }

        const cs = model.coarseGrained;
        if (cs.isDefined) {
            if (cs.spheres.count > 0) {
                const elements = SortedArray.ofBounds(0, cs.spheres.count);
                builder.addUnit(Unit.Kind.Spheres, model, SymmetryOperator.Default, elements);
            }
            if (cs.gaussians.count > 0) {
                const elements = SortedArray.ofBounds(0, cs.spheres.count);
                builder.addUnit(Unit.Kind.Gaussians, model, SymmetryOperator.Default, elements);
            }
        }

        return builder.getStructure();
    }

    export class StructureBuilder {
        private units: Unit[] = [];

        addUnit(kind: Unit.Kind, model: Model, operator: SymmetryOperator, elements: SortedArray): Unit {
            const unit = Unit.create(this.units.length, kind, model, operator, elements);
            this.units.push(unit);
            return unit;
        }

        addWithOperator(unit: Unit, operator: SymmetryOperator): Unit {
            const newUnit = unit.applyOperator(this.units.length, operator);
            this.units.push(newUnit);
            return newUnit;
        }

        getStructure(): Structure {
            return create(this.units);
        }

        get isEmpty() {
            return this.units.length === 0;
        }
    }

    export function Builder() { return new StructureBuilder(); }

    export class SubsetBuilder {
        private ids: number[] = [];
        private unitMap = IntMap.Mutable<number[]>();
        private parentId = -1;
        private currentUnit: number[] = [];
        elementCount = 0;

        addToUnit(parentId: number, e: number) {
            const unit = this.unitMap.get(parentId);
            if (!!unit) { unit[unit.length] = e; }
            else {
                this.unitMap.set(parentId, [e]);
                this.ids[this.ids.length] = parentId;
            }
            this.elementCount++;
        }

        beginUnit(parentId: number) {
            this.parentId = parentId;
            this.currentUnit = this.currentUnit.length > 0 ? [] : this.currentUnit;
        }

        addElement(e: number) {
            this.currentUnit[this.currentUnit.length] = e;
            this.elementCount++;
        }

        commitUnit() {
            if (this.currentUnit.length === 0) return;
            this.ids[this.ids.length] = this.parentId;
            this.unitMap.set(this.parentId, this.currentUnit);
            this.parentId = -1;
        }

        setUnit(parentId: number, elements: ArrayLike<number>) {
            this.ids[this.ids.length] = parentId;
            this.unitMap.set(parentId, elements as number[]);
            this.elementCount += elements.length;
        }

        getStructure(): Structure {
            if (this.isEmpty) return Structure.Empty;

            const newUnits: Unit[] = [];
            sortArray(this.ids);

            const symmGroups = StructureSymmetry.UnitEquivalenceBuilder();

            for (let i = 0, _i = this.ids.length; i < _i; i++) {
                const id = this.ids[i];
                const parent = this.parent.unitMap.get(id);

                const unit = this.unitMap.get(id);
                const l = unit.length;

                // if the length is the same, just copy the old unit.
                if (unit.length === parent.elements.length) {
                    newUnits[newUnits.length] = parent;
                    symmGroups.add(parent.id, parent);
                    continue;
                }

                if (!this.isSorted && l > 1) sortArray(unit);

                let child = parent.getChild(SortedArray.ofSortedArray(unit));
                const pivot = symmGroups.add(child.id, child);
                if (child !== pivot) child = pivot.applyOperator(child.id, child.conformation.operator, true);
                newUnits[newUnits.length] = child;
            }

            return create(newUnits);
        }

        setSingletonLocation(location: Element.Location) {
            const id = this.ids[0];
            location.unit = this.parent.unitMap.get(id);
            location.element = this.unitMap.get(id)[0];
        }

        get isEmpty() {
            return this.elementCount === 0;
        }

        constructor(private parent: Structure, private isSorted: boolean) {

        }
    }

    export function getModels(s: Structure) {
        const { units } = s;
        const arr = UniqueArray.create<Model['id'], Model>();
        for (const u of units) {
            UniqueArray.add(arr, u.model.id, u.model);
        }
        return arr.array;
    }

    export function getLookup3d(s: Structure): StructureLookup3D {
        return 0 as any;
    }

    export function getBoundary(s: Structure) {
        return getLookup3d(s).boundary;
    }

    export function hashCode(s: Structure) {
        return s.hashCode;
    }

    export function areEqual(a: Structure, b: Structure) {
        if (a.elementCount !== b.elementCount) return false;
        const len = a.units.length;
        if (len !== b.units.length) return false;

        for (let i = 0; i < len; i++) {
            if (a.units[i].id !== b.units[i].id) return false;
        }

        for (let i = 0; i < len; i++) {
            if (!SortedArray.areEqual(a.units[i].elements, b.units[i].elements)) return false;
        }

        return true;
    }

    export class ElementLocationIterator implements Iterator<Element.Location> {
        private current = Element.Location();
        private unitIndex = 0;
        private elements: SortedArray;
        private len = 0;
        private idx = 0;

        hasNext: boolean;
        move(): Element.Location {
            this.current.element = this.elements[this.idx];
            this.advance();
            return this.current;
        }

        private advance() {
            if (this.idx < this.len - 1) {
                this.idx++;
                return;
            }

            this.idx = 0;
            this.unitIndex++;
            if (this.unitIndex >= this.structure.units.length) {
                this.hasNext = false;
                return;
            }

            this.current.unit = this.structure.units[this.unitIndex];
            this.elements = this.current.unit.elements;
            this.len = this.elements.length;
        }

        constructor(private structure: Structure) {
            this.hasNext = structure.elementCount > 0;
            if (this.hasNext) {
                this.elements = structure.units[0].elements;
                this.len = this.elements.length;
                this.current.unit = structure.units[0];
            }
        }
    }
}

export default Structure