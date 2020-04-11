/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IntMap, SortedArray } from '../../../../mol-data/int';
import { sortArray } from '../../../../mol-data/util';
import StructureSymmetry from '../symmetry';
import Unit from '../unit';
import Structure from '../structure';
import { UniqueArray } from '../../../../mol-data/generic';

type UArray = UniqueArray<number, number>

export class StructureUniqueSubsetBuilder {
    private ids: number[] = [];
    private unitMap = IntMap.Mutable<UArray>();
    private parentId = -1;
    private currentUnit: UArray = UniqueArray.create();
    elementCount = 0;

    addToUnit(parentId: number, e: number) {
        const unit = this.unitMap.get(parentId);
        if (!!unit) {
            if (UniqueArray.add(unit, e, e)) this.elementCount++;
        } else {
            const arr: UArray = UniqueArray.create();
            UniqueArray.add(arr, e, e);
            this.unitMap.set(parentId, arr);
            this.ids[this.ids.length] = parentId;
            this.elementCount++;
        }
    }

    has(parentId: number, e: number) {
        const unit = this.unitMap.get(parentId);
        if (!unit) return false;
        return UniqueArray.has(unit, e);
    }

    beginUnit(parentId: number) {
        this.parentId = parentId;
        if (this.unitMap.has(parentId)) {
            this.currentUnit = this.unitMap.get(parentId);
        } else {
            this.currentUnit = this.currentUnit.array.length > 0 ? UniqueArray.create() : this.currentUnit;
        }
    }

    addElement(e: number) {
        if (UniqueArray.add(this.currentUnit, e, e)) this.elementCount++;
    }

    commitUnit() {
        if (this.currentUnit.array.length === 0 || this.unitMap.has(this.parentId)) return;
        this.ids[this.ids.length] = this.parentId;
        this.unitMap.set(this.parentId, this.currentUnit);
        this.parentId = -1;
    }

    getStructure(): Structure {
        if (this.isEmpty) return Structure.Empty;

        const newUnits: Unit[] = [];
        sortArray(this.ids);

        const symmGroups = StructureSymmetry.UnitEquivalenceBuilder();

        for (let i = 0, _i = this.ids.length; i < _i; i++) {
            const id = this.ids[i];
            const parent = this.parent.unitMap.get(id);

            let unit: ArrayLike<number> = this.unitMap.get(id).array;

            const l = unit.length;

            // if the length is the same, just copy the old unit.
            if (unit.length === parent.elements.length) {
                newUnits[newUnits.length] = parent;
                symmGroups.add(parent.id, parent);
                continue;
            }

            if (l > 1) sortArray(unit);

            let child = parent.getChild(SortedArray.ofSortedArray(unit));
            const pivot = symmGroups.add(child.id, child);
            if (child !== pivot) child = pivot.applyOperator(child.id, child.conformation.operator, true);
            newUnits[newUnits.length] = child;
        }

        return Structure.create(newUnits, { parent: this.parent });
    }

    get isEmpty() {
        return this.elementCount === 0;
    }

    constructor(private parent: Structure) {

    }
}
