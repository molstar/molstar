/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IntMap, SortedArray } from '../../../../mol-data/int';
import { sortArray } from '../../../../mol-data/util';
import StructureElement from '../element';
import StructureSymmetry from '../symmetry';
import Unit from '../unit';
import Structure from '../structure';
import { ElementIndex } from '../../model';

export class StructureSubsetBuilder {
    private ids: number[] = [];
    private unitMap = IntMap.Mutable<ElementIndex[]>();
    private parentId = -1;
    private currentUnit: ElementIndex[] = [];
    elementCount = 0;

    addToUnit(parentId: number, e: ElementIndex) {
        const unit = this.unitMap.get(parentId);
        if (!!unit) {
            unit[unit.length] = e;
        } else {
            this.unitMap.set(parentId, [e]);
            this.ids[this.ids.length] = parentId;
        }
        this.elementCount++;
    }

    beginUnit(parentId: number) {
        this.parentId = parentId;
        this.currentUnit = this.currentUnit.length > 0 ? [] : this.currentUnit;
    }

    addElement(e: ElementIndex) {
        this.currentUnit[this.currentUnit.length] = e;
        this.elementCount++;
    }

    commitUnit() {
        if (this.currentUnit.length === 0) return;
        this.ids[this.ids.length] = this.parentId;
        this.unitMap.set(this.parentId, this.currentUnit);
        this.parentId = -1;
    }

    setUnit(parentId: number, elements: ArrayLike<ElementIndex>) {
        this.ids[this.ids.length] = parentId;
        this.unitMap.set(parentId, elements as ElementIndex[]);
        this.elementCount += elements.length;
    }

    private _getStructure(deduplicateElements: boolean): Structure {
        if (this.isEmpty) return Structure.Empty;

        const newUnits: Unit[] = [];
        sortArray(this.ids);

        const symmGroups = StructureSymmetry.UnitEquivalenceBuilder();

        for (let i = 0, _i = this.ids.length; i < _i; i++) {
            const id = this.ids[i];
            const parent = this.parent.unitMap.get(id);

            let unit: ArrayLike<number> = this.unitMap.get(id);
            let sorted = false;

            if (deduplicateElements) {
                if (!this.isSorted) sortArray(unit);
                unit = SortedArray.deduplicate(SortedArray.ofSortedArray(this.currentUnit));
                sorted = true;
            }

            const l = unit.length;

            // if the length is the same, just copy the old unit.
            if (unit.length === parent.elements.length) {
                newUnits[newUnits.length] = parent;
                symmGroups.add(parent.id, parent);
                continue;
            }

            if (!this.isSorted && !sorted && l > 1) sortArray(unit);

            let child = parent.getChild(SortedArray.ofSortedArray(unit));
            const pivot = symmGroups.add(child.id, child);
            if (child !== pivot) child = pivot.applyOperator(child.id, child.conformation.operator, true);
            newUnits[newUnits.length] = child;
        }

        return Structure.create(newUnits, { parent: this.parent });
    }

    getStructure() {
        return this._getStructure(false);
    }

    getStructureDeduplicate() {
        return this._getStructure(true);
    }

    setSingletonLocation(location: StructureElement.Location) {
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
