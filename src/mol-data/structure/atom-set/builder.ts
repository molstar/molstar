/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import AtomSet from '../atom-set'
import OrderedSet from '../../../mol-base/collections/integer/ordered-set'
import { sortArray } from '../../../mol-base/collections/sort'

class Builder {
    private keys: number[] = [];
    private units: number[][] = Object.create(null);
    private currentUnit: number[] = [];

    add(u: number, a: number) {
        const unit = this.units[u];
        if (!!unit) { unit[unit.length] = a; }
        else {
            this.units[u] = [a];
            this.keys[this.keys.length] = u;
        }
    }

    beginUnit() { this.currentUnit = this.currentUnit.length > 0 ? [] : this.currentUnit; }
    addToUnit(a: number) { this.currentUnit[this.currentUnit.length] = a; }
    commitUnit(u: number) {
        if (this.currentUnit.length === 0) return;
        this.keys[this.keys.length] = u;
        this.units[u] = this.currentUnit;
    }

    getSet(): AtomSet {
        const sets: { [key: number]: OrderedSet } = Object.create(null);

        for (let i = 0, _i = this.keys.length; i < _i; i++) {
            const k = this.keys[i];
            const unit = this.units[k];
            const l = unit.length;
            if (!this.sorted && l > 1) sortArray(unit);
            if (l === 1) {
                sets[k] = OrderedSet.ofSingleton(unit[0]);
            } else {
                const set = OrderedSet.ofSortedArray(unit);
                const parentSet = AtomSet.unitGetById(this.parent, k);
                sets[k] = OrderedSet.areEqual(set, parentSet) ? parentSet : set;
            }
        }
        return AtomSet.create(sets);
    }

    constructor(private parent: AtomSet, private sorted: boolean) { }
}

export default function createBuilder(parent: AtomSet, sorted: boolean) {
    return new Builder(parent, sorted);
}