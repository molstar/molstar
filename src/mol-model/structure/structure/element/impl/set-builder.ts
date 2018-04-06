/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import ElementSet from '../set'
import Element from '../../element'
import { OrderedSet, IntMap } from 'mol-data/int'
import { sortArray } from 'mol-data/util/sort'

export class Builder {
    private keys: number[] = [];
    private units = IntMap.Mutable<number[]>();
    private currentUnit: number[] = [];

    elementCount = 0;

    add(u: number, e: number) {
        const unit = this.units.get(u);
        if (!!unit) { unit[unit.length] = e; }
        else {
            this.units.set(u, [e]);
            this.keys[this.keys.length] = u;
        }
        this.elementCount++;
    }

    beginUnit() { this.currentUnit = this.currentUnit.length > 0 ? [] : this.currentUnit; }
    addToUnit(a: number) { this.currentUnit[this.currentUnit.length] = a; this.elementCount++; }
    commitUnit(u: number) {
        if (this.currentUnit.length === 0) return;
        this.keys[this.keys.length] = u;
        this.units.set(u, this.currentUnit);
    }

    getSet(): ElementSet {
        const generator = ElementSet.TemplateGenerator(this.parent);

        for (let i = 0, _i = this.keys.length; i < _i; i++) {
            const k = this.keys[i];
            const unit = this.units.get(k);
            const l = unit.length;
            if (!this.sorted && l > 1) sortArray(unit);
            generator.add(k, OrderedSet.ofSortedArray(unit));
        }

        return generator.getSet();
    }

    singleton(): Element {
        const u = this.keys[0];
        return Element.create(u, this.units.get(u)[0]);
    }

    constructor(private parent: ElementSet, private sorted: boolean) { }
}

export function LinearBuilder(parent: ElementSet) {
    return new Builder(parent, true);
}

export function UnsortedBuilder(parent: ElementSet) {
    return new Builder(parent, false);
}