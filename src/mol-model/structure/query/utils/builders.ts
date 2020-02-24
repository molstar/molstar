/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement, Structure } from '../../structure';
import { StructureSelection } from '../selection';
import { HashSet } from '../../../../mol-data/generic';
import { structureUnion } from './structure-set';
import { StructureSubsetBuilder } from '../../structure/util/subset-builder';
import { ElementIndex } from '../../model';

export class UniqueStructuresBuilder {
    private set = HashSet(Structure.hashCode, Structure.areUnitAndIndicesEqual);
    private structures: Structure[] = [];
    private allSingletons = true;

    add(s: Structure) {
        if (!s.elementCount) return;
        if (s.elementCount !== 1) this.allSingletons = false;
        if (this.set.add(s)) {
            this.structures[this.structures.length] = s;
        }
    }

    getSelection() {
        if (this.allSingletons) return StructureSelection.Singletons(this.source, structureUnion(this.source, this.structures));
        return StructureSelection.Sequence(this.source, this.structures);
    }

    constructor(private source: Structure) {
    }
}

export class LinearGroupingBuilder {
    private builders: StructureSubsetBuilder[] = [];
    private builderMap = new Map<string, StructureSubsetBuilder>();

    add(key: any, unit: number, element: ElementIndex) {
        let b = this.builderMap.get(key);
        if (!b) {
            b = this.source.subsetBuilder(true);
            this.builders[this.builders.length] = b;
            this.builderMap.set(key, b);
        }
        b.addToUnit(unit, element);
    }

    private allSingletons() {
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            if (this.builders[i].elementCount > 1) return false;
        }
        return true;
    }

    private singletonSelection(): StructureSelection {
        const builder = this.source.subsetBuilder(true);
        const loc = StructureElement.Location.create(this.source);
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            this.builders[i].setSingletonLocation(loc);
            builder.addToUnit(loc.unit.id, loc.element);
        }
        return StructureSelection.Singletons(this.source, builder.getStructure());
    }

    private fullSelection() {
        const structures: Structure[] = new Array(this.builders.length);
        for (let i = 0, _i = this.builders.length; i < _i; i++) {
            structures[i] = this.builders[i].getStructure();
        }
        return StructureSelection.Sequence(this.source, structures);
    }

    getSelection(): StructureSelection {
        const len = this.builders.length;
        if (len === 0) return StructureSelection.Empty(this.source);
        if (this.allSingletons()) return this.singletonSelection();
        return this.fullSelection();
    }

    constructor(private source: Structure) { }
}