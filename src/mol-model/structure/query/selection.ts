/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HashSet } from 'mol-data/generic'
import { Structure } from '../structure'
import { SortedArray } from 'mol-data/int';

// A selection is a pair of a Structure and a sequence of unique AtomSets
type Selection = Selection.Singletons | Selection.Sequence

namespace Selection {
    // If each element of the selection is a singleton, we can use a more efficient representation.
    export interface Singletons { readonly kind: 'singletons', readonly source: Structure, readonly structure: Structure }
    export interface Sequence { readonly kind: 'sequence', readonly source: Structure, readonly structures: Structure[] }

    export function Singletons(source: Structure, structure: Structure): Singletons { return { kind: 'singletons', source, structure } }
    export function Sequence(source: Structure, structures: Structure[]): Sequence { return { kind: 'sequence', source, structures } }
    export function Empty(source: Structure): Selection { return Singletons(source, Structure.Empty); };

    export function isSingleton(s: Selection): s is Singletons { return s.kind === 'singletons'; }
    export function isEmpty(s: Selection) { return isSingleton(s) ? s.structure.units.length === 0 : s.structures.length === 0; }

    export function structureCount(sel: Selection) {
        if (isSingleton(sel)) return sel.structure.elementCount;
        return sel.structures.length;
    }

    export function unionStructure(sel: Selection): Structure {
        if (isEmpty(sel)) return Structure.Empty;
        if (isSingleton(sel)) return sel.structure;
        return union(sel.source, sel.structures);
    }

    export interface Builder {
        add(structure: Structure): void,
        getSelection(): Selection
    }

    function getSelection(source: Structure, structures: Structure[], allSingletons: boolean) {
        const len = structures.length;
        if (len === 0) return Empty(source);
        if (allSingletons) return Singletons(source, union(source, structures));
        return Sequence(source, structures);
    }

    class LinearBuilderImpl implements Builder {
        private structures: Structure[] = [];
        private allSingletons = true;

        add(structure: Structure) {
            const elementCount = structure.elementCount;
            if (elementCount === 0) return;
            this.structures[this.structures.length] = structure;
            if (elementCount !== 1) this.allSingletons = false;
        }

        getSelection() { return getSelection(this.source, this.structures, this.allSingletons); }

        constructor(private source: Structure) { }
    }

    class HashBuilderImpl implements Builder {
        private structures: Structure[] = [];
        private allSingletons = true;
        private uniqueSets = HashSet(Structure.hashCode, Structure.areEqual);

        add(structure: Structure) {
            const atomCount = structure.elementCount;
            if (atomCount === 0 || !this.uniqueSets.add(structure)) return;
            this.structures[this.structures.length] = structure;
            if (atomCount !== 1) this.allSingletons = false;
        }

        getSelection() { return getSelection(this.structure, this.structures, this.allSingletons); }

        constructor(private structure: Structure) { }
    }

    export function LinearBuilder(structure: Structure): Builder { return new LinearBuilderImpl(structure); }
    export function UniqueBuilder(structure: Structure): Builder { return new HashBuilderImpl(structure); }

    // TODO: spatial lookup

    function union(source: Structure, structures: Structure[]) {
        if (structures.length === 0) return Structure.Empty;
        if (structures.length === 1) return structures[0];

        const unitMap = new Map<number, SortedArray>();
        const fullUnits = new Set<number>();

        for (const { units } of structures) {
            for (let i = 0, _i = units.length; i < _i; i++) {
                const u = units[i];
                if (unitMap.has(u.id)) {
                    // check if there is anything more to union in this particual unit.
                    if (fullUnits.has(u.id)) continue;
                    const merged = SortedArray.union(unitMap.get(u.id)!, u.elements);
                    unitMap.set(u.id, merged);
                    if (merged.length === source.unitMap.get(u.id).elements.length) fullUnits.add(u.id);
                } else {
                    unitMap.set(u.id, u.elements);
                    if (u.elements.length === source.unitMap.get(u.id).elements.length) fullUnits.add(u.id);
                }
            }
        }

        const builder = source.subsetBuilder(true);
        unitMap.forEach(buildUnion, builder);
        return builder.getStructure();
    }

    function buildUnion(this: Structure.SubsetBuilder, elements: SortedArray, id: number) {
        this.setUnit(id, elements);
    }
}

export default Selection