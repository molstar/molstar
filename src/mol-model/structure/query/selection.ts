/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HashSet } from 'mol-data/generic'
import { Structure, ElementSet } from '../structure'

// A selection is a pair of a Structure and a sequence of unique AtomSets
type Selection = Selection.Singletons | Selection.Sequence

namespace Selection {
    // If each element of the selection is a singleton, we can use a more efficient representation.
    export interface Singletons { readonly kind: 'singletons', readonly structure: Structure, readonly set: ElementSet }
    export interface Sequence { readonly kind: 'sequence', readonly structure: Structure, readonly sets: ReadonlyArray<ElementSet> }

    export function Singletons(structure: Structure, set: ElementSet): Singletons { return { kind: 'singletons', structure, set } }
    export function Sequence(structure: Structure, sets: ElementSet[]): Sequence { return { kind: 'sequence', structure, sets } }
    export function Empty(structure: Structure): Selection { return Sequence(structure, []); };

    export function isSingleton(s: Selection): s is Singletons { return s.kind === 'singletons'; }
    export function isEmpty(s: Selection) { return isSingleton(s) ? ElementSet.elementCount(s.set) === 0 : s.sets.length === 0; }

    export function structureCount(sel: Selection) {
        if (isSingleton(sel)) return ElementSet.elementCount(sel.set);
        return sel.sets.length;
    }

    export function unionStructure(sel: Selection): Structure {
        if (isEmpty(sel)) return Structure.Empty(sel.structure.units);
        if (isSingleton(sel)) return Structure.create(sel.structure.units, sel.set);
        return Structure.create(sel.structure.units, ElementSet.union(sel.sets, sel.structure.elements));
    }

    export function getAt(sel: Selection, i: number): Structure {
        if (isSingleton(sel)) {
            const atom = ElementSet.elementAt(sel.set, i);
            return Structure.create(sel.structure.units, ElementSet.singleton(atom, sel.structure.elements));
        }
        return Structure.create(sel.structure.units, sel.sets[i]);
    }

    export function toStructures(sel: Selection): Structure[] {
        const { units } = sel.structure;
        if (isSingleton(sel)) {
            const ret: Structure[] = new Array(ElementSet.elementCount(sel.set));
            const atoms = ElementSet.elements(sel.set);
            let offset = 0;
            while (atoms.hasNext) {
                const atom = atoms.move();
                ret[offset++] = Structure.create(units, ElementSet.singleton(atom, sel.structure.elements))
            }
            return ret;
        } else {
            const { sets } = sel;
            const ret: Structure[] = new Array(sets.length);
            for (let i = 0, _i = sets.length; i < _i; i++) ret[i] = Structure.create(units, sets[i]);
            return ret;
        }
    }

    export interface Builder {
        add(set: ElementSet): void,
        getSelection(): Selection
    }

    function getSelection(structure: Structure, sets: ElementSet[], allSingletons: boolean) {
        const len = sets.length;
        if (len === 0) return Empty(structure);
        if (allSingletons) return Singletons(structure, ElementSet.union(sets, structure.elements));
        return Sequence(structure, sets);
    }

    class LinearBuilderImpl implements Builder {
        private sets: ElementSet[] = [];
        private allSingletons = true;

        add(atoms: ElementSet) {
            const atomCount = ElementSet.elementCount(atoms);
            if (atomCount === 0) return;
            this.sets[this.sets.length] = atoms;
            if (atomCount !== 1) this.allSingletons = false;
        }

        getSelection() { return getSelection(this.structure, this.sets, this.allSingletons); }

        constructor(private structure: Structure) { }
    }

    class HashBuilderImpl implements Builder {
        private sets: ElementSet[] = [];
        private allSingletons = true;
        private uniqueSets = HashSet(ElementSet.hashCode, ElementSet.areEqual);

        add(atoms: ElementSet) {
            const atomCount = ElementSet.elementCount(atoms);
            if (atomCount === 0 || !this.uniqueSets.add(atoms)) return;
            this.sets[this.sets.length] = atoms;
            if (atomCount !== 1) this.allSingletons = false;
        }

        getSelection() { return getSelection(this.structure, this.sets, this.allSingletons); }

        constructor(private structure: Structure) { }
    }

    export function LinearBuilder(structure: Structure): Builder { return new LinearBuilderImpl(structure); }
    export function UniqueBuilder(structure: Structure): Builder { return new HashBuilderImpl(structure); }

    // TODO: spatial lookup
}

export default Selection