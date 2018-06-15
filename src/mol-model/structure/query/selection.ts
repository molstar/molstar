/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HashSet } from 'mol-data/generic'
import { Structure, Element, Unit } from '../structure'
import { structureUnion } from './utils/structure';
import { OrderedSet, SortedArray } from 'mol-data/int';

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
        return structureUnion(sel.source, sel.structures);
    }

    export function toLoci(sel: Selection): Element.Loci {
        const loci: { unit: Unit, indices: OrderedSet }[] = [];
        const { unitMap } = sel.source;

        for (const unit of unionStructure(sel).units) {
            if (unit === unitMap.get(unit.id)) {
                loci[loci.length] = { unit, indices: OrderedSet.ofBounds(0, unit.elements.length) };
            } else {
                loci[loci.length] = {
                    unit,
                    indices: OrderedSet.ofSortedArray(SortedArray.indicesOf(sel.source.unitMap.get(unit.id).elements, unit.elements))
                };
            }
        }

        return Element.Loci(loci);
    }

    export interface Builder {
        add(structure: Structure): void,
        getSelection(): Selection
    }

    function getSelection(source: Structure, structures: Structure[], allSingletons: boolean) {
        const len = structures.length;
        if (len === 0) return Empty(source);
        if (allSingletons) return Singletons(source, structureUnion(source, structures));
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
}

export default Selection