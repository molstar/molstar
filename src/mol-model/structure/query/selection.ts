/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from 'mol-data/iterator'
import HashSet from 'mol-data/util/hash-set'
import { Structure, Atom, AtomSet } from '../structure'

type Selection =
    | Structure // each atom is interpreted as a singleton structure
    | Structure[]

namespace Selection {
    export const Empty: Selection = [];

    function isStructure(x: Selection): x is Structure { return !!(x as Structure).units && !!(x as Structure).atoms; }

    export function structureCount(sel: Selection) {
        if (isStructure(sel)) return AtomSet.atomCount(sel.atoms);
        return sel.length;
    }

    export function union(sel: Selection): Structure {
        if (isStructure(sel)) return sel;
        if (!sel.length) return Structure.Empty;
        const sets = [];
        for (let i = 0, _i = sel.length; i < _i; i++) sets[sets.length] = sel[i].atoms;
        return Structure.create(unionUnits(sel), AtomSet.unionMany(sets));
    }

    export function structures(sel: Selection): Iterator<Structure> {
        if (isStructure(sel)) {
            const units = sel.units;
            return Iterator.map<Atom, Structure>(AtomSet.atoms(sel.atoms), atoms => Structure.create(units, atoms));
        }
        return Iterator.Array(sel);
    }

    export function getAt(sel: Selection, i: number): Structure {
        if (isStructure(sel)) {
            return Structure.create(sel.units, AtomSet.atomGetAt(sel.atoms, i));
        }
        return sel[i];
    }

    export interface Builder {
        add(s: Structure): void,
        getSelection(): Selection
    }

    class LinearBuilderImpl implements Builder {
        private structures: Structure[] = [];
        private allSingletons = true;

        add(s: Structure) {
            const atomCount = AtomSet.atomCount(s.atoms);
            if (atomCount === 0) return;
            this.structures[this.structures.length] = s;
            if (atomCount !== 1) this.allSingletons = false;
        }

        getSelection() {
            const len = this.structures.length;
            if (len === 0) return Empty;
            if (len === 1) return this.structures[0];
            if (this.allSingletons) return union(this.structures);
            return this.structures;
        }

        constructor() { }
    }

    class HashBuilderImpl implements Builder {
        private structures: Structure[] = [];
        private allSingletons = true;
        private sets = HashSet(AtomSet.hashCode, AtomSet.areEqual);

        add(s: Structure) {
            const atomCount = AtomSet.atomCount(s.atoms);
            if (atomCount === 0 || !this.sets.add(s.atoms)) return;
            this.structures[this.structures.length] = s;
            if (atomCount !== 1) this.allSingletons = false;
        }

        getSelection() {
            const len = this.structures.length;
            if (len === 0) return Empty;
            if (len === 1) return this.structures[0];
            if (this.allSingletons) return union(this.structures);
            return this.structures;
        }

        constructor() { }
    }

    export function LinearBuilder(): Builder { return new LinearBuilderImpl(); }
    export function UniqueBuilder(): Builder { return new HashBuilderImpl(); }

    // TODO: spatial lookup
}

export default Selection

function unionUnits(xs: Structure[]): Structure['units'] {
    let prev = xs[0].units;
    let sameModel = true;
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i].units !== prev) sameModel = false;
    }
    if (sameModel) return prev;

    let ret: any = { ...prev };
    for (let i = 1, _i = xs.length; i < _i; i++) {
        const units = xs[i].units;
        if (units !== prev) {
            const keys = Object.keys(units);
            for (let j = 0; j < keys.length; j++) ret[keys[j]] = (units as any)[keys[j]];
        }
        prev = xs[i];
    }
    return ret;
}
