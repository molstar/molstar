/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int'
import Unit from '../unit'

interface AtomGroup {
    atoms: OrderedSet,
    id: number
}

namespace AtomGroup {
    export const Empty = createNew(OrderedSet.Empty)

    export function singleton(idx: number) {
        return createNew(OrderedSet.ofSingleton(idx));
    }

    export function createNew(atoms: OrderedSet): AtomGroup {
        return { id: nextId(), atoms };
    }

    export function create(unit: Unit, atoms: OrderedSet): AtomGroup {
        if (OrderedSet.areEqual(atoms, unit.naturalGroup.atoms)) return unit.naturalGroup;
        return createNew(atoms);
    }

    export function size(group: AtomGroup) { return OrderedSet.size(group.atoms); }
    export function has(group: AtomGroup, atom: number) { return OrderedSet.has(group.atoms, atom); }
    export function getAt(group: AtomGroup, i: number) { return OrderedSet.getAt(group.atoms, i); }
    export function indexOf(group: AtomGroup, atom: number) { return OrderedSet.indexOf(group.atoms, atom); }
    export function hashCode(group: AtomGroup) { return OrderedSet.hashCode(group.atoms); }
    export function areEqual(a: AtomGroup, b: AtomGroup) { return OrderedSet.areEqual(a.atoms, b.atoms); }

    export function intersect(a: AtomGroup, b: AtomGroup) {
        const set = OrderedSet.intersect(a.atoms, b.atoms);
        if (set === a.atoms) return a;
        if (set === b.atoms) return b;
        return createNew(set);
    }

    export function union(a: AtomGroup, b: AtomGroup) {
        const set = OrderedSet.union(a.atoms, b.atoms);
        if (set === a.atoms) return a;
        if (set === b.atoms) return b;
        return createNew(set);
    }

    export function subtract(a: AtomGroup, b: AtomGroup) {
        const set = OrderedSet.subtract(a.atoms, b.atoms);
        if (set === a.atoms) return a;
        return createNew(set);
    }

    let _id = 0;
    function nextId() {
        const ret = _id;
        _id = (_id + 1) % 0x3fffffff;
        return ret;
    }
}

export default AtomGroup