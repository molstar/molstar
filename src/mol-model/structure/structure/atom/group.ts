/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int'
import Unit from '../unit'

interface AtomGroup {
    set: OrderedSet,
    id: number
}

namespace AtomGroup {
    export const Empty = createNew(OrderedSet.Empty)

    // export function singleton(idx: number) {
    //     return create(OrderedSet.ofSingleton(idx));
    // }

    export function createNew(set: OrderedSet) {
        return { id: nextId(), set };
    }

    export function create(unit: Unit, set: OrderedSet): AtomGroup {
        if (OrderedSet.areEqual(set, unit.naturalGroup.set)) return unit.naturalGroup;
        return createNew(set);
    }

    // export function size(group: AtomGroup) { return OrderedSet.size(group.set); }
    // export function has(group: AtomGroup, atom: number) { return OrderedSet.has(group.set, atom); }
    // export function getAt(group: AtomGroup, i: number) { return OrderedSet.getAt(group.set, i); }
    // export function indexOf(group: AtomGroup, atom: number) { return OrderedSet.indexOf(group.set, atom); }
    // export function hashCode(group: AtomGroup) { return OrderedSet.hashCode(group.set); }
    // export function areEqual(a: AtomGroup, b: AtomGroup) { return OrderedSet.areEqual(a.set, b.set); }

    let _id = 0;
    function nextId() {
        const ret = _id;
        _id = (_id + 1) % 0x3fffffff;
        return ret;
    }
}

export default AtomGroup