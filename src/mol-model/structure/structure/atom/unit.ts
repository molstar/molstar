/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO : wrap OrderedSet with and include "id" to identify unique (sub)sets in AtomSet and the ability to cache data (also include "ancestor" field for effictient subsets).

import { OrderedSet } from 'mol-data/int'

interface AtomUnit {
    set: OrderedSet,
    id: number,
    ancestor?: AtomUnit
}

namespace AtomUnit {
    export function create(set: OrderedSet, ancestor?: AtomUnit): AtomUnit {
        if (!ancestor) {
            return { id: nextId(), set };
        }
        if (OrderedSet.areEqual(set, ancestor.set)) return ancestor;
        return { id: nextId(), set, ancestor: ancestor.ancestor || ancestor };
    }

    let _id = 0;
    function nextId() {
        const ret = _id;
        _id = (_id + 1) % 0x3fffffff;
        return ret;
    }
}

export default AtomUnit