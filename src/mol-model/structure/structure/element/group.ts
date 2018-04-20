/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int'
import { Lookup3D } from 'mol-math/geometry'
import Unit from '../unit'

interface ElementGroup {
    elements: OrderedSet,
    // Unique identifier of the group, usable as partial key for various "caches".
    key: number,

    __lookup3d__?: Lookup3D
}

namespace ElementGroup {
    export const Empty = createNew(OrderedSet.Empty)

    export function singleton(idx: number) {
        return createNew(OrderedSet.ofSingleton(idx));
    }

    export function createNew(elements: OrderedSet): ElementGroup {
        return { key: nextKey(), elements };
    }

    export function create(unit: Unit, elements: OrderedSet): ElementGroup {
        if (OrderedSet.areEqual(elements, unit.fullGroup.elements)) return unit.fullGroup;
        return createNew(elements);
    }

    export function createChild(parent: ElementGroup, elements: OrderedSet): ElementGroup {
        if (OrderedSet.areEqual(elements, parent.elements)) return parent;
        return createNew(elements);
    }

    export function size(group: ElementGroup) { return OrderedSet.size(group.elements); }
    export function has(group: ElementGroup, element: number) { return OrderedSet.has(group.elements, element); }
    export function getAt(group: ElementGroup, i: number) { return OrderedSet.getAt(group.elements, i); }
    export function indexOf(group: ElementGroup, element: number) { return OrderedSet.indexOf(group.elements, element); }
    export function hashCode(group: ElementGroup) { return OrderedSet.hashCode(group.elements); }
    export function areEqual(a: ElementGroup, b: ElementGroup) { return OrderedSet.areEqual(a.elements, b.elements); }

    export function intersect(a: ElementGroup, b: ElementGroup) {
        const set = OrderedSet.intersect(a.elements, b.elements);
        if (set === a.elements) return a;
        if (set === b.elements) return b;
        return createNew(set);
    }

    export function union(a: ElementGroup, b: ElementGroup) {
        const set = OrderedSet.union(a.elements, b.elements);
        if (set === a.elements) return a;
        if (set === b.elements) return b;
        return createNew(set);
    }

    export function subtract(a: ElementGroup, b: ElementGroup) {
        const set = OrderedSet.subtract(a.elements, b.elements);
        if (set === a.elements) return a;
        return createNew(set);
    }

    let _id = 0;
    function nextKey() {
        const ret = _id;
        _id = (_id + 1) % 0x3fffffff;
        return ret;
    }
}

export default ElementGroup