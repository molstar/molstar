/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from '../mol-base/collections/ordered-set'
import Iterator from '../mol-base/collections/iterator'
import IntTuple from '../mol-base/collections/int-tuple'
import * as Base from './atom-set/base'
import createBuilder from './atom-set/builder'

/** A map-like representation of integer set */
namespace AtomSet {
    export const Empty: AtomSet = Base.Empty as any;

    export const create: (data: IntTuple | ArrayLike<IntTuple> | IntTuple | { [id: number]: OrderedSet }) => AtomSet = Base.create as any;

    export const keys: (set: AtomSet) => OrderedSet = Base.getKeys as any;
    export const keyCount: (set: AtomSet) => number = Base.keyCount as any;
    export const hasKey: (set: AtomSet, key: number) => boolean = Base.hasKey as any;
    export const geyKey: (set: AtomSet, i: number) => number = Base.getKey as any;
    export const getByKey: (set: AtomSet, key: number) => OrderedSet = Base.getByKey as any;
    export const getByIndex: (set: AtomSet, i: number) => OrderedSet = Base.getByIndex as any;
    export const has: (set: AtomSet, x: IntTuple) => boolean = Base.hasTuple as any;
    export const indexOf: (set: AtomSet, x: IntTuple) => number = Base.indexOf as any;
    export const getAt: (set: AtomSet, i: number) => IntTuple = Base.getAt as any;
    export const values: (set: AtomSet) => Iterator<IntTuple.Unpacked> = Base.values as any;

    export const size: (set: AtomSet) => number = Base.size as any;
    export const hashCode: (set: AtomSet) => number = Base.hashCode as any;

    export const areEqual: (a: AtomSet, b: AtomSet) => boolean = Base.areEqual as any;
    export const areIntersecting: (a: AtomSet, b: AtomSet) => boolean = Base.areIntersecting as any;

    export const union: (a: AtomSet, b: AtomSet) => AtomSet = Base.union as any;
    export const unionMany: (sets: AtomSet[]) => AtomSet = Base.unionMany as any;
    export const intersect: (a: AtomSet, b: AtomSet) => AtomSet = Base.intersect as any;
    export const subtract: (a: AtomSet, b: AtomSet) => AtomSet = Base.subtract as any;

    export function SortedBuilder(parent: AtomSet) { return createBuilder(parent, true); }
    export function Builder(parent: AtomSet) { return createBuilder(parent, false); }
}

interface AtomSet { '@type': 'atom-set' }

export default AtomSet