/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, SortedArray, Iterator } from 'mol-base/collections/integer'
import Atom from '../atom'
import * as Impl from './set/impl'
import createBuilder, { Builder as AtomSetBuilder } from './set/builder'

/** A map-like representation of grouped atom set */
namespace AtomSet {
    export const Empty: AtomSet = Impl.Empty as any;

    export const create: (data: Atom | ArrayLike<Atom> | { [unitId: number]: OrderedSet }) => AtomSet = Impl.create as any;

    export const unitCount: (set: AtomSet) => number = Impl.keyCount as any;
    export const unitIds: (set: AtomSet) => SortedArray = Impl.getKeys as any;
    export const unitHas: (set: AtomSet, id: number) => boolean = Impl.hasKey as any;
    export const unitGetId: (set: AtomSet, i: number) => number = Impl.getKey as any;

    export const unitGetById: (set: AtomSet, key: number) => OrderedSet = Impl.getByKey as any;
    export const unitGetByIndex: (set: AtomSet, i: number) => OrderedSet = Impl.getByIndex as any;

    export const atomCount: (set: AtomSet) => number = Impl.size as any;
    export const atomHas: (set: AtomSet, x: Atom) => boolean = Impl.hasAtom as any;
    export const atomIndexOf: (set: AtomSet, x: Atom) => number = Impl.indexOf as any;
    export const atomGetAt: (set: AtomSet, i: number) => Atom = Impl.getAt as any;
    export const atoms: (set: AtomSet) => Iterator<Atom> = Impl.values as any;

    export const hashCode: (set: AtomSet) => number = Impl.hashCode as any;
    export const areEqual: (a: AtomSet, b: AtomSet) => boolean = Impl.areEqual as any;
    export const areIntersecting: (a: AtomSet, b: AtomSet) => boolean = Impl.areIntersecting as any;

    export const union: (a: AtomSet, b: AtomSet) => AtomSet = Impl.union as any;
    export const unionMany: (sets: AtomSet[]) => AtomSet = Impl.unionMany as any;
    export const intersect: (a: AtomSet, b: AtomSet) => AtomSet = Impl.intersect as any;
    export const subtract: (a: AtomSet, b: AtomSet) => AtomSet = Impl.subtract as any;

    export type Builder = AtomSetBuilder
    export function LinearBuilder(parent: AtomSet): Builder { return createBuilder(parent, true); }
    export function UnsortedBuilder(parent: AtomSet): Builder { return createBuilder(parent, false); }

    // TODO: bounding sphere
    // TODO: distance, areWithIn?
    // TODO: check connected
    // TODO: add "parent" property? how to avoid using too much memory? Transitive parents? Parent unlinking?
}

interface AtomSet { '@type': 'atom-set' | Atom['@type'] }

export default AtomSet