/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray, Iterator, OrderedSet } from 'mol-data/int'
import Atom from '../atom'
import AtomGroup from './group'
import * as Impl from './impl/set'
import * as Builders from './impl/set-builder'

/** A map-like representation of grouped atom set */
namespace AtomSet {
    export const Empty: AtomSet = Impl.Empty as any;

    export const ofAtoms: (atoms: ArrayLike<Atom>, template: AtomSet) => AtomSet = Impl.ofAtoms as any;

    export const unitCount: (set: AtomSet) => number = Impl.keyCount as any;
    export const unitIds: (set: AtomSet) => SortedArray = Impl.getKeys as any;
    export const unitHas: (set: AtomSet, id: number) => boolean = Impl.hasKey as any;
    export const unitGetId: (set: AtomSet, i: number) => number = Impl.getKey as any;

    export const unitGetById: (set: AtomSet, key: number) => AtomGroup = Impl.getByKey as any;
    export const unitGetByIndex: (set: AtomSet, i: number) => AtomGroup = Impl.getByIndex as any;

    export const atomCount: (set: AtomSet) => number = Impl.size as any;
    export const atomHas: (set: AtomSet, x: Atom) => boolean = Impl.hasAtom as any;
    export const atomIndexOf: (set: AtomSet, x: Atom) => number = Impl.indexOf as any;
    export const atomGetAt: (set: AtomSet, i: number) => Atom = Impl.getAt as any;
    export const atoms: (set: AtomSet) => Iterator<Atom> = Impl.values as any;

    export const hashCode: (set: AtomSet) => number = Impl.hashCode as any;
    export const areEqual: (a: AtomSet, b: AtomSet) => boolean = Impl.areEqual as any;
    export const areIntersecting: (a: AtomSet, b: AtomSet) => boolean = Impl.areIntersecting as any;

    export const union: (sets: ArrayLike<AtomSet>, template: AtomSet) => AtomSet = Impl.unionMany as any;
    export const intersect: (a: AtomSet, b: AtomSet) => AtomSet = Impl.intersect as any;
    export const subtract: (a: AtomSet, b: AtomSet) => AtomSet = Impl.subtract as any;

    export type Builder = Builders.Builder
    export const LinearBuilder = Builders.LinearBuilder
    export const UnsortedBuilder = Builders.UnsortedBuilder

    export interface Generator { add(unit: number, set: AtomGroup): void, getSet(): AtomSet }
    export const Generator: () => Generator = Impl.Generator as any

    export interface TemplateGenerator { add(unit: number, set: OrderedSet): void, getSet(): AtomSet }
    export const TemplateGenerator: (template: AtomSet) => TemplateGenerator = Impl.TemplateGenerator as any

    // TODO: bounding sphere
    // TODO: distance, areWithIn?
    // TODO: check connected
    // TODO: add "parent" property? how to avoid using too much memory? Transitive parents? Parent unlinking?
}

interface AtomSet { '@type': 'atom-set' | Atom['@type'] }

export default AtomSet