/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray, Iterator, OrderedSet } from 'mol-data/int'
import Element from '../element'
import ElementGroup from './group'
import * as Impl from './impl/set'
import * as Builders from './impl/set-builder'
import { StructureLookup3D } from '../util/lookup3d';
import Structure from '../structure';

/**
 * A map-like representation of grouped atom set
 *
 * Essentially corresponds to the type { [unitId: number]: ElementGroup }.
 */
namespace ElementSet {
    export const Empty: ElementSet = Impl.Empty as any;

    export const ofAtoms: (elements: ArrayLike<Element>, template: ElementSet) => ElementSet = Impl.ofElements as any;
    export const singleton: (element: Element, template: ElementSet) => ElementSet = Impl.singleton as any;

    export const unitIndices: (set: ElementSet) => SortedArray = Impl.getKeys as any;
    export const unitHas: (set: ElementSet, index: number) => boolean = Impl.hasKey as any;

    export const groupCount: (set: ElementSet) => number = Impl.keyCount as any;
    export const groupUnitIndex: (set: ElementSet, index: number) => number = Impl.getKey as any;
    export const groupFromUnitIndex: (set: ElementSet, unitId: number) => ElementGroup = Impl.getByKey as any;
    export const groupAt: (set: ElementSet, index: number) => ElementGroup = Impl.getByIndex as any;

    export const elementCount: (set: ElementSet) => number = Impl.size as any;
    export const elementHas: (set: ElementSet, x: Element) => boolean = Impl.hasAtom as any;
    export const elementIndexOf: (set: ElementSet, x: Element) => number = Impl.indexOf as any;
    export const elementAt: (set: ElementSet, i: number) => Element = Impl.getAt as any;
    export const elements: (set: ElementSet) => Iterator<Element> = Impl.values as any;

    export const hashCode: (set: ElementSet) => number = Impl.hashCode as any;
    export const areEqual: (a: ElementSet, b: ElementSet) => boolean = Impl.areEqual as any;
    export const areIntersecting: (a: ElementSet, b: ElementSet) => boolean = Impl.areIntersecting as any;

    export const union: (sets: ArrayLike<ElementSet>, template: ElementSet) => ElementSet = Impl.unionMany as any;
    export const intersect: (a: ElementSet, b: ElementSet) => ElementSet = Impl.intersect as any;
    export const subtract: (a: ElementSet, b: ElementSet) => ElementSet = Impl.subtract as any;

    export type Builder = Builders.Builder
    export const LinearBuilder = Builders.LinearBuilder
    export const UnsortedBuilder = Builders.UnsortedBuilder

    export interface Generator { add(unit: number, set: ElementGroup): void, getSet(): ElementSet }
    export const Generator: () => Generator = Impl.Generator as any

    export interface TemplateGenerator { add(unit: number, set: OrderedSet): void, getSet(): ElementSet }
    export const TemplateGenerator: (template: ElementSet) => TemplateGenerator = Impl.TemplateGenerator as any

    // TODO: distance, areWithIn?
    // TODO: check connected
    // TODO: add "parent" property? how to avoid using too much memory? Transitive parents? Parent unlinking?
}

interface ElementSet { '@type': 'element-set' }

export default ElementSet