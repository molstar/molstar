/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Tuple, SortedArray } from 'mol-data/int'
import Unit from './unit'
import Structure from './structure'

/** Atom pointer */
interface Element { '@type': Tuple['@type'] }

namespace Element {
    export const Zero: Element = Tuple.Zero;
    export const create: (unit: number, index: number) => Element = Tuple.create;
    export const is: (x: any) => x is Element = Tuple.is;
    export const unitId: (e: Element) => number = Tuple.fst;
    export const elementIndex: (e: Element) => number = Tuple.snd;
    export const areEqual: (e: Element, b: Element) => boolean = Tuple.areEqual;
    export const hashCode: (e: Element) => number = Tuple.hashCode;

    export function createEmptyArray(n: number): Element[] { return new Float64Array(n) as any; }

    /** All the information required to access element properties */
    export interface Location { unit: Unit, element: number }
    export function Location(): Location { return { unit: {} as any, element: 0 }; }
    export interface Property<T> { (location: Location): T }
    export interface Predicate extends Property<boolean> { }

    export function updateLocation(structure: Structure, l: Location, element: Element) {
        l.unit = structure.units[unitId(element)];
        l.element = elementIndex(element);
        return l;
    }

    export function property<T>(p: Property<T>) { return p; }

    /** Represents multiple element locations */
    export type Loci = ReadonlyArray<{ unit: Unit, elements: SortedArray }>
}

export default Element