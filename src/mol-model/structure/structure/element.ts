/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Tuple } from 'mol-data/int'
import Unit from './unit'
import Structure from './structure'

/** Atom pointer */
interface Element { '@type': Tuple['@type'] }

namespace Element {
    export const Zero: Element = Tuple.Zero;
    export const create: (unit: number, index: number) => Element = Tuple.create;
    export const is: (x: any) => x is Element = Tuple.is;
    export const unit: (e: Element) => number = Tuple.fst;
    export const index: (e: Element) => number = Tuple.snd;
    export const areEqual: (e: Element, b: Element) => boolean = Tuple.areEqual;
    export const hashCode: (e: Element) => number = Tuple.hashCode;

    export function createEmptyArray(n: number): Element[] { return new Float64Array(n) as any; }

    /** All the information required to access atom properties */
    export interface Location { unit: Unit, atom: number }
    export function Location(): Location { return { unit: {} as any, atom: 0 }; }
    export interface Property<T> { (location: Location): T }
    export interface Predicate extends Property<boolean> { }

    export function updateLocation(structure: Structure, l: Location, atom: Element) {
        l.unit = structure.units[unit(atom)];
        l.atom = index(atom);
        return l;
    }

    export function property<T>(p: Property<T>) { return p; }
}

export default Element