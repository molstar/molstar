/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Tuple } from '../../../mol-base/collections/integer'
import Unit from './unit'
import Structure from './structure'

/** Atom pointer */
interface Atom { '@type': Tuple['@type'] }

namespace Atom {
    export const Zero: Atom = Tuple.Zero;
    export const create: (unit: number, index: number) => Atom = Tuple.create;
    export const is: (x: any) => x is Atom = Tuple.is;
    export const unit: (a: Atom) => number = Tuple.fst;
    export const index: (a: Atom) => number = Tuple.snd;
    export const areEqual: (a: Atom, b: Atom) => boolean = Tuple.areEqual;
    export const hashCode: (a: Atom) => number = Tuple.hashCode;

    /** All the information required to access atom properties */
    export interface Location { unit: Unit, atom: number }
    export function Location(): Location { return { unit: {} as any, atom: 0 }; }
    export interface Property<T> { (location: Location): T }
    export interface Predicate extends Property<boolean> { }

    export function setLocation(structure: Structure, l: Location, atom: Atom) {
        l.unit = structure.units[unit(atom)];
        l.atom = index(atom);
    }

    export function property<T>(p: Property<T>) { return p; }
}

export default Atom