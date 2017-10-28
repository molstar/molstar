/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Tuple from '../mol-base/collections/integer/tuple'

interface Atom { '@type': Tuple['@type'] }

namespace Atom {
    export const Zero: Atom = Tuple.Zero;
    export const create: (unit: number, index: number) => Atom = Tuple.create;
    export const is: (x: any) => x is Atom = Tuple.is;
    export const unit: (a: Atom) => number = Tuple.fst;
    export const index: (a: Atom) => number = Tuple.snd;
    export const areEqual: (a: Atom, b: Atom) => boolean = Tuple.areEqual;
    export const hashCode: (a: Atom) => number = Tuple.hashCode;
}

export default Atom