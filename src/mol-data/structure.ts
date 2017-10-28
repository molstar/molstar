/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

//import { Vec3 } from '../mol-base/math/linear-algebra'
import AtomSet from './structure/atom-set'
import Unit from './structure/unit'
import * as Base from './structure/base'
//import Model from './model'
//import Operator from './structure/operator'

// TODO: do "single model" version of the structure?
export interface Structure extends Readonly<{
    units: Readonly<{ [id: number]: Unit }>,
    atoms: AtomSet
}> { }

export namespace Structure {
    export const Empty = Base.Empty;
    export const ofModel = Base.ofModel;

    // TODO: "lift" atom set operators
    // TODO: "diff"
}

export default Structure