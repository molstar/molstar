/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, RuntimeContext } from 'mol-task'
import { Structure } from '../structure'
import Selection from './selection'

// TODO: Query { (s: Structure): Computation<Selection> }

interface Query { (s: Structure): Task<Selection>, provider: Query.Provider }
function Query(q: Query.Provider): Query {
    const ret = (s => Task.create('Query', ctx => q(s, ctx))) as Query;
    ret.provider = q;
    return ret;
}

namespace Query {
    export interface Provider { (s: Structure, ctx: RuntimeContext): Promise<Selection> }
}

export default Query