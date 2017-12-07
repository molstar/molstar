/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, ExecuteSynchronous } from 'mol-task'

async function test() {
    const t = Task.create('test', async () => 1);
    const r = await ExecuteSynchronous(t);
    console.log(r);
}

test();