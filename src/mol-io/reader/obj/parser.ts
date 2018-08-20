/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Result from '../result'
import { Task, RuntimeContext } from 'mol-task'
import { Mesh } from 'mol-geo/shape/mesh';

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<Mesh>> {
    // TODO
    const result: Mesh = Mesh.createEmpty();
    return Result.success(result);
}

export function parse(data: string) {
    return Task.create<Result<Mesh>>('Parse OBJ', async ctx => {
        return await parseInternal(data, ctx);
    });
}

export default parse;