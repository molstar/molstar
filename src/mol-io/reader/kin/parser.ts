/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author ReliaSolve <russ@reliasolve.com>
 */

import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext } from '../../../mol-task';
import { Kinemage } from './schema';
import KinParser from './ngl-based-parser';

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<Kinemage>> {

    const NGLParser = new KinParser(data);
    const kinData = NGLParser.kinemage;
    return Result.success(kinData);
}

export function parseKin(data: string) {
    return Task.create<Result<Kinemage>>('Parse KIN', async ctx => {
        return await parseInternal(data, ctx);
    });
}
