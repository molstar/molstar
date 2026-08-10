/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { StringLike } from '../../common/string-like';
import { Tokenizer } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';


// from https://chanzuckerberg.github.io/cryoet-data-portal/stable/cryoet_data_portal_docsite_data.html#annotations

interface OrientedPoint {
    readonly type: 'orientedPoint'
    readonly location: {
        readonly x: number
        readonly y: number
        readonly z: number
    }
    readonly xyz_rotation_matrix: [
        [number, number, number],
        [number, number, number],
        [number, number, number]
    ]
}

interface Point {
    readonly type: 'point'
    readonly location: {
        readonly x: number
        readonly y: number
        readonly z: number
    }
}

export type CryoEtDataPortalNdjsonRecord = OrientedPoint | Point;

function isSupportedRecord(record: any): record is CryoEtDataPortalNdjsonRecord {
    return record?.type === 'orientedPoint' || record?.type === 'point';
}

export interface CryoEtDataPortalNdjsonFile {
    readonly records: ReadonlyArray<CryoEtDataPortalNdjsonRecord>
}

async function parseInternal(data: StringLike, ctx: RuntimeContext) {
    const tokenizer = Tokenizer(data);
    const records: CryoEtDataPortalNdjsonRecord[] = [];
    let prevPosition = 0;

    while (tokenizer.tokenEnd < tokenizer.length) {
        if (tokenizer.position - prevPosition > 100000 && ctx.shouldUpdate) {
            prevPosition = tokenizer.position;
            await ctx.update({ current: tokenizer.position, max: tokenizer.length });
        }

        const line = Tokenizer.readLine(tokenizer).trim();
        if (!line) continue;

        const parsed = JSON.parse(line);
        if (isSupportedRecord(parsed)) {
            records.push(parsed);
        }
    }

    if (records.length === 0) {
        return Result.error<CryoEtDataPortalNdjsonFile>('No readable CryoET Data Portal ndjson records were found.');
    }

    return Result.success({ records });
}

export function parseCryoEtDataPortalNdjson(data: StringLike) {
    return Task.create<Result<CryoEtDataPortalNdjsonFile>>('Parse CryoET Data Portal ndjson', async ctx => {
        return await parseInternal(data, ctx);
    });
}
