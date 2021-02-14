/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeMsgPack } from '../../mol-io/common/msgpack/decode';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { inflate } from '../../mol-util/zip/zip';

export interface G3dHeader {
    magic: 'G3D',
    version: number,
    genome: string,
    name: string,
    offsets: { [resolution: string]: { offset: number, size: number } },
    resolutions: number[]
}

export type G3dDataBlock = {
    header: G3dHeader,
    resolution: number,
    data: {
        [haplotype: string]: {
            [ch: string]: {
                start: number[]
                x: number[],
                y: number[],
                z: number[],
            }
        }
    }
}

const HEADER_SIZE = 64000;

export async function getG3dHeader(ctx: PluginContext, urlOrData: string | Uint8Array): Promise<G3dHeader> {
    const data: Uint8Array = await getRawData(ctx, urlOrData, { offset: 0, size: HEADER_SIZE });
    let last = data.length - 1;
    for (; last >= 0; last--) {
        if (data[last] !== 0) break;
    }
    const header = decodeMsgPack(data.slice(0, last + 1));
    return header;
}

export async function getG3dDataBlock(ctx: PluginContext, header: G3dHeader, urlOrData: string | Uint8Array, resolution: number): Promise<G3dDataBlock> {
    if (!header.offsets[resolution]) throw new Error(`Resolution ${resolution} not available.`);
    const data = await getRawData(ctx, urlOrData, header.offsets[resolution]);
    const unzipped = await ctx.runTask(Task.create('Unzip', ctx => inflate(ctx, data)));

    return {
        header,
        resolution,
        data: decodeMsgPack(unzipped)
    };
}

async function getRawData(ctx: PluginContext, urlOrData: string | Uint8Array, range: { offset: number, size: number }) {
    if (typeof urlOrData === 'string') {
        return await ctx.runTask(ctx.fetch({ url: urlOrData, headers: [['Range', `bytes=${range.offset}-${range.offset + range.size - 1}`]], type: 'binary' }));
    } else {
        return urlOrData.slice(range.offset, range.offset + range.size);
    }
}