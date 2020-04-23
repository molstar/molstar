/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { Dsn6File, Dsn6Header } from './schema';
import { ReaderResult as Result } from '../result';
import { FileHandle } from '../../common/file-handle';
import { SimpleBuffer } from '../../../mol-io/common/simple-buffer';

export const dsn6HeaderSize = 512;

function parseBrixHeader(str: string): Dsn6Header {
    return {
        xStart: parseInt(str.substr(10, 5)),
        yStart: parseInt(str.substr(15, 5)),
        zStart: parseInt(str.substr(20, 5)),
        xExtent: parseInt(str.substr(32, 5)),
        yExtent: parseInt(str.substr(38, 5)),
        zExtent: parseInt(str.substr(42, 5)),
        xRate: parseInt(str.substr(52, 5)),
        yRate: parseInt(str.substr(58, 5)),
        zRate: parseInt(str.substr(62, 5)),
        xlen: parseFloat(str.substr(73, 10)),
        ylen: parseFloat(str.substr(83, 10)),
        zlen: parseFloat(str.substr(93, 10)),
        alpha: parseFloat(str.substr(103, 10)),
        beta: parseFloat(str.substr(113, 10)),
        gamma: parseFloat(str.substr(123, 10)),
        divisor: parseFloat(str.substr(138, 12)),
        summand: parseInt(str.substr(155, 8)),
        sigma: parseFloat(str.substr(170, 12))
    };
}

function parseDsn6Header(buffer: SimpleBuffer, littleEndian: boolean): Dsn6Header {
    const readInt = littleEndian ? (o: number) => buffer.readInt16LE(o * 2) : (o: number) => buffer.readInt16BE(o * 2);
    const factor = 1 / readInt(17);
    return {
        xStart: readInt(0),
        yStart: readInt(1),
        zStart: readInt(2),
        xExtent: readInt(3),
        yExtent: readInt(4),
        zExtent: readInt(5),
        xRate: readInt(6),
        yRate: readInt(7),
        zRate: readInt(8),
        xlen: readInt(9) * factor,
        ylen: readInt(10) * factor,
        zlen: readInt(11) * factor,
        alpha: readInt(12) * factor,
        beta: readInt(13) * factor,
        gamma: readInt(14) * factor,
        divisor: readInt(15) / 100,
        summand: readInt(16),
        sigma: undefined
    };
}

function getBlocks(header: Dsn6Header) {
    const { xExtent, yExtent, zExtent } = header;
    const xBlocks = Math.ceil(xExtent / 8);
    const yBlocks = Math.ceil(yExtent / 8);
    const zBlocks = Math.ceil(zExtent / 8);
    return { xBlocks, yBlocks, zBlocks };
}

export async function readDsn6Header(file: FileHandle): Promise<{ header: Dsn6Header, littleEndian: boolean }> {
    const { buffer } = await file.readBuffer(0, dsn6HeaderSize);
    const brixStr = String.fromCharCode.apply(null, buffer) as string;
    const isBrix = brixStr.startsWith(':-)');
    const littleEndian = isBrix || buffer.readInt16LE(18 * 2) === 100;
    const header = isBrix ? parseBrixHeader(brixStr) : parseDsn6Header(buffer, littleEndian);
    return { header, littleEndian };
}

export async function parseDsn6Values(header: Dsn6Header, source: Uint8Array, target: Float32Array, littleEndian: boolean) {
    if (!littleEndian) {
        // even though the values are one byte they need to be swapped like they are 2
        SimpleBuffer.flipByteOrderInPlace2(source.buffer);
    }

    const { divisor, summand, xExtent, yExtent, zExtent } = header;
    const { xBlocks, yBlocks, zBlocks } = getBlocks(header);

    let offset = 0;
    // loop over blocks
    for (let zz = 0; zz < zBlocks; ++zz) {
        for (let yy = 0; yy < yBlocks; ++yy) {
            for (let xx = 0; xx < xBlocks; ++xx) {
                // loop inside block
                for (let k = 0; k < 8; ++k) {
                    const z = 8 * zz + k;
                    for (let j = 0; j < 8; ++j) {
                        const y = 8 * yy + j;
                        for (let i = 0; i < 8; ++i) {
                            const x = 8 * xx + i;
                            // check if remaining slice-part contains values
                            if (x < xExtent && y < yExtent && z < zExtent) {
                                const idx = ((((x * yExtent) + y) * zExtent) + z);
                                target[idx] = (source[offset] - summand) / divisor;
                                ++offset;
                            } else {
                                offset += 8 - i;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

export function getDsn6Counts(header: Dsn6Header) {
    const { xExtent, yExtent, zExtent } = header;
    const { xBlocks, yBlocks, zBlocks } = getBlocks(header);
    const valueCount = xExtent * yExtent * zExtent;
    const count = xBlocks * 8 * yBlocks * 8 * zBlocks * 8;
    const elementByteSize = 1;
    const byteCount = count * elementByteSize;
    return { count, byteCount, valueCount };
}

async function parseInternal(file: FileHandle, size: number, ctx: RuntimeContext): Promise<Dsn6File> {
    await ctx.update({ message: 'Parsing DSN6/BRIX file...' });
    const { header, littleEndian } = await readDsn6Header(file);
    const { buffer } = await file.readBuffer(dsn6HeaderSize, size - dsn6HeaderSize);
    const { valueCount } = getDsn6Counts(header);

    const values = new Float32Array(valueCount);
    await parseDsn6Values(header, buffer, values, littleEndian);

    const result: Dsn6File = { header, values, name: file.name };
    return result;
}

export function parseFile(file: FileHandle, size: number) {
    return Task.create<Result<Dsn6File>>('Parse DSN6/BRIX', async ctx => {
        try {
            return Result.success(await parseInternal(file, size, ctx));
        } catch (e) {
            return Result.error(e);
        }
    });
}

export function parse(buffer: Uint8Array, name: string) {
    return parseFile(FileHandle.fromBuffer(SimpleBuffer.fromUint8Array(buffer), name), buffer.length);
}