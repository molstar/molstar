/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task';
import { Dsn6File, Dsn6Header } from './schema'
import Result from '../result'
import { FileHandle } from '../../common/file-handle';

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
    }
}

function parseDsn6Header(int: Int16Array): Dsn6Header {
    const factor = 1 / int[ 17 ]
    return {
        xStart: int[ 0 ],
        yStart: int[ 1 ],
        zStart: int[ 2 ],
        xExtent: int[ 3 ],
        yExtent: int[ 4 ],
        zExtent: int[ 5 ],
        xRate: int[ 6 ],
        yRate: int[ 7 ],
        zRate: int[ 8 ],
        xlen: int[ 9 ] * factor,
        ylen: int[ 10 ] * factor,
        zlen: int[ 11 ] * factor,
        alpha: int[ 12 ] * factor,
        beta: int[ 13 ] * factor,
        gamma: int[ 14 ] * factor,
        divisor: int[ 15 ] / 100,
        summand: int[ 16 ],
        sigma: undefined
    }
}

async function parseInternal(file: FileHandle, ctx: RuntimeContext): Promise<Result<Dsn6File>> {
    await ctx.update({ message: 'Parsing DSN6/BRIX file...' });

    const { buffer } = await file.readBuffer(0, file.length)
    const bin = buffer.buffer

    const intView = new Int16Array(bin)
    const byteView = new Uint8Array(bin)
    const brixStr = String.fromCharCode.apply(null, byteView.subarray(0, 512))
    const isBrix = brixStr.startsWith(':-)')

    if (!isBrix) {
        // for DSN6, swap byte order when big endian
        if (intView[18] !== 100) {
            for (let i = 0, n = intView.length; i < n; ++i) {
                const val = intView[i]
                intView[i] = ((val & 0xff) << 8) | ((val >> 8) & 0xff)
            }
        }
    }

    const header = isBrix ? parseBrixHeader(brixStr) : parseDsn6Header(intView)
    const { divisor, summand } = header

    const values = new Float32Array(header.xExtent * header.yExtent * header.zExtent)

    let offset = 512
    const xBlocks = Math.ceil(header.xExtent / 8)
    const yBlocks = Math.ceil(header.yExtent / 8)
    const zBlocks = Math.ceil(header.zExtent / 8)

    // loop over blocks
    for (let zz = 0; zz < zBlocks; ++zz) {
        for (let yy = 0; yy < yBlocks; ++yy) {
            for (let xx = 0; xx < xBlocks; ++xx) {
                // loop inside block
                for (let k = 0; k < 8; ++k) {
                    const z = 8 * zz + k
                    for (let j = 0; j < 8; ++j) {
                        const y = 8 * yy + j
                        for (let i = 0; i < 8; ++i) {
                            const x = 8 * xx + i
                            // check if remaining slice-part contains values
                            if (x < header.xExtent && y < header.yExtent && z < header.zExtent) {
                                const idx = ((((x * header.yExtent) + y) * header.zExtent) + z)
                                values[ idx ] = (byteView[ offset ] - summand) / divisor
                                ++offset
                            } else {
                                offset += 8 - i
                                break
                            }
                        }
                    }
                }
            }
        }
    }

    const result: Dsn6File = { header, values };
    return Result.success(result);
}

export function parseFile(file: FileHandle) {
    return Task.create<Result<Dsn6File>>('Parse DSN6/BRIX', ctx => parseInternal(file, ctx));
}

export function parse(buffer: Uint8Array) {
    return parseFile(FileHandle.fromBuffer(buffer))
}