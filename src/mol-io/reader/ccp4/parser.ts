/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task';
import * as Schema from './schema'
import Result from '../result'
import { FileHandle } from '../../common/file-handle';
import { flipByteOrder } from '../../common/binary';

async function parseInternal(file: FileHandle, ctx: RuntimeContext): Promise<Result<Schema.Ccp4File>> {
    await ctx.update({ message: 'Parsing CCP4 file...' });

    const { buffer } = await file.readBuffer(0, file.length)
    const bin = buffer.buffer

    const intView = new Int32Array(bin, 0, 56)
    const floatView = new Float32Array(bin, 0, 56)
    const dv = new DataView(bin)

    // 53  MAP         Character string 'MAP ' to identify file type
    const MAP = String.fromCharCode(
        dv.getUint8(52 * 4), dv.getUint8(52 * 4 + 1),
        dv.getUint8(52 * 4 + 2), dv.getUint8(52 * 4 + 3)
    )
    if (MAP !== 'MAP ') {
        return Result.error('ccp4 format error, missing "MAP " string');
    }

    // 54  MACHST      Machine stamp indicating machine type which wrote file
    //                 17 and 17 for big-endian or 68 and 65 for little-endian
    const MACHST = [ dv.getUint8(53 * 4), dv.getUint8(53 * 4 + 1) ]

    if (MACHST[ 0 ] === 17 && MACHST[ 1 ] === 17) {
        flipByteOrder(buffer, buffer.length)
    }

    const header: Schema.Ccp4Header = {
        NC: intView[0],
        NR: intView[1],
        NS: intView[2],

        MODE: intView[3],

        NCSTART: intView[4],
        NRSTART: intView[5],
        NSSTART: intView[6],

        NX: intView[7],
        NY: intView[8],
        NZ: intView[9],

        xLength: floatView[10],
        yLength: floatView[11],
        zLength: floatView[12],

        alpha: floatView[13],
        beta: floatView[14],
        gamma: floatView[15],

        MAPC: intView[16],
        MAPR: intView[17],
        MAPS: intView[18],

        AMIN: floatView[19],
        AMAX: floatView[20],
        AMEAN: floatView[21],

        ISPG: intView[22],

        NSYMBT: intView[23],

        LSKFLG: intView[24],

        SKWMAT: [], // TODO bytes 26-34
        SKWTRN: [], // TODO bytes 35-37

        // bytes 50-52 origin in X,Y,Z used for transforms
        originX: floatView[49],
        originY: floatView[50],
        originZ: floatView[51],

        MAP, // bytes 53 MAP
        MACHST, // bytes 54 MACHST

        ARMS: floatView[54],

        // TODO bytes 56 NLABL
        // TODO bytes 57-256 LABEL
    }

    let values
    if (header.MODE === 2) {
        values = new Float32Array(
            bin, 256 * 4 + header.NSYMBT,
            header.NX * header.NY * header.NZ
        )
    } else if (header.MODE === 0) {
        values = new Float32Array(new Int8Array(
            bin, 256 * 4 + header.NSYMBT,
            header.NX * header.NY * header.NZ
        ))
    } else {
        return Result.error(`ccp4 mode '${header.MODE}' unsupported`);
    }

    const result: Schema.Ccp4File = { header, values };
    return Result.success(result);
}

export function parse(file: FileHandle) {
    return Task.create<Result<Schema.Ccp4File>>('Parse CCP4', async ctx => {
        return await parseInternal(file, ctx);
    });
}

export default parse;