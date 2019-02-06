/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task';
import { Ccp4File, Ccp4Header } from './schema'
import Result from '../result'
import { FileHandle } from '../../common/file-handle';

async function parseInternal(file: FileHandle, ctx: RuntimeContext): Promise<Result<Ccp4File>> {
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
    // found MRC files that don't have the MACHST stamp set and are big-endian
    if (MACHST[0] !== 68 && MACHST[1] !== 65) {
        // flip byte order in-place
        for (let i = 0, il = bin.byteLength; i < il; i += 4) {
            dv.setFloat32(i, dv.getFloat32(i), true)
        }
    }

    const header: Ccp4Header = {
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

    const offset = 256 * 4 + header.NSYMBT
    const count = header.NC * header.NR * header.NS
    let values
    if (header.MODE === 2) {
        values = new Float32Array(bin, offset, count)
    } else if (header.MODE === 0) {
        values = new Int8Array(bin, offset, count)
    } else {
        return Result.error(`ccp4 mode '${header.MODE}' unsupported`);
    }

    // if the file was converted by mapmode2to0 - scale the data
    // based on uglymol (https://github.com/uglymol/uglymol) by Marcin Wojdyr (wojdyr)
    if (intView[39] === -128 && intView[40] === 127) {
        values = new Float32Array(values)
        // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
        const b1 = (header.AMAX - header.AMIN) / 255.0
        const b0 = 0.5 * (header.AMIN + header.AMAX + b1)
        for (let j = 0, jl = values.length; j < jl; ++j) {
            values[j] = b1 * values[j] + b0
        }
    }

    const result: Ccp4File = { header, values };
    return Result.success(result);
}

export function parseFile(file: FileHandle) {
    return Task.create<Result<Ccp4File>>('Parse CCP4', ctx => parseInternal(file, ctx));
}

export function parse(buffer: Uint8Array) {
    const file = FileHandle.fromBuffer(buffer)
    return Task.create<Result<Ccp4File>>('Parse CCP4', ctx => parseInternal(file, ctx));
}