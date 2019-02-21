/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task';
import { Ccp4File, Ccp4Header } from './schema'
import { ReaderResult as Result } from '../result'
import { FileHandle } from '../../common/file-handle';
import { SimpleBuffer } from 'mol-io/common/simple-buffer';

export async function readCcp4Header(file: FileHandle): Promise<{ header: Ccp4Header, littleEndian: boolean }> {
    const headerSize = 1024;
    const { buffer } = await file.readBuffer(0, headerSize)

    // 53  MAP         Character string 'MAP ' to identify file type
    const MAP = String.fromCharCode(
        buffer.readUInt8(52 * 4), buffer.readUInt8(52 * 4 + 1),
        buffer.readUInt8(52 * 4 + 2), buffer.readUInt8(52 * 4 + 3)
    )
    if (MAP !== 'MAP ') {
        throw new Error('ccp4 format error, missing "MAP " string');
    }

    // 54  MACHST      Machine stamp indicating machine type which wrote file
    //                 17 and 17 for big-endian or 68 and 65 for little-endian
    const MACHST = [ buffer.readUInt8(53 * 4), buffer.readUInt8(53 * 4 + 1) ]
    let littleEndian = true
    // found MRC files that don't have the MACHST stamp set and are big-endian
    if (MACHST[0] !== 68 && MACHST[1] !== 65) {
        littleEndian = false;
    }

    const readInt = littleEndian ? (o: number) => buffer.readInt32LE(o * 4) : (o: number) => buffer.readInt32BE(o * 4)
    const readFloat = littleEndian ? (o: number) => buffer.readFloatLE(o * 4) : (o: number) => buffer.readFloatBE(o * 4)

    const header: Ccp4Header = {
        NC: readInt(0),
        NR: readInt(1),
        NS: readInt(2),

        MODE: readInt(3),

        NCSTART: readInt(4),
        NRSTART: readInt(5),
        NSSTART: readInt(6),

        NX: readInt(7),
        NY: readInt(8),
        NZ: readInt(9),

        xLength: readFloat(10),
        yLength: readFloat(11),
        zLength: readFloat(12),

        alpha: readFloat(13),
        beta: readFloat(14),
        gamma: readFloat(15),

        MAPC: readInt(16),
        MAPR: readInt(17),
        MAPS: readInt(18),

        AMIN: readFloat(19),
        AMAX: readFloat(20),
        AMEAN: readFloat(21),

        ISPG: readInt(22),

        NSYMBT: readInt(23),

        LSKFLG: readInt(24),

        SKWMAT: [], // TODO bytes 26-34
        SKWTRN: [], // TODO bytes 35-37

        userFlag1: readInt(39),
        userFlag2: readInt(40),

        // bytes 50-52 origin in X,Y,Z used for transforms
        originX: readFloat(49),
        originY: readFloat(50),
        originZ: readFloat(51),

        MAP, // bytes 53 MAP
        MACHST, // bytes 54 MACHST

        ARMS: readFloat(54),

        // TODO bytes 56 NLABL
        // TODO bytes 57-256 LABEL
    }

    return { header, littleEndian }
}

function getElementByteSize(mode: number) {
    switch (mode) {
        case 2: return 4
        case 1: return 2
        case 0: return 1
    }
    throw new Error(`ccp4 mode '${mode}' unsupported`);
}

async function parseInternal(file: FileHandle, size: number, ctx: RuntimeContext): Promise<Ccp4File> {
    await ctx.update({ message: 'Parsing CCP4/MRC file...' });

    const { header, littleEndian } = await readCcp4Header(file)

    const offset = 256 * 4 + header.NSYMBT
    const { buffer, bytesRead } = await file.readBuffer(offset, size - offset)

    const count = header.NC * header.NR * header.NS
    const elementByteSize = getElementByteSize(header.MODE)
    const byteCount = count * elementByteSize

    if (byteCount !== bytesRead) {
        console.warn(`byteCount ${byteCount} and bytesRead ${bytesRead} differ`)
    }

    let values
    if (header.MODE === 2) {
        values = new Float32Array(buffer.buffer, offset, count)
    } else if (header.MODE === 1) {
        values = new Int16Array(buffer.buffer, offset, count)
    } else if (header.MODE === 0) {
        values = new Int8Array(buffer.buffer, offset, count)
    } else {
        throw new Error(`ccp4 mode '${header.MODE}' unsupported`);
    }

    if (!littleEndian) {
        SimpleBuffer.flipByteOrder(buffer, new Uint8Array(values.buffer), byteCount, elementByteSize, 0)
    }

    // if the file was converted by mapmode2to0 - scale the data
    // based on uglymol (https://github.com/uglymol/uglymol) by Marcin Wojdyr (wojdyr)
    if (header.userFlag1 === -128 && header.userFlag2 === 127) {
        values = new Float32Array(values)
        // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
        const b1 = (header.AMAX - header.AMIN) / 255.0
        const b0 = 0.5 * (header.AMIN + header.AMAX + b1)
        for (let j = 0, jl = values.length; j < jl; ++j) {
            values[j] = b1 * values[j] + b0
        }
    }

    const result: Ccp4File = { header, values };
    return result
}

export function parseFile(file: FileHandle, size: number) {
    return Task.create<Result<Ccp4File>>('Parse CCP4/MRC', async ctx => {
        try {
            return Result.success(await parseInternal(file, size, ctx));
        } catch (e) {
            return Result.error(e);
        }
    })
}

export function parse(buffer: Uint8Array) {
    return parseFile(FileHandle.fromBuffer(SimpleBuffer.fromUint8Array(buffer)), buffer.length)
}