/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { Ccp4File, Ccp4Header } from './schema';
import { ReaderResult as Result } from '../result';
import { FileHandle } from '../../common/file-handle';
import { SimpleBuffer } from '../../../mol-io/common/simple-buffer';
import { TypedArrayValueType, getElementByteSize, TypedArrayBufferContext, readTypedArray, createTypedArrayBufferContext } from '../../../mol-io/common/typed-array';

export async function readCcp4Header(file: FileHandle): Promise<{ header: Ccp4Header, littleEndian: boolean }> {
    const headerSize = 1024;
    const { buffer } = await file.readBuffer(0, headerSize);

    // 53  MAP         Character string 'MAP ' to identify file type
    const MAP = String.fromCharCode(
        buffer.readUInt8(52 * 4), buffer.readUInt8(52 * 4 + 1),
        buffer.readUInt8(52 * 4 + 2), buffer.readUInt8(52 * 4 + 3)
    );
    if (MAP !== 'MAP ') {
        throw new Error('ccp4 format error, missing "MAP " string');
    }

    // 54  MACHST      Machine stamp indicating machine type which wrote file
    //                 17 and 17 for big-endian or 68 and 65 for little-endian
    const MACHST = [ buffer.readUInt8(53 * 4), buffer.readUInt8(53 * 4 + 1) ];
    let littleEndian = false;
    if (MACHST[0] === 68 && MACHST[1] === 65) {
        littleEndian = true;
    } else if (MACHST[0] === 17 && MACHST[1] === 17) {
        littleEndian = false;
    } else {
        const modeLE = buffer.readInt32LE(3 * 4);
        if (modeLE <= 16) littleEndian = true;
    }

    const readInt = littleEndian ? (o: number) => buffer.readInt32LE(o * 4) : (o: number) => buffer.readInt32BE(o * 4);
    const readFloat = littleEndian ? (o: number) => buffer.readFloatLE(o * 4) : (o: number) => buffer.readFloatBE(o * 4);

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
    };

    return { header, littleEndian };
}

export async function readCcp4Slices(header: Ccp4Header, buffer: TypedArrayBufferContext, file: FileHandle, byteOffset: number, length: number, littleEndian: boolean) {
    if (isMapmode2to0(header)) {
        // data from mapmode2to0 is in MODE 0 (Int8) and needs to be scaled and written as float32
        const valueByteOffset = 3 * length;
        // read int8 data to last quarter of the read buffer
        await file.readBuffer(byteOffset, buffer.readBuffer, length, valueByteOffset);
        // get int8 view of last quarter of the read buffer
        const int8 = new Int8Array(buffer.valuesBuffer.buffer, valueByteOffset);
        // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
        const b1 = (header.AMAX - header.AMIN) / 255.0;
        const b0 = 0.5 * (header.AMIN + header.AMAX + b1);
        for (let j = 0, jl = length; j < jl; ++j) {
            buffer.values[j] = b1 * int8[j] + b0;
        }
    } else {
        await readTypedArray(buffer, file, byteOffset, length, 0, littleEndian);
    }
}

function getCcp4DataType(mode: number) {
    switch (mode) {
        case 0: return TypedArrayValueType.Int8;
        case 1: return TypedArrayValueType.Int16;
        case 2: return TypedArrayValueType.Float32;
        case 3: throw new Error('mode 3 unsupported, complex 16-bit integers');
        case 4: throw new Error('mode 4 unsupported, complex 32-bit reals');
        case 6: TypedArrayValueType.Uint16;
        case 16: throw new Error('mode 16 unsupported, unsigned char * 3 (for rgb data, non-standard)');
    }
    throw new Error(`unknown mode '${mode}'`);
}

/** check if the file was converted by mapmode2to0, see https://github.com/uglymol/uglymol */
function isMapmode2to0(header: Ccp4Header) {
    return header.userFlag1 === -128 && header.userFlag2 === 127;
}

export function getCcp4ValueType(header: Ccp4Header) {
    return isMapmode2to0(header) ? TypedArrayValueType.Float32 : getCcp4DataType(header.MODE);
}

export function getCcp4DataOffset(header: Ccp4Header) {
    return 256 * 4 + header.NSYMBT;
}

async function parseInternal(file: FileHandle, size: number, ctx: RuntimeContext): Promise<Ccp4File> {
    await ctx.update({ message: 'Parsing CCP4/MRC/MAP file...' });

    const { header, littleEndian } = await readCcp4Header(file);
    const offset = getCcp4DataOffset(header);
    const dataType = getCcp4DataType(header.MODE);
    const valueType = getCcp4ValueType(header);

    const count = header.NC * header.NR * header.NS;
    const elementByteSize = getElementByteSize(dataType);
    const byteCount = count * elementByteSize;

    const buffer = createTypedArrayBufferContext(count, valueType);
    readCcp4Slices(header, buffer, file, offset, byteCount, littleEndian);

    const result: Ccp4File = { header, values: buffer.values, name: file.name };
    return result;
}

export function parseFile(file: FileHandle, size: number) {
    return Task.create<Result<Ccp4File>>('Parse CCP4/MRC/MAP', async ctx => {
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