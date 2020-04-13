/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Schema from './binary-schema';
import { FileHandle } from '../../../mol-io/common/file-handle';
import { TypedArrayValueType } from '../../../mol-io/common/typed-array';

export interface Spacegroup {
    number: number,
    size: number[],
    angles: number[],
    /** Determine if the data should be treated as periodic or not. (e.g. X-ray = periodic, EM = not periodic) */
    isPeriodic: boolean,
}

export interface ValuesInfo {
    mean: number,
    sigma: number,
    min: number,
    max: number
}

export interface Sampling {
    byteOffset: number,

    /** How many values along each axis were collapsed into 1 */
    rate: number,
    valuesInfo: ValuesInfo[],

    /** Number of samples along each axis, in axisOrder */
    sampleCount: number[]
}

export interface Header {
    /** Format version number  */
    formatVersion: string,

    /** Axis order from the slowest to fastest moving, same as in CCP4 */
    axisOrder: number[],

    /** Origin in fractional coordinates, in axisOrder */
    origin: number[],

    /** Dimensions in fractional coordinates, in axisOrder */
    dimensions: number[],

    spacegroup: Spacegroup,
    channels: string[],

    /** Determines the data type of the values */
    valueType: TypedArrayValueType,

    /** The value are stored in blockSize^3 cubes */
    blockSize: number,
    sampling: Sampling[]
}

namespace _schema {
    const { array, obj, int, bool, float, str } = Schema;

    export const schema = obj<Header>([
        ['formatVersion', str],
        ['axisOrder', array(int)],
        ['origin', array(float)],
        ['dimensions', array(float)],
        ['spacegroup', obj<Spacegroup>([
            ['number', int],
            ['size', array(float)],
            ['angles', array(float)],
            ['isPeriodic', bool],
        ])],
        ['channels', array(str)],
        ['valueType', str],
        ['blockSize', int],
        ['sampling', array(obj<Sampling>([
            ['byteOffset', float],
            ['rate', int],
            ['valuesInfo', array(obj<ValuesInfo>([
                ['mean', float],
                ['sigma', float],
                ['min', float],
                ['max', float]
            ]))],
            ['sampleCount', array(int)]
        ]))]
    ]);
}

const headerSchema = _schema.schema;

export function encodeHeader(header: Header) {
    return Schema.encode(headerSchema, header);
}

export async function readHeader(file: FileHandle): Promise<{ header: Header, dataOffset: number }> {
    let { buffer } = await file.readBuffer(0, 4 * 4096);
    const headerSize = buffer.readInt32LE(0);

    if (headerSize > buffer.byteLength - 4) {
        buffer = (await file.readBuffer(0, headerSize + 4)).buffer;
    }

    const header = Schema.decode<Header>(headerSchema, buffer, 4);
    return { header, dataOffset: headerSize + 4 };
}