/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import * as Format from './format';
import * as DataFormat from '../common/data-format';
import { FileHandle } from '../../../mol-io/common/file-handle';
import { SimpleBuffer } from '../../../mol-io/common/simple-buffer';
import { TypedArrayValueArray, TypedArrayValueType } from '../../../mol-io/common/typed-array';

const FORMAT_VERSION = '1.0.0';

export interface Progress {
    current: number,
    max: number
}

export interface ValuesInfo {
    sum: number,
    sqSum: number,
    min: number,
    max: number
}

export interface BlockBuffer {
    values: TypedArrayValueArray[],
    buffers: SimpleBuffer[],
    slicesWritten: number
}

export interface DownsamplingBuffer {
    /** dimensions (sampleCount[1], sampleCount[0] / 2, 1), axis order (K, H, L) */
    downsampleH: TypedArrayValueArray,
    /** "Cyclic" (in the 1st dimensions) buffer with dimensions (5, sampleCount[0] / 2, sampleCount[1] / 2), axis order (L, H, K),  */
    downsampleHK: TypedArrayValueArray,

    slicesWritten: number,
    startSliceIndex: number
}

export interface Sampling {
    /** How many values along each axis are collapsed into 1 */
    rate: number,

    sampleCount: number[],

    /** One per channel, same indexing */
    blocks: BlockBuffer,
    valuesInfo: ValuesInfo[],
    downsampling?: DownsamplingBuffer[],

    /** Info about location in the output file, 0 offset is where the header ends */
    byteOffset: number,
    byteSize: number,
    /** where to write the next block relative to the byteoffset */
    writeByteOffset: number
}

/** Kernel used for downsampling */
export interface Kernel {
    /** The kernel size is curently fixed at 5 */
    size: 5,
    /** Compute new sample as c[0] * data[i - 2] + ... + c[4] * data[i + 2] */
    coefficients: number[],
    /** Precomputed coefficients.sum() */
    coefficientSum: number
}

export interface Context {
    /** Output file handle  */
    file: FileHandle,

    /** Periodic are x-ray density files that cover the entire grid and have [0,0,0] origin */
    isPeriodic: boolean,

    channels: Format.Context[],
    valueType: TypedArrayValueType,
    blockSize: number,
    /** Able to store channels.length * blockSize^3 values. */
    cubeBuffer: SimpleBuffer,
    /** All values are stored in little endian format which might not be the native endian of the system  */
    litteEndianCubeBuffer: SimpleBuffer,

    kernel: Kernel,
    sampling: Sampling[],
    dataByteOffset: number,
    totalByteSize: number,

    progress: Progress
}

export function createHeader(ctx: Context): DataFormat.Header {
    const header = ctx.channels[0].data.header;
    const grid = header.grid;

    function normalize(data: number[]) {
        return [data[0] / grid[0], data[1] / grid[1], data[2] / grid[2]];
    }

    return {
        formatVersion: FORMAT_VERSION,
        valueType: header.valueType,
        blockSize: ctx.blockSize,
        axisOrder: header.axisOrder,
        origin: normalize(header.origin),
        dimensions: normalize(header.extent),
        spacegroup: { number: header.spacegroupNumber, size: header.cellSize, angles: header.cellAngles, isPeriodic: ctx.isPeriodic },
        channels: ctx.channels.map(c => c.data.header.name),
        sampling: ctx.sampling.map(s => {
            const N = s.sampleCount[0] * s.sampleCount[1] * s.sampleCount[2];
            const valuesInfo = [];
            for (const { sum, sqSum, min, max } of s.valuesInfo) {
                const mean = sum / N;
                const sigma = Math.sqrt(Math.max(0, sqSum / N - mean * mean));
                valuesInfo.push({ mean, sigma, min, max });
            }
            return {
                byteOffset: s.byteOffset,
                rate: s.rate,
                valuesInfo,
                sampleCount: s.sampleCount,
            };
        })
    };
}

export function samplingBlockCount(sampling: Sampling, blockSize: number) {
    return sampling.sampleCount.map(c => Math.ceil(c / blockSize)).reduce((c, v) => c * v, 1);
}