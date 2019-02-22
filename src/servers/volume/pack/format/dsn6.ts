/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { FileHandle } from 'mol-io/common/file-handle';
import { Header, Provider, Data } from '../format';
import { readDsn6Header, dsn6HeaderSize, parseDsn6Values } from 'mol-io/reader/dsn6/parser';
import { TypedArrayValueType } from 'mol-io/common/typed-array';
import { Dsn6Header } from 'mol-io/reader/dsn6/schema';

async function readHeader(name: string, file: FileHandle) {
    const { header: dsn6Header, littleEndian } = await readDsn6Header(file)

    const header: Header = {
        name,
        valueType: TypedArrayValueType.Int16,
        grid: [dsn6Header.xRate, dsn6Header.yRate, dsn6Header.zRate],
        axisOrder: [0, 1, 2],
        extent: [dsn6Header.xExtent, dsn6Header.yExtent, dsn6Header.zExtent],
        origin: [dsn6Header.xStart, dsn6Header.yStart, dsn6Header.zStart],
        spacegroupNumber: 1, // P 1
        cellSize: [dsn6Header.xlen, dsn6Header.ylen, dsn6Header.zlen],
        cellAngles: [dsn6Header.alpha, dsn6Header.beta, dsn6Header.gamma],
        littleEndian,
        dataOffset: dsn6HeaderSize,
        originalHeader: dsn6Header
    };
    return header;
}

export async function readSlices(data: Data) {
    // TODO due to the dsn6 data layout, the file must be read a a whole, need check if the file is too big for that

    const { slices, header, file } = data;
    if (slices.isFinished) {
        return;
    }

    const { extent, originalHeader } = header;
    const sliceCount = extent[2]

    const { xExtent, yExtent, zExtent } = originalHeader as Dsn6Header
    const xBlocks = Math.ceil(xExtent / 8)
    const yBlocks = Math.ceil(yExtent / 8)
    const zBlocks = Math.ceil(zExtent / 8)
    const count = xBlocks * 8 * yBlocks * 8 * zBlocks * 8
    const elementByteSize = 1
    const byteCount = count * elementByteSize

    const { buffer } = await file.readBuffer(dsn6HeaderSize, byteCount)
    await parseDsn6Values(originalHeader as Dsn6Header, buffer, slices.values as Float32Array) // TODO fix cast

    slices.slicesRead += sliceCount;
    slices.sliceCount = sliceCount;

    if (slices.slicesRead >= extent[2]) {
        slices.isFinished = true;
    }
}

export const Dsn6Provider: Provider = { readHeader, readSlices }