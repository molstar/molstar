/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { FileHandle } from '../../../../mol-io/common/file-handle';
import { Header, Provider, Data } from '../format';
import { readDsn6Header, dsn6HeaderSize, parseDsn6Values, getDsn6Counts } from '../../../../mol-io/reader/dsn6/parser';
import { TypedArrayValueType } from '../../../../mol-io/common/typed-array';
import { Dsn6Header } from '../../../../mol-io/reader/dsn6/schema';

async function readHeader(name: string, file: FileHandle) {
    const { header: dsn6Header, littleEndian } = await readDsn6Header(file);

    const header: Header = {
        name,
        valueType: TypedArrayValueType.Float32,
        grid: [dsn6Header.xRate, dsn6Header.yRate, dsn6Header.zRate].reverse(),
        axisOrder: [0, 1, 2].reverse(),
        extent: [dsn6Header.xExtent, dsn6Header.yExtent, dsn6Header.zExtent].reverse(),
        origin: [dsn6Header.xStart, dsn6Header.yStart, dsn6Header.zStart].reverse(),
        spacegroupNumber: 1, // set as P 1, since it is not available in DSN6 files
        cellSize: [dsn6Header.xlen, dsn6Header.ylen, dsn6Header.zlen],
        cellAngles: [dsn6Header.alpha, dsn6Header.beta, dsn6Header.gamma],
        littleEndian,
        dataOffset: dsn6HeaderSize,
        originalHeader: dsn6Header
    };
    return header;
}

export async function readSlices(data: Data) {
    // TODO due to the dsn6 data layout we the read file into one big buffer
    //      to avoid this, either change the sampling algoritm to work with this layout or
    //      read the data into a collection of buffers that can be access like one big buffer
    //      => for now not worth putting time in, for big files better use another file format

    const { slices, header, file } = data;
    if (slices.isFinished) {
        return;
    }

    const { extent, dataOffset, originalHeader } = header;
    const sliceCount = extent[2];

    const { byteCount } = getDsn6Counts(originalHeader as Dsn6Header);
    if (byteCount > slices.maxBlockBytes) {
        throw new Error(`dsn6 file to large, can't read ${byteCount} bytes at once, increase block size or use another file format`);
    }

    const { buffer } = await file.readBuffer(dataOffset, byteCount);
    if (!(slices.values instanceof Float32Array)) {
        throw new Error(`dsn6 reader only supports Float32Array for output values`);
    }
    await parseDsn6Values(originalHeader as Dsn6Header, buffer, slices.values, header.littleEndian);

    slices.slicesRead += sliceCount;
    slices.sliceCount = sliceCount;

    if (slices.slicesRead >= extent[2]) {
        slices.isFinished = true;
    }
}

export const Dsn6Provider: Provider = { readHeader, readSlices };