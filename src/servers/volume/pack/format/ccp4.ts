/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FileHandle } from '../../../../mol-io/common/file-handle';
import { readCcp4Header, readCcp4Slices, getCcp4DataOffset, getCcp4ValueType } from '../../../../mol-io/reader/ccp4/parser';
import { Header, Provider, Data } from '../format';
import { getCcp4Origin } from '../../../../mol-model-formats/volume/ccp4';
import { Ccp4Header } from '../../../../mol-io/reader/ccp4/schema';

async function readHeader(name: string, file: FileHandle) {
    const { header: ccp4Header, littleEndian } = await readCcp4Header(file);

    const header: Header = {
        name,
        valueType: getCcp4ValueType(ccp4Header),
        grid: [ccp4Header.NX, ccp4Header.NY, ccp4Header.NZ],
        axisOrder: [ccp4Header.MAPC, ccp4Header.MAPR, ccp4Header.MAPS].map(i => i - 1),
        extent: [ccp4Header.NC, ccp4Header.NR, ccp4Header.NS],
        origin: getCcp4Origin(ccp4Header),
        spacegroupNumber: ccp4Header.ISPG,
        cellSize: [ccp4Header.xLength, ccp4Header.yLength, ccp4Header.zLength],
        cellAngles: [ccp4Header.alpha, ccp4Header.beta, ccp4Header.gamma],
        littleEndian,
        dataOffset: getCcp4DataOffset(ccp4Header),
        originalHeader: ccp4Header
    };
    // "normalize" the grid axis order
    header.grid = [header.grid[header.axisOrder[0]], header.grid[header.axisOrder[1]], header.grid[header.axisOrder[2]]];
    return header;
}

export async function readSlices(data: Data) {
    const { slices, header } = data;
    if (slices.isFinished) {
        return;
    }

    const { extent, originalHeader } = header;
    const sliceSize = extent[0] * extent[1];
    const sliceByteOffset = slices.buffer.elementByteSize * sliceSize * slices.slicesRead;
    const sliceCount = Math.min(slices.sliceCapacity, extent[2] - slices.slicesRead);
    const sliceByteCount = slices.buffer.elementByteSize * sliceCount * sliceSize;

    await readCcp4Slices(originalHeader as Ccp4Header, slices.buffer, data.file, header.dataOffset + sliceByteOffset, sliceByteCount, header.littleEndian);
    slices.slicesRead += sliceCount;
    slices.sliceCount = sliceCount;

    if (slices.slicesRead >= extent[2]) {
        slices.isFinished = true;
    }
}

export const Ccp4Provider: Provider = { readHeader, readSlices };