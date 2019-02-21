/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Data } from '../format';
import { readTypedArray } from 'mol-io/common/typed-array';

export async function readSlices(data: Data) {
    const { slices, header } = data;
    if (slices.isFinished) {
        return;
    }

    const { extent } = header;
    const sliceSize = extent[0] * extent[1];
    const sliceByteOffset = slices.buffer.elementByteSize * sliceSize * slices.slicesRead;
    const sliceCount = Math.min(slices.sliceCapacity, extent[2] - slices.slicesRead);
    const sliceByteCount = sliceCount * sliceSize;
    console.log('sliceByteOffset', sliceByteOffset, 'sliceSize', sliceSize, 'sliceCount', sliceCount)

    await readTypedArray(slices.buffer, data.file, header.dataOffset + sliceByteOffset, sliceByteCount, 0, header.littleEndian);
    slices.slicesRead += sliceCount;
    slices.sliceCount = sliceCount;

    if (slices.slicesRead >= extent[2]) {
        slices.isFinished = true;
    }
}