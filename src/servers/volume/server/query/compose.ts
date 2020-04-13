/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data-model';
import * as Box from '../algebra/box';
import * as Coords from '../algebra/coordinate';
import { createTypedArrayBufferContext, getElementByteSize, readTypedArray } from '../../../../mol-io/common/typed-array';

export default async function compose(query: Data.QueryContext.Data) {
    for (const block of query.samplingInfo.blocks) {
        await fillBlock(query, block);
    }
}

async function readBlock(query: Data.QueryContext.Data, coord: Coords.Grid<'Block'>, blockBox: Box.Fractional): Promise<Data.BlockData> {
    const { valueType, blockSize } = query.data.header;
    const elementByteSize = getElementByteSize(valueType);
    const numChannels = query.data.header.channels.length;
    const blockSampleCount = Box.dimensions(Box.fractionalToGrid(blockBox, query.samplingInfo.sampling.dataDomain));
    const size = numChannels * blockSampleCount[0] * blockSampleCount[1] * blockSampleCount[2];
    const byteSize = elementByteSize * size;
    const dataSampleCount = query.data.header.sampling[query.samplingInfo.sampling.index].sampleCount;
    const buffer = createTypedArrayBufferContext(size, valueType);
    const byteOffset = query.samplingInfo.sampling.byteOffset
        + elementByteSize * numChannels * blockSize
        * (blockSampleCount[1] * blockSampleCount[2] * coord[0]
            + dataSampleCount[0] * blockSampleCount[2] * coord[1]
            + dataSampleCount[0] * dataSampleCount[1] * coord[2]);

    const values = await readTypedArray(buffer, query.data.file, byteOffset, byteSize, 0);
    return {
        sampleCount: blockSampleCount,
        values
    };
}

function fillData(query: Data.QueryContext.Data, blockData: Data.BlockData, blockGridBox: Box.Grid<'BlockGrid'>, queryGridBox: Box.Grid<'Query'>) {
    const source = blockData.values;

    const { sizeX: tSizeH, sizeXY: tSizeHK } = Coords.gridMetrics(query.samplingInfo.gridDomain.sampleCount);
    const { sizeX: sSizeH, sizeXY: sSizeHK } = Coords.gridMetrics(blockData.sampleCount);

    const offsetTarget = queryGridBox.a[0] + queryGridBox.a[1] * tSizeH + queryGridBox.a[2] * tSizeHK;

    const [maxH, maxK, maxL] = Box.dimensions(blockGridBox);

    for (let channelIndex = 0, _ii = query.data.header.channels.length; channelIndex < _ii; channelIndex++) {
        const target = query.values[channelIndex];
        const offsetSource = channelIndex * blockGridBox.a.domain.sampleVolume
            + blockGridBox.a[0] + blockGridBox.a[1] * sSizeH + blockGridBox.a[2] * sSizeHK;

        for (let l = 0; l < maxL; l++) {
            for (let k = 0; k < maxK; k++) {
                for (let h = 0; h < maxH; h++) {
                    target[offsetTarget + h + k * tSizeH + l * tSizeHK]
                        = source[offsetSource + h + k * sSizeH + l * sSizeHK];
                }
            }
        }
    }
}

function createBlockGridDomain(block: Coords.Grid<'Block'>, grid: Coords.GridDomain<'Data'>): Coords.GridDomain<'BlockGrid'> {
    const blockBox = Box.fractionalFromBlock(block);
    const origin = blockBox.a;
    const dimensions = Coords.sub(blockBox.b, blockBox.a);
    const sampleCount = Coords.sampleCounts(dimensions, grid.delta);
    return Coords.domain<'BlockGrid'>('BlockGrid', { origin, dimensions, delta: grid.delta, sampleCount });
}

/** Read the block data and fill all the overlaps with the query region. */
async function fillBlock(query: Data.QueryContext.Data, block: Data.QueryBlock) {
    const baseBox = Box.fractionalFromBlock(block.coord);
    const blockGridDomain = createBlockGridDomain(block.coord, query.samplingInfo.sampling.dataDomain);

    const blockData: Data.BlockData = await readBlock(query, block.coord, baseBox);

    for (const offset of block.offsets) {
        const offsetQueryBox = Box.shift(query.samplingInfo.fractionalBox, offset);
        const dataBox = Box.intersect(baseBox, offsetQueryBox);
        if (!dataBox) continue;

        const offsetDataBox = Box.shift(dataBox, Coords.invert(offset));

        const blockGridBox = Box.clampGridToSamples(Box.fractionalToGrid(dataBox, blockGridDomain));
        const queryGridBox = Box.clampGridToSamples(Box.fractionalToGrid(offsetDataBox, query.samplingInfo.gridDomain));

        fillData(query, blockData, blockGridBox, queryGridBox);
    }
}