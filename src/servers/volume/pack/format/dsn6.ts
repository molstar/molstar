/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { FileHandle } from 'mol-io/common/file-handle';
import { Header, Provider, Mode } from '../format';
import { readSlices } from './common';
import { readDsn6Header, dsn6HeaderSize } from 'mol-io/reader/dsn6/parser';

async function readHeader(name: string, file: FileHandle) {
    const { header: dsn6Header, littleEndian } = await readDsn6Header(file)

    const header: Header = {
        name,
        mode: Mode.Int16,
        grid: [dsn6Header.xRate, dsn6Header.yRate, dsn6Header.zRate],
        axisOrder: [0, 1, 2],
        extent: [dsn6Header.xExtent, dsn6Header.yExtent, dsn6Header.zExtent],
        origin: [dsn6Header.xStart, dsn6Header.yStart, dsn6Header.zStart],
        spacegroupNumber: 0, // P 1
        cellSize: [dsn6Header.xlen, dsn6Header.ylen, dsn6Header.zlen],
        cellAngles: [dsn6Header.alpha, dsn6Header.beta, dsn6Header.gamma],
        littleEndian,
        dataOffset: dsn6HeaderSize
    };
    return header;
}

export const Dsn6Provider: Provider = { readHeader, readSlices }