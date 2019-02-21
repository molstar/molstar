/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FileHandle } from 'mol-io/common/file-handle';
import { readCcp4Header } from 'mol-io/reader/ccp4/parser';
import { Header, Provider } from '../format';
import { readSlices } from './common';
import { TypedArrayValueType } from 'mol-io/common/typed-array';

function getTypedArrayValueType(mode: number) {
    switch (mode) {
        case 2: return TypedArrayValueType.Float32
        case 1: return TypedArrayValueType.Int16
        case 0: return TypedArrayValueType.Int8
    }
    throw new Error(`ccp4 mode '${mode}' unsupported`);
}

async function readHeader(name: string, file: FileHandle) {
    const { header: ccp4Header, littleEndian } = await readCcp4Header(file)

    const origin2k = [ccp4Header.originX, ccp4Header.originY, ccp4Header.originZ];
    const nxyzStart = [ccp4Header.NCSTART, ccp4Header.NRSTART, ccp4Header.NSSTART];
    const header: Header = {
        name,
        valueType: getTypedArrayValueType(ccp4Header.MODE),
        grid: [ccp4Header.NX, ccp4Header.NY, ccp4Header.NZ],
        axisOrder: [ccp4Header.MAPC, ccp4Header.MAPR, ccp4Header.MAPS].map(i => i - 1),
        extent: [ccp4Header.NC, ccp4Header.NR, ccp4Header.NS],
        origin: origin2k[0] === 0.0 && origin2k[1] === 0.0 && origin2k[2] === 0.0 ? nxyzStart : origin2k,
        spacegroupNumber: ccp4Header.ISPG,
        cellSize: [ccp4Header.xLength, ccp4Header.yLength, ccp4Header.zLength],
        cellAngles: [ccp4Header.alpha, ccp4Header.beta, ccp4Header.gamma],
        littleEndian,
        dataOffset: 256 * 4 + ccp4Header.NSYMBT /* symBytes */
    };
    // "normalize" the grid axis order
    header.grid = [header.grid[header.axisOrder[0]], header.grid[header.axisOrder[1]], header.grid[header.axisOrder[2]]];
    return header;
}

export const Ccp4Provider: Provider = { readHeader, readSlices }