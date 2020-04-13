/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReaderResult as Result } from '../result';
import { Task } from '../../../mol-task';
import { Cell } from '../../../mol-math/geometry/spacegroup/cell';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Mutable } from '../../../mol-util/type-helpers';
import { uint8ToString } from '../../common/binary';

export interface DcdHeader {
    readonly NSET: number,
    readonly ISTART: number,
    readonly NSAVC: number,
    readonly NAMNF: number,
    readonly DELTA: number,
    readonly TITLE: string,
    readonly NATOM: number
}

export interface DcdFrame {
    readonly cell: Cell
    readonly elementCount: number

    // positions
    readonly x: ArrayLike<number>
    readonly y: ArrayLike<number>
    readonly z: ArrayLike<number>
}

export interface DcdFile {
    readonly header: DcdHeader
    readonly frames: DcdFrame[]
}

export function _parseDcd(data: Uint8Array): DcdFile {
    // http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

    // The DCD format is structured as follows
    //   (FORTRAN UNFORMATTED, with Fortran data type descriptions):
    // HDR     NSET    ISTRT   NSAVC   5-ZEROS NATOM-NFREAT    DELTA   9-ZEROS
    // `CORD'  #files  step 1  step    zeroes  (zero)          timestep  (zeroes)
    //                         interval
    // C*4     INT     INT     INT     5INT    INT             DOUBLE  9INT
    // ==========================================================================
    // NTITLE          TITLE
    // INT (=2)        C*MAXTITL
    //                 (=32)
    // ==========================================================================
    // NATOM
    // #atoms
    // INT
    // ==========================================================================
    // X(I), I=1,NATOM         (DOUBLE)
    // Y(I), I=1,NATOM
    // Z(I), I=1,NATOM
    // ==========================================================================

    const dv = new DataView(data.buffer);

    const header: Mutable<DcdHeader> = Object.create(null);
    const frames: DcdFrame[] = [];

    let nextPos = 0;

    // header block

    const intView = new Int32Array(data.buffer, 0, 23);
    const ef = intView[0] !== dv.getInt32(0); // endianess flag
    // swap byte order when big endian (84 indicates little endian)
    if (intView[0] !== 84) {
        const n = data.byteLength;
        for (let i = 0; i < n; i += 4) {
            dv.setFloat32(i, dv.getFloat32(i), true);
        }
    }
    if (intView[0] !== 84) {
        throw new Error('dcd bad format, header block start');
    }

    // format indicator, should read 'CORD'
    const formatString = String.fromCharCode(
        dv.getUint8(4), dv.getUint8(5),
        dv.getUint8(6), dv.getUint8(7)
    );
    if (formatString !== 'CORD') {
        throw new Error('dcd bad format, format string');
    }
    let isCharmm = false;
    let extraBlock = false;
    let fourDims = false;
    // version field in charmm, unused in X-PLOR
    if (intView[22] !== 0) {
        isCharmm = true;
        if (intView[12] !== 0) extraBlock = true;
        if (intView[13] === 1) fourDims = true;
    }
    header.NSET = intView[2];
    header.ISTART = intView[3];
    header.NSAVC = intView[4];
    header.NAMNF = intView[10];

    if (isCharmm) {
        header.DELTA = dv.getFloat32(44, ef);
    } else {
        header.DELTA = dv.getFloat64(44, ef);
    }

    if (intView[22] !== 84) {
        throw new Error('dcd bad format, header block end');
    }
    nextPos = nextPos + 21 * 4 + 8;

    // title block

    const titleEnd = dv.getInt32(nextPos, ef);
    const titleStart = nextPos + 1;
    if ((titleEnd - 4) % 80 !== 0) {
        throw new Error('dcd bad format, title block start');
    }
    header.TITLE = uint8ToString(data.subarray(titleStart, titleEnd));
    if (dv.getInt32(titleStart + titleEnd + 4 - 1, ef) !== titleEnd) {
        throw new Error('dcd bad format, title block end');
    }

    nextPos = nextPos + titleEnd + 8;

    // natom block

    if (dv.getInt32(nextPos, ef) !== 4) {
        throw new Error('dcd bad format, natom block start');
    }
    header.NATOM = dv.getInt32(nextPos + 4, ef);
    if (dv.getInt32(nextPos + 8, ef) !== 4) {
        throw new Error('dcd bad format, natom block end');
    }
    nextPos = nextPos + 4 + 8;

    // fixed atoms block

    if (header.NAMNF > 0) {
        // TODO read coordinates and indices of fixed atoms
        throw new Error('dcd format with fixed atoms unsupported, aborting');
    }

    // frames

    const natom = header.NATOM;
    const natom4 = natom * 4;

    for (let i = 0, n = header.NSET; i < n; ++i) {
        const frame: Mutable<DcdFrame> = Object.create({
            elementCount: natom
        });

        if (extraBlock) {
            nextPos += 4; // block start
            // cell: A, alpha, B, beta, gamma, C (doubles)
            const size = Vec3.create(
                dv.getFloat64(nextPos, ef),
                dv.getFloat64(nextPos + 2 * 8, ef),
                dv.getFloat64(nextPos + 5 * 8, ef)
            );
            const anglesInRadians = Vec3.create(
                dv.getFloat64(nextPos + 1, ef),
                dv.getFloat64(nextPos + 3 * 8, ef),
                dv.getFloat64(nextPos + 4 * 8, ef)
            );
            frame.cell = Cell.create(size, anglesInRadians);
            nextPos += 48;
            nextPos += 4; // block end
        }

        // xyz coordinates
        for (let j = 0; j < 3; ++j) {
            if (dv.getInt32(nextPos, ef) !== natom4) {
                throw new Error(`dcd bad format, coord block start: ${i}, ${j}`);
            }
            nextPos += 4; // block start
            const c = new Float32Array(data.buffer, nextPos, natom);
            if (j === 0) frame.x = c;
            else if (j === 1) frame.y = c;
            else frame.z = c;

            nextPos += natom4;
            if (dv.getInt32(nextPos, ef) !== natom4) {
                throw new Error(`dcd bad format, coord block end: ${i}, ${j}`);
            }
            nextPos += 4; // block end
        }

        if (fourDims) {
            const bytes = dv.getInt32(nextPos, ef);
            nextPos += 4 + bytes + 4; // block start + skip + block end
        }

        frames.push(frame);
    }

    return { header, frames };
}

export function parseDcd(data: Uint8Array) {
    return Task.create<Result<DcdFile>>('Parse DCD', async ctx => {
        try {
            const dcdFile = _parseDcd(data);
            return Result.success(dcdFile);
        } catch (e) {
            return Result.error(e);
        }
    });
}