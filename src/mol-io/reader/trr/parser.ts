/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from NGL.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';

export interface TrrFile {
    frames: { count: number, x: Float32Array, y: Float32Array, z: Float32Array }[],
    boxes: number[][],
    times: number[],
    timeOffset: number,
    deltaTime: number
}

async function parseInternal(data: Uint8Array) {
    // https://github.com/gromacs/gromacs/blob/master/src/gromacs/fileio/trrio.cpp

    const dv = new DataView(data.buffer);

    const f: TrrFile = {
        frames: [],
        boxes: [],
        times: [],
        timeOffset: 0,
        deltaTime: 0
    };
    const coordinates = f.frames;
    const boxes = f.boxes;
    const times = f.times;

    let offset = 0;

    while (true) {
        // const magicnum = dv.getInt32(offset)
        // const i1 = dv.getFloat32(offset + 4)
        offset += 8;

        const versionSize = dv.getInt32(offset);
        offset += 4;
        offset += versionSize;

        // const irSize = dv.getInt32(offset)
        // const eSize = dv.getInt32(offset + 4)
        const boxSize = dv.getInt32(offset + 8);
        const virSize = dv.getInt32(offset + 12);
        const presSize = dv.getInt32(offset + 16);
        // const topSize = dv.getInt32(offset + 20)
        // const symSize = dv.getInt32(offset + 24)
        const coordSize = dv.getInt32(offset + 28);
        const velocitySize = dv.getInt32(offset + 32);
        const forceSize = dv.getInt32(offset + 36);
        const natoms = dv.getInt32(offset + 40);
        // const step = dv.getInt32(offset + 44)
        // const nre = dv.getInt32(offset + 48)
        offset += 52;

        const floatSize = boxSize / 9;
        const natoms3 = natoms * 3;

        // let lambda
        if (floatSize === 8) {
            times.push(dv.getFloat64(offset));
            // lambda = dv.getFloat64(offset + 8)
        } else {
            times.push(dv.getFloat32(offset));
            // lambda = dv.getFloat32(offset + 4)
        }
        offset += 2 * floatSize;

        if (boxSize) {
            const box = new Float32Array(9);
            if (floatSize === 8) {
                for (let i = 0; i < 9; ++i) {
                    box[i] = dv.getFloat64(offset) * 10;
                    offset += 8;
                }
            } else {
                for (let i = 0; i < 9; ++i) {
                    box[i] = dv.getFloat32(offset) * 10;
                    offset += 4;
                }
            }
            boxes.push(box as unknown as number[]);
        }

        // ignore, unused
        offset += virSize;

        // ignore, unused
        offset += presSize;

        if (coordSize) {
            const x = new Float32Array(natoms);
            const y = new Float32Array(natoms);
            const z = new Float32Array(natoms);
            if (floatSize === 8) {
                for (let i = 0; i < natoms; ++i) {
                    x[i] = dv.getFloat64(offset) * 10;
                    y[i] = dv.getFloat64(offset + 8) * 10;
                    z[i] = dv.getFloat64(offset + 16) * 10;
                    offset += 24;
                }
            } else {
                const tmp = new Uint32Array(data.buffer, offset, natoms3);
                for (let i = 0; i < natoms3; ++i) {
                    const value = tmp[i];
                    tmp[i] = (
                        ((value & 0xFF) << 24) | ((value & 0xFF00) << 8) |
                        ((value >> 8) & 0xFF00) | ((value >> 24) & 0xFF)
                    );
                }
                const frameCoords = new Float32Array(data.buffer, offset, natoms3);
                for (let i = 0; i < natoms; ++i) {
                    x[i] = frameCoords[i * 3] * 10;
                    y[i] = frameCoords[i * 3 + 1] * 10;
                    z[i] = frameCoords[i * 3 + 2] * 10;
                    offset += 12;
                }
            }
            coordinates.push({ count: natoms, x, y, z });
        }

        // ignore, unused
        offset += velocitySize;

        // ignore, unused
        offset += forceSize;

        if (offset >= data.byteLength) break;
    }

    if (times.length >= 1) {
        f.timeOffset = times[0];
    }
    if (times.length >= 2) {
        f.deltaTime = times[1] - times[0];
    }

    return f;
}

export function parseTrr(data: Uint8Array) {
    return Task.create<Result<TrrFile>>('Parse TRR', async ctx => {
        try {
            ctx.update({ canAbort: true, message: 'Parsing trajectory...' });
            const file = await parseInternal(data);
            return Result.success(file);
        } catch (e) {
            return Result.error('' + e);
        }
    });
}