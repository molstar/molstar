/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from NGL.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext, Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';

export interface XtcFile {
    frames: { count: number, x: Float32Array, y: Float32Array, z: Float32Array }[],
    boxes: number[][],
    times: number[],
    timeOffset: number,
    deltaTime: number
}

const MagicInts = new Uint32Array([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
    16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031,
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216
]);
const FirstIdx = 9;
// const LastIdx = MagicInts.length

namespace Decoder {
    export function sizeOfInt(size: number) {
        let num = 1;
        let numOfBits = 0;
        while (size >= num && numOfBits < 32) {
            numOfBits++;
            num <<= 1;
        }
        return numOfBits;
    }

    const _tmpBytes = new Uint8Array(32);

    export function sizeOfInts(numOfInts: number, sizes: number[]) {
        let numOfBytes = 1;
        let numOfBits = 0;
        _tmpBytes[0] = 1;
        for (let i = 0; i < numOfInts; i++) {
            let bytecnt;
            let tmp = 0;
            for (bytecnt = 0; bytecnt < numOfBytes; bytecnt++) {
                tmp += _tmpBytes[bytecnt] * sizes[i];
                _tmpBytes[bytecnt] = tmp & 0xff;
                tmp >>= 8;
            }
            while (tmp !== 0) {
                _tmpBytes[bytecnt++] = tmp & 0xff;
                tmp >>= 8;
            }
            numOfBytes = bytecnt;
        }
        let num = 1;
        numOfBytes--;
        while (_tmpBytes[numOfBytes] >= num) {
            numOfBits++;
            num *= 2;
        }
        return numOfBits + numOfBytes * 8;
    }

    const _buffer = new ArrayBuffer(8 * 3);
    export const buf = new Int32Array(_buffer);
    const uint32view = new Uint32Array(_buffer);

    export function decodeBits(cbuf: Uint8Array, offset: number, numOfBits1: number) {
        let numOfBits = numOfBits1;
        const mask = (1 << numOfBits) - 1;
        let lastBB0 = uint32view[1];
        let lastBB1 = uint32view[2];
        let cnt = buf[0];
        let num = 0;

        while (numOfBits >= 8) {
            lastBB1 = (lastBB1 << 8) | cbuf[offset + cnt++];
            num |= (lastBB1 >> lastBB0) << (numOfBits - 8);
            numOfBits -= 8;
        }

        if (numOfBits > 0) {
            if (lastBB0 < numOfBits) {
                lastBB0 += 8;
                lastBB1 = (lastBB1 << 8) | cbuf[offset + cnt++];
            }
            lastBB0 -= numOfBits;
            num |= (lastBB1 >> lastBB0) & ((1 << numOfBits) - 1);
        }

        num &= mask;
        buf[0] = cnt;
        buf[1] = lastBB0;
        buf[2] = lastBB1;

        return num;
    }

    function decodeByte(cbuf: Uint8Array, offset: number) {
        // special version of decodeBits with numOfBits = 8

        // const mask = 0xff; // (1 << 8) - 1;
        // let lastBB0 = uint32view[1];
        let lastBB1 = uint32view[2];
        let cnt = buf[0];

        lastBB1 = (lastBB1 << 8) | cbuf[offset + cnt];

        buf[0] = cnt + 1;
        // buf[1] = lastBB0;
        buf[2] = lastBB1;

        return (lastBB1 >> uint32view[1]) & 0xff;
    }

    const intBytes = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    // new Int32Array(32);

    export function decodeInts(cbuf: Uint8Array, offset: number, numOfBits1: number, sizes: number[], nums: number[]) {
        let numOfBits = numOfBits1;
        let numOfBytes = 0;

        intBytes[0] = 0;
        intBytes[1] = 0;
        intBytes[2] = 0;
        intBytes[3] = 0;

        while (numOfBits > 8) {
            // this is inversed??? why??? because of the endiannness???
            intBytes[numOfBytes++] = decodeByte(cbuf, offset);
            numOfBits -= 8;
        }

        if (numOfBits > 0) {
            intBytes[numOfBytes++] = decodeBits(cbuf, offset, numOfBits);
        }

        for (let i = 2; i > 0; i--) {
            let num = 0;
            const s = sizes[i];
            for (let j = numOfBytes - 1; j >= 0; j--) {
                num = (num << 8) | intBytes[j];
                const t = (num / s) | 0;
                intBytes[j] = t;
                num = num - t * s;
            }
            nums[i] = num;
        }
        nums[0] = intBytes[0] | (intBytes[1] << 8) | (intBytes[2] << 16) | (intBytes[3] << 24);
    }
}

function undefinedError() {
    throw new Error('(xdrfile error) Undefined error.');
}

async function parseInternal(ctx: RuntimeContext, data: Uint8Array) {
    // https://github.com/gromacs/gromacs/blob/master/src/gromacs/fileio/xtcio.cpp
    // https://github.com/gromacs/gromacs/blob/master/src/gromacs/fileio/libxdrf.cpp

    const dv = new DataView(data.buffer, data.byteOffset);

    const f: XtcFile = {
        frames: [],
        boxes: [],
        times: [],
        timeOffset: 0,
        deltaTime: 0
    };
    const coordinates = f.frames;
    const boxes = f.boxes;
    const times = f.times;

    const minMaxInt = [0, 0, 0, 0, 0, 0];
    const sizeint = [0, 0, 0];
    const bitsizeint = [0, 0, 0];
    const sizesmall = [0, 0, 0];
    const thiscoord = [0.1, 0.1, 0.1];
    const prevcoord = [0.1, 0.1, 0.1];

    let offset = 0;
    const buf = Decoder.buf;
    // const buf2 = new Uint32Array(buf.buffer);

    while (true) {
        let frameCoords: XtcFile['frames'][0];

        // const magicnum = dv.getInt32(offset)
        const natoms = dv.getInt32(offset + 4);
        // const step = dv.getInt32(offset + 8)
        offset += 12;

        times.push(dv.getFloat32(offset));
        offset += 4;

        const box = new Float32Array(9);
        for (let i = 0; i < 9; ++i) {
            box[i] = dv.getFloat32(offset) * 10;
            offset += 4;
        }
        boxes.push(box as unknown as number[]);

        if (natoms <= 9) { // no compression
            frameCoords = { count: natoms / 3, x: new Float32Array(natoms / 3), y: new Float32Array(natoms / 3), z: new Float32Array(natoms / 3) };
            for (let i = 0; i < natoms / 3; ++i) {
                frameCoords.x[i] = dv.getFloat32(offset);
                frameCoords.y[i] = dv.getFloat32(offset);
                frameCoords.z[i] = dv.getFloat32(offset);
                offset += 4;
            }
        } else {
            buf[0] = buf[1] = buf[2] = 0;
            sizeint[0] = sizeint[1] = sizeint[2] = 0;
            sizesmall[0] = sizesmall[1] = sizesmall[2] = 0;
            bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
            thiscoord[0] = thiscoord[1] = thiscoord[2] = 0;
            prevcoord[0] = prevcoord[1] = prevcoord[2] = 0;

            frameCoords = { count: natoms, x: new Float32Array(natoms), y: new Float32Array(natoms), z: new Float32Array(natoms) };
            let lfp = 0;

            const lsize = dv.getInt32(offset);
            offset += 4;
            const precision = dv.getFloat32(offset);
            offset += 4;

            minMaxInt[0] = dv.getInt32(offset);
            minMaxInt[1] = dv.getInt32(offset + 4);
            minMaxInt[2] = dv.getInt32(offset + 8);
            minMaxInt[3] = dv.getInt32(offset + 12);
            minMaxInt[4] = dv.getInt32(offset + 16);
            minMaxInt[5] = dv.getInt32(offset + 20);
            sizeint[0] = minMaxInt[3] - minMaxInt[0] + 1;
            sizeint[1] = minMaxInt[4] - minMaxInt[1] + 1;
            sizeint[2] = minMaxInt[5] - minMaxInt[2] + 1;
            offset += 24;

            let bitsize;
            if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
                bitsizeint[0] = Decoder.sizeOfInt(sizeint[0]);
                bitsizeint[1] = Decoder.sizeOfInt(sizeint[1]);
                bitsizeint[2] = Decoder.sizeOfInt(sizeint[2]);
                bitsize = 0; // flag the use of large sizes
            } else {
                bitsize = Decoder.sizeOfInts(3, sizeint);
            }

            let smallidx = dv.getInt32(offset);
            offset += 4;

            let tmpIdx = smallidx - 1;
            tmpIdx = (FirstIdx > tmpIdx) ? FirstIdx : tmpIdx;
            let smaller = (MagicInts[tmpIdx] / 2) | 0;
            let smallnum = (MagicInts[smallidx] / 2) | 0;

            sizesmall[0] = sizesmall[1] = sizesmall[2] = MagicInts[smallidx];

            let adz = Math.ceil(dv.getInt32(offset) / 4) * 4;
            offset += 4;

            const invPrecision = 1.0 / precision;
            let run = 0;
            let i = 0;

            // const buf8 = new Uint8Array(data.buffer, data.byteOffset + offset, 32 * 4); // 229...

            thiscoord[0] = thiscoord[1] = thiscoord[2] = 0;

            while (i < lsize) {
                if (bitsize === 0) {
                    thiscoord[0] = Decoder.decodeBits(data, offset, bitsizeint[0]);
                    thiscoord[1] = Decoder.decodeBits(data, offset, bitsizeint[1]);
                    thiscoord[2] = Decoder.decodeBits(data, offset, bitsizeint[2]);
                } else {
                    Decoder.decodeInts(data, offset, bitsize, sizeint, thiscoord);
                }

                i++;

                thiscoord[0] += minMaxInt[0];
                thiscoord[1] += minMaxInt[1];
                thiscoord[2] += minMaxInt[2];

                prevcoord[0] = thiscoord[0];
                prevcoord[1] = thiscoord[1];
                prevcoord[2] = thiscoord[2];

                const flag = Decoder.decodeBits(data, offset, 1);
                let isSmaller = 0;

                if (flag === 1) {
                    run = Decoder.decodeBits(data, offset, 5);
                    isSmaller = run % 3;
                    run -= isSmaller;
                    isSmaller--;
                }

                // if ((lfp-ptrstart)+run > size3){
                //   fprintf(stderr, "(xdrfile error) Buffer overrun during decompression.\n");
                //   return 0;
                // }

                if (run > 0) {
                    thiscoord[0] = thiscoord[1] = thiscoord[2] = 0;

                    for (let k = 0; k < run; k += 3) {
                        Decoder.decodeInts(data, offset, smallidx, sizesmall, thiscoord);
                        i++;

                        thiscoord[0] += prevcoord[0] - smallnum;
                        thiscoord[1] += prevcoord[1] - smallnum;
                        thiscoord[2] += prevcoord[2] - smallnum;

                        if (k === 0) {
                            // interchange first with second atom for
                            // better compression of water molecules
                            let tmpSwap = thiscoord[0];
                            thiscoord[0] = prevcoord[0];
                            prevcoord[0] = tmpSwap;

                            tmpSwap = thiscoord[1];
                            thiscoord[1] = prevcoord[1];
                            prevcoord[1] = tmpSwap;

                            tmpSwap = thiscoord[2];
                            thiscoord[2] = prevcoord[2];
                            prevcoord[2] = tmpSwap;

                            frameCoords.x[lfp] = prevcoord[0] * invPrecision;
                            frameCoords.y[lfp] = prevcoord[1] * invPrecision;
                            frameCoords.z[lfp] = prevcoord[2] * invPrecision;
                            lfp++;
                        } else {
                            prevcoord[0] = thiscoord[0];
                            prevcoord[1] = thiscoord[1];
                            prevcoord[2] = thiscoord[2];
                        }
                        frameCoords.x[lfp] = thiscoord[0] * invPrecision;
                        frameCoords.y[lfp] = thiscoord[1] * invPrecision;
                        frameCoords.z[lfp] = thiscoord[2] * invPrecision;
                        lfp++;
                    }
                } else {
                    frameCoords.x[lfp] = thiscoord[0] * invPrecision;
                    frameCoords.y[lfp] = thiscoord[1] * invPrecision;
                    frameCoords.z[lfp] = thiscoord[2] * invPrecision;
                    lfp++;
                }

                smallidx += isSmaller;

                if (isSmaller < 0) {
                    smallnum = smaller;
                    if (smallidx > FirstIdx) {
                        smaller = (MagicInts[smallidx - 1] / 2) | 0;
                    } else {
                        smaller = 0;
                    }
                } else if (isSmaller > 0) {
                    smaller = smallnum;
                    smallnum = (MagicInts[smallidx] / 2) | 0;
                }
                sizesmall[0] = sizesmall[1] = sizesmall[2] = MagicInts[smallidx];

                if (sizesmall[0] === 0 || sizesmall[1] === 0 || sizesmall[2] === 0) {
                    undefinedError();
                }
            }
            offset += adz;
        }

        for (let c = 0; c < natoms; c++) {
            frameCoords.x[c] *= 10;
            frameCoords.y[c] *= 10;
            frameCoords.z[c] *= 10;
        }

        coordinates.push(frameCoords);

        if (ctx.shouldUpdate) {
            await ctx.update({ current: offset, max: data.length });
        }

        if (offset >= data.length) break;
    }

    if (times.length >= 1) {
        f.timeOffset = times[0];
    }
    if (times.length >= 2) {
        f.deltaTime = times[1] - times[0];
    }

    return f;
}

export function parseXtc(data: Uint8Array) {
    return Task.create<Result<XtcFile>>('Parse XTC', async ctx => {
        try {
            ctx.update({ canAbort: true, message: 'Parsing trajectory...' });
            const file = await parseInternal(ctx, data);
            return Result.success(file);
        } catch (e) {
            return Result.error('' + e);
        }
    });
}