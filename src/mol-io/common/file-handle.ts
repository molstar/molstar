/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { defaults } from 'mol-util';

export interface FileHandle {
    /** The number of bytes in the file */
    length: number
    /**
     * @param position The offset from the beginning of the file from which data should be read.
     * @param sizeOrBuffer The buffer the data will be written to. If a number a buffer of that size will be created.
     * @param size The number of bytes to read.
     * @param byteOffset The offset in the buffer at which to start writing.
     */
    readBuffer(position: number, sizeOrBuffer: Uint8Array | number, size?: number, byteOffset?: number): Promise<{ bytesRead: number, buffer: Uint8Array }>
}

export namespace FileHandle {
    export function fromBuffer(buffer: Uint8Array): FileHandle {
        return {
            length: buffer.length,
            readBuffer: (position: number, sizeOrBuffer: Uint8Array | number, size?: number, byteOffset?: number) => {
                if (typeof sizeOrBuffer === 'number') {
                    const start = position
                    const end = Math.min(buffer.length, start + (defaults(size, sizeOrBuffer)))
                    return Promise.resolve({
                        bytesRead: end - start,
                        buffer: buffer.subarray(start, end),
                    })
                } else {
                    if (size === void 0) {
                        return Promise.reject('readBuffer: Specify size.');
                    }
                    const start = position
                    const end = Math.min(buffer.length, start + defaults(size, sizeOrBuffer.length))
                    sizeOrBuffer.set(buffer.subarray(start, end), byteOffset)
                    return Promise.resolve({
                        bytesRead: end - start,
                        buffer: sizeOrBuffer,
                    })
                }
            }
        }
    }
}