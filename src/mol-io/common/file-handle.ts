/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaults, noop } from 'mol-util';
import { SimpleBuffer } from './simple-buffer';
// only import 'fs' in node.js
const fs = typeof document === 'undefined' ? require('fs') as typeof import('fs') : void 0;

export interface FileHandle {
    /**
     * Asynchronously reads data, returning buffer and number of bytes read
     *
     * @param position The offset from the beginning of the file from which data should be read.
     * @param sizeOrBuffer The buffer the data will be read from.
     * @param length The number of bytes to read.
     * @param byteOffset The offset in the buffer at which to start writing.
     */
    readBuffer(position: number, sizeOrBuffer: SimpleBuffer | number, length?: number, byteOffset?: number): Promise<{ bytesRead: number, buffer: SimpleBuffer }>

    /**
     * Asynchronously writes buffer, returning the number of bytes written.
     *
     * @param position — The offset from the beginning of the file where this data should be written.
     * @param buffer - The buffer data to be written.
     * @param length — The number of bytes to write. If not supplied, defaults to buffer.length - offset.
     */
    writeBuffer(position: number, buffer: SimpleBuffer, length?: number): Promise<number>

    /**
     * Synchronously writes buffer, returning the number of bytes written.
     *
     * @param position — The offset from the beginning of the file where this data should be written.
     * @param buffer - The buffer data to be written.
     * @param length — The number of bytes to write. If not supplied, defaults to buffer.length - offset.
     */
    writeBufferSync(position: number, buffer: SimpleBuffer, length?: number): number

    /** Closes a file handle */
    close(): void
}

export namespace FileHandle {
    export function fromBuffer(buffer: SimpleBuffer): FileHandle {
        return {
            readBuffer: (position: number, sizeOrBuffer: SimpleBuffer | number, size?: number, byteOffset?: number) => {
                if (typeof sizeOrBuffer === 'number') {
                    const start = position
                    const end = Math.min(buffer.length, start + (defaults(size, sizeOrBuffer)))
                    return Promise.resolve({ bytesRead: end - start, buffer: SimpleBuffer.fromUint8Array(buffer.subarray(start, end)) })
                } else {
                    if (size === void 0) {
                        return Promise.reject('readBuffer: Specify size.');
                    }
                    const start = position
                    const end = Math.min(buffer.length, start + defaults(size, sizeOrBuffer.length))
                    sizeOrBuffer.set(buffer.subarray(start, end), byteOffset)
                    return Promise.resolve({ bytesRead: end - start, buffer: sizeOrBuffer })
                }
            },
            writeBuffer: (position: number, buffer: SimpleBuffer, length?: number) => {
                length = defaults(length, buffer.length)
                console.warn('FileHandle.writeBuffer not implemented')
                return Promise.resolve(0)
            },
            writeBufferSync: (position: number, buffer: SimpleBuffer, length?: number, ) => {
                length = defaults(length, buffer.length)
                console.warn('FileHandle.writeSync not implemented')
                return 0
            },
            close: noop
        }
    }

    export function fromDescriptor(file: number): FileHandle {
        if (fs === undefined) throw new Error('fs module not available')
        return {
            readBuffer: (position: number, sizeOrBuffer: SimpleBuffer | number, length?: number, byteOffset?: number) => {
                return new Promise((res, rej) => {
                    if (typeof sizeOrBuffer === 'number') {
                        let buff = new Buffer(new ArrayBuffer(sizeOrBuffer));
                        fs.read(file, buff, 0, sizeOrBuffer, position, (err, bytesRead, buffer) => {
                            if (err) {
                                rej(err);
                                return;
                            }
                            res({ bytesRead, buffer });
                        });
                    } else {
                        if (length === void 0) {
                            rej('readBuffer: Specify size.');
                            return;
                        }
                        fs.read(file, sizeOrBuffer, byteOffset ? +byteOffset : 0, length, position, (err, bytesRead, buffer) => {
                            if (err) {
                                rej(err);
                                return;
                            }
                            res({ bytesRead, buffer });
                        });
                    }
                })
            },
            writeBuffer: (position: number, buffer: Buffer, length?: number) => {
                return new Promise<number>((res, rej) => {
                    fs.write(file, buffer, 0, length !== void 0 ? length : buffer.length, position, (err, written) => {
                        if (err) rej(err);
                        else res(written);
                    })
                })
            },
            writeBufferSync: (position: number, buffer: Uint8Array, length?: number) => {
                return fs.writeSync(file, buffer, 0, length, position);
            },
            close: () => {
                try {
                    if (file !== void 0) fs.close(file, noop);
                } catch (e) {

                }
            }
        }
    }
}