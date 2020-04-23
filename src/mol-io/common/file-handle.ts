/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaults, noop } from '../../mol-util';
import { SimpleBuffer } from './simple-buffer';
// only import 'fs' in node.js
const fs = typeof document === 'undefined' ? require('fs') as typeof import('fs') : void 0;

export interface FileHandle {
    name: string
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
     * @param length — The number of bytes to write. If not supplied, defaults to buffer.length
     */
    writeBuffer(position: number, buffer: SimpleBuffer, length?: number): Promise<number>

    /**
     * Synchronously writes buffer, returning the number of bytes written.
     *
     * @param position — The offset from the beginning of the file where this data should be written.
     * @param buffer - The buffer data to be written.
     * @param length — The number of bytes to write. If not supplied, defaults to buffer.length
     */
    writeBufferSync(position: number, buffer: SimpleBuffer, length?: number): number

    /** Closes a file handle */
    close(): void
}

export namespace FileHandle {
    export function fromBuffer(buffer: SimpleBuffer, name: string): FileHandle {
        return {
            name,
            readBuffer: (position: number, sizeOrBuffer: SimpleBuffer | number, size?: number, byteOffset?: number) => {
                let bytesRead: number;
                let outBuffer: SimpleBuffer;
                if (typeof sizeOrBuffer === 'number') {
                    size = defaults(size, sizeOrBuffer);
                    const start = position;
                    const end = Math.min(buffer.length, start + size);
                    bytesRead = end - start;
                    outBuffer = SimpleBuffer.fromUint8Array(new Uint8Array(buffer.buffer, start, end - start));
                } else {
                    size = defaults(size, sizeOrBuffer.length);
                    const start = position;
                    const end = Math.min(buffer.length, start + size);
                    sizeOrBuffer.set(buffer.subarray(start, end), byteOffset);
                    bytesRead = end - start;
                    outBuffer = sizeOrBuffer;
                }
                if (size !== bytesRead) {
                    console.warn(`byteCount ${size} and bytesRead ${bytesRead} differ`);
                }
                return Promise.resolve({ bytesRead, buffer: outBuffer });
            },
            writeBuffer: (position: number, buffer: SimpleBuffer, length?: number) => {
                length = defaults(length, buffer.length);
                console.error('.writeBuffer not implemented for FileHandle.fromBuffer');
                return Promise.resolve(0);
            },
            writeBufferSync: (position: number, buffer: SimpleBuffer, length?: number, ) => {
                length = defaults(length, buffer.length);
                console.error('.writeSync not implemented for FileHandle.fromBuffer');
                return 0;
            },
            close: noop
        };
    }

    export function fromDescriptor(file: number, name: string): FileHandle {
        if (fs === undefined) throw new Error('fs module not available');
        return {
            name,
            readBuffer: (position: number, sizeOrBuffer: SimpleBuffer | number, length?: number, byteOffset?: number) => {
                return new Promise((res, rej) => {
                    let outBuffer: SimpleBuffer;
                    if (typeof sizeOrBuffer === 'number') {
                        byteOffset = defaults(byteOffset, 0);
                        length = defaults(length, sizeOrBuffer);
                        outBuffer = SimpleBuffer.fromArrayBuffer(new ArrayBuffer(sizeOrBuffer));
                    } else {
                        byteOffset = defaults(byteOffset, 0);
                        length = defaults(length, sizeOrBuffer.length);
                        outBuffer = sizeOrBuffer;
                    }
                    fs.read(file, outBuffer, byteOffset, length, position, (err, bytesRead, buffer) => {
                        if (err) {
                            rej(err);
                            return;
                        }
                        if (length !== bytesRead) {
                            console.warn(`byteCount ${length} and bytesRead ${bytesRead} differ`);
                        }
                        res({ bytesRead, buffer });
                    });
                });
            },
            writeBuffer: (position: number, buffer: SimpleBuffer, length?: number) => {
                length = defaults(length, buffer.length);
                return new Promise<number>((res, rej) => {
                    fs.write(file, buffer, 0, length, position, (err, written) => {
                        if (err) rej(err);
                        else res(written);
                    });
                });
            },
            writeBufferSync: (position: number, buffer: Uint8Array, length?: number) => {
                length = defaults(length, buffer.length);
                return fs.writeSync(file, buffer, 0, length, position);
            },
            close: () => {
                try {
                    if (file !== void 0) fs.close(file, noop);
                } catch (e) {

                }
            }
        };
    }
}