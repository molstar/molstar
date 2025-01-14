/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import * as fs from 'fs';
import fetch from 'node-fetch';
import { FileHandle } from '../../mol-io/common/file-handle';
import { SimpleBuffer } from '../../mol-io/common/simple-buffer';
import { defaults, noop } from '../../mol-util';
import { openRead } from '../volume/common/file';
import { downloadGs } from './google-cloud-storage';


/** Create a file handle from either a file path or a URL (supports file://, http(s)://, gs:// protocols).  */
export async function fileHandleFromPathOrUrl(pathOrUrl: string, name: string): Promise<FileHandle> {
    if (pathOrUrl.startsWith('gs://')) {
        return fileHandleFromGS(pathOrUrl, name);
    } else if (pathOrUrl.startsWith('http://') || pathOrUrl.startsWith('https://')) {
        return fileHandleFromHTTP(pathOrUrl, name);
    } else if (pathOrUrl.startsWith('file://')) {
        return fileHandleFromDescriptor(await openRead(pathOrUrl.slice('file://'.length)), name);
    } else {
        return fileHandleFromDescriptor(await openRead(pathOrUrl), name);
    }
}

export function fileHandleFromDescriptor(file: number, name: string): FileHandle {
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

/** Create a read-only file handle from a Google Cloud Storage URL (gs://bucket-name/file-name).  */
export function fileHandleFromGS(url: string, name: string): FileHandle {
    return {
        name,
        readBuffer: async (position: number, sizeOrBuffer: SimpleBuffer | number, length?: number, byteOffset?: number) => {
            let outBuffer: SimpleBuffer;
            if (typeof sizeOrBuffer === 'number') {
                length = defaults(length, sizeOrBuffer);
                outBuffer = SimpleBuffer.fromArrayBuffer(new ArrayBuffer(sizeOrBuffer));
            } else {
                length = defaults(length, sizeOrBuffer.length);
                outBuffer = sizeOrBuffer;
            }
            let data: Buffer;
            try {
                data = await downloadGs(url, { decompress: false, start: position, end: position + length - 1 });
            } catch (err) {
                err.isFileNotFound = true;
                throw err;
            }
            const bytesRead = data.copy(outBuffer, byteOffset);
            if (length !== bytesRead) {
                console.warn(`byteCount ${length} and bytesRead ${bytesRead} differ`);
            }
            return { bytesRead, buffer: outBuffer };
        },
        writeBuffer: () => {
            throw new Error('Writing to Google Cloud Storage file handle not supported');
        },
        writeBufferSync: () => {
            throw new Error('Writing to Google Cloud Storage file handle not supported');
        },
        close: () => { },
    };
}

/** Create a read-only file handle from a HTTP or HTTPS URL.  */
export function fileHandleFromHTTP(url: string, name: string): FileHandle {
    let innerHandle: FileHandle | undefined = undefined;

    return {
        name,
        readBuffer: async (position: number, sizeOrBuffer: SimpleBuffer | number, length?: number, byteOffset?: number) => {
            if (!innerHandle) {
                const response = await fetch(url);
                if (response.ok) {
                    const buffer = await response.buffer();
                    innerHandle = FileHandle.fromBuffer(buffer, name);
                } else {
                    const error = new Error(`fileHandleFromHTTP: Fetch failed with status code ${response.status}`);
                    (error as any).isFileNotFound = true;
                    throw error;
                }
            }
            return innerHandle.readBuffer(position, sizeOrBuffer, length, byteOffset);
        },
        writeBuffer: () => {
            throw new Error('Writing to HTTP(S) file handle not supported');
        },
        writeBufferSync: () => {
            throw new Error('Writing to HTTP(S) file handle not supported');
        },
        close: () => { },
    };
}
