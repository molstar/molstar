/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import { FileHandle } from '../../mol-io/common/file-handle';
import { SimpleBuffer } from '../../mol-io/common/simple-buffer';
import { defaults, noop } from '../../mol-util';

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