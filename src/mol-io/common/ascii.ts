/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

export function asciiWrite(data: Uint8Array, str: string) {
    for (let i = 0, il = str.length; i < il; ++i) {
        data[i] = str.charCodeAt(i);
    }
}