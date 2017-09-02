/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * This ensures there is only 1 instance of a short string.
 * Also known as string interning, see https://en.wikipedia.org/wiki/String_interning
 */
export type ShortStringPool = { [key: string]: string }
export namespace ShortStringPool {
    export function create(): ShortStringPool { return Object.create(null); }
    export function get(pool: ShortStringPool, str: string) {
        if (str.length > 6) return str;
        const value = pool[str];
        if (value !== void 0) return value;
        pool[str] = str;
        return str;
    }
}
