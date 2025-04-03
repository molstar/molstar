/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */


interface CustomString extends String {
    __string_like__: true,
    // toString()!!! // TODO
}

export type StringLike = string | String | CustomString

export function isStringLike(s: unknown): s is StringLike {
    return typeof s === 'string' || s instanceof String || (s as CustomString).__string_like__;
}

/** Maximum allow string length (might be bigger for some engines, but in Chrome and Node it is this). */
const MAX_STRING_LENGTH = 536_870_888;

/** Try to convert `StringLike` to a primitive `string`. Might fail if the contents is longer that max allowed string length. */
export function stringLikeToString(s: StringLike): string {
    try {
        return s.toString();
    } catch (err) {
        throw new Error(`Failed to convert StringLike object into string. This might be because the length ${s.length} exceeds maximum allowed string length ${MAX_STRING_LENGTH}. (${err})`);
    }
}
