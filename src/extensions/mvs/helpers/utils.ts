/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { hashString } from '../../../mol-data/util';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';


/** Represents either the result or the reason of failure of an operation that might have failed */
export type Maybe<T> = { ok: true, value: T } | { ok: false, error: any }

/** Try to await a promise and return an object with its result (if resolved) or with the error (if rejected) */
export async function safePromise<T>(promise: T): Promise<Maybe<Awaited<T>>> {
    try {
        const value = await promise;
        return { ok: true, value };
    } catch (error) {
        return { ok: false, error };
    }
}


/** A map where values are arrays. Handles missing keys when adding values. */
export class MultiMap<K, V> implements Mapping<K, V[]> {
    private _map = new Map();

    /** Return the array of values assidned to a key (or `undefined` if no such values) */
    get(key: K): V[] | undefined {
        return this._map.get(key);
    }
    /** Append value to a key (handles missing keys) */
    add(key: K, value: V) {
        if (!this._map.has(key)) {
            this._map.set(key, []);
        }
        this._map.get(key)!.push(value);
    }
}

/** Basic subset of `Map<K, V>`, only needs to have `get` method */
export type Mapping<K, V> = Pick<Map<K, V>, 'get'>

/** Implementation of `Map` where keys are integers
 * and most keys are expected to be from interval `[0, limit)`.
 * For the keys within this interval, performance is better than `Map` (implemented by array).
 * For the keys out of this interval, performance is slightly worse than `Map`. */
export class NumberMap<K extends number, V> implements Mapping<K, V> {
    private array: V[];
    private map: Map<K, V>;
    constructor(public readonly limit: K) {
        this.array = new Array(limit);
        this.map = new Map();
    }
    get(key: K): V | undefined {
        if (0 <= key && key < this.limit) return this.array[key];
        else return this.map.get(key);
    }
    set(key: K, value: V): void {
        if (0 <= key && key < this.limit) this.array[key] = value;
        else this.map.set(key, value);
    }
}


/** Return `true` if `value` is not `undefined` or `null`.
 * Prefer this over `value !== undefined`
 * (for maybe if we want to allow `null` in `AnnotationRow` in the future) */
export function isDefined<T>(value: T | undefined | null): value is T {
    return value !== undefined && value !== null;
}
/** Return `true` if at least one of `values` is not `undefined` or `null`. */
export function isAnyDefined(...values: any[]): boolean {
    return values.some(v => isDefined(v));
}
/** Return filtered array containing all original elements except `undefined` or `null`. */
export function filterDefined<T>(elements: (T | undefined | null)[]): T[] {
    return elements.filter(x => x !== undefined && x !== null) as T[];
}

/** Create an 8-hex-character hash for a given input string, e.g. 'spanish inquisition' -> '7f9ac4be' */
function stringHash32(input: string): string {
    const uint32hash = hashString(input) >>> 0; // >>>0 converts to uint32, LOL
    return uint32hash.toString(16).padStart(8, '0');
}
/** Create an 16-hex-character hash for a given input string, e.g. 'spanish inquisition' -> '7f9ac4be544330be'*/
export function stringHash(input: string): string {
    const reversed = input.split('').reverse().join('');
    return stringHash32(input) + stringHash32(reversed);
}

/** Return type of elements in a set */
export type ElementOfSet<S> = S extends Set<infer T> ? T : never


/** Convert `colorString` (either X11 color name like 'magenta' or hex code like '#ff00ff') to Color.
 * Return `undefined` if `colorString` cannot be converted. */
export function decodeColor(colorString: string | undefined): Color | undefined {
    if (colorString === undefined) return undefined;
    let result: Color | undefined;
    if (isHexColorString(colorString)) {
        if (colorString.length === 4) {
            // convert short form to full form (#f0f -> #ff00ff)
            colorString = `#${colorString[1]}${colorString[1]}${colorString[2]}${colorString[2]}${colorString[3]}${colorString[3]}`;
        }
        result = Color.fromHexStyle(colorString);
        if (result !== undefined && !isNaN(result)) return result;
    }
    result = ColorNames[colorString.toLowerCase() as keyof typeof ColorNames];
    if (result !== undefined) return result;
    return undefined;
}

/** Hexadecimal color string, e.g. '#FF1100' */
export type HexColor = string & { '@type': 'HexColorString' }
export function HexColor(str: string) {
    if (!isHexColorString(str)) {
        throw new Error(`ValueError: "${str}" is not a valid hex color string`);
    }
    return str as HexColor;
}

/** Regular expression matching a hexadecimal color string, e.g. '#FF1100' or '#f10' */
const hexColorRegex = /^#([0-9A-F]{3}){1,2}$/i;

/** Decide if a string is a valid hexadecimal color string (6-digit or 3-digit, e.g. '#FF1100' or '#f10') */
export function isHexColorString(str: any): str is HexColor {
    return typeof str === 'string' && hexColorRegex.test(str);
}
