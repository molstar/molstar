/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '../../../mol-data/util';
import { StateObject } from '../../../mol-state';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { decodeColor as _decodeColor } from '../../../mol-util/color/utils';


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
    return values.some(isDefined);
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
export function decodeColor(colorString: string | number | undefined | null): Color | undefined {
    if (typeof colorString === 'number') {
        return Color(colorString);
    }
    return _decodeColor(colorString);
}

/** Regular expression matching a hexadecimal color string, e.g. '#FF1100' or '#f10' */
const hexColorRegex = /^#([0-9A-F]{3}){1,2}$/i;

/** Hexadecimal color string, e.g. '#FF1100' (the type matches more than just valid HexColor strings) */
export type HexColor = `#${string}`

export const HexColor = {
    /** Decide if a string is a valid hexadecimal color string (6-digit or 3-digit, e.g. '#FF1100' or '#f10') */
    is(str: any): str is HexColor {
        return typeof str === 'string' && hexColorRegex.test(str);
    },
};

/** Named color string, e.g. 'red' */
export type ColorName = keyof ColorNames

export const ColorName = {
    /** Decide if a string is a valid named color string */
    is(str: any): str is ColorName {
        return str in ColorNames;
    },
};

export function collectMVSReferences<T extends StateObject.Ctor>(type: T[], dependencies: Record<string, StateObject>): Record<string, StateObject.From<T>['data']> {
    const ret: any = {};

    for (const key of Object.keys(dependencies)) {
        const o = dependencies[key];
        let okType = false;
        for (const t of type) {
            if (t.is(o)) {
                okType = true;
                break;
            }
        }
        if (!okType || !o.tags) continue;
        for (const tag of o.tags) {
            if (tag.startsWith('mvs-ref:')) {
                ret[tag.substring(8)] = o.data;
                break;
            }
        }
    }

    return ret;
}

export function getMVSReferenceObject<T extends StateObject.Ctor>(type: T[], dependencies: Record<string, StateObject> | undefined, ref: string): StateObject | undefined {
    if (!dependencies) return undefined;

    for (const key of Object.keys(dependencies)) {
        const o = dependencies[key];
        let okType = false;
        for (const t of type) {
            if (t.is(o)) {
                okType = true;
                break;
            }
        }
        if (!okType || !o.tags) continue;
        for (const tag of o.tags) {
            if (tag.startsWith('mvs-ref:')) {
                if (tag.substring(8) === ref) return o;
            }
        }
    }
}