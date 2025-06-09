/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorNames } from './color';

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
export type ColorName = keyof typeof ColorNames

export const ColorName = {
    /** Decide if a string is a valid named color string */
    is(str: any): str is ColorName {
        return str in ColorNames;
    },
};

/** Convert `colorString` (either X11 color name like 'magenta' or hex code like '#ff00ff') to numeric color value.
 * Return `undefined` if `colorString` cannot be converted. */
export function decodeColor(colorString: string | undefined | null): number | undefined {
    if (colorString === undefined || colorString === null) return undefined;
    let result: number | undefined;
    if (HexColor.is(colorString)) {
        if (colorString.length === 4) {
            // convert short form to full form (#f0f -> #ff00ff)
            colorString = `#${colorString[1]}${colorString[1]}${colorString[2]}${colorString[2]}${colorString[3]}${colorString[3]}`;
        }
        result = parseInt(colorString.slice(1), 16);
        if (!isNaN(result)) return result;
    }
    result = ColorNames[colorString.toLowerCase() as keyof typeof ColorNames];
    if (result !== undefined) return result;
    return undefined;
}

/** Create a simple hash for a given input string */
export function stringHash(input: string): string {
    let hash = 0;
    if (input.length === 0) return hash.toString(16);
    for (let i = 0; i < input.length; i++) {
        const char = input.charCodeAt(i);
        hash = ((hash << 5) - hash) + char;
        hash = hash & hash; // Convert to 32bit integer
    }
    return Math.abs(hash).toString(16);
}
