/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

/** Check if value is a plain object (not array, null, etc.) */
export function isPlainObject(value: any): value is object {
    return value !== null && typeof value === 'object' && !Array.isArray(value);
}

/** Deep clone an object */
export function deepClone<T>(obj: T): T {
    if (obj === null || typeof obj !== 'object') return obj;
    if (obj instanceof Date) return new Date(obj.getTime()) as any;
    if (Array.isArray(obj)) return obj.map(item => deepClone(item)) as any;

    const cloned = {} as T;
    for (const key in obj) {
        if (obj.hasOwnProperty(key)) {
            cloned[key] = deepClone(obj[key]);
        }
    }
    return cloned;
}

/** Pick specific keys from an object */
export function pickObjectKeys<T extends object, K extends keyof T>(obj: T, keys: readonly K[]): Pick<T, K> {
    const result = {} as Pick<T, K>;
    for (const key of keys) {
        if (key in obj) {
            result[key] = obj[key];
        }
    }
    return result;
}

/** Omit specific keys from an object */
export function omitObjectKeys<T, K extends keyof T>(obj: T, keys: readonly K[]): Omit<T, K> {
    const result = {} as Omit<T, K>;
    const keysSet = new Set(keys);
    for (const key in obj) {
        if (!keysSet.has(key as any)) {
            (result as any)[key] = obj[key];
        }
    }
    return result;
}

/** Map over object values while preserving keys */
export function mapObjectMap<T, U>(obj: Record<string, T>, fn: (value: T, key: string) => U): Record<string, U> {
    const result: Record<string, U> = {};
    for (const key in obj) {
        if (obj.hasOwnProperty(key)) {
            result[key] = fn(obj[key], key);
        }
    }
    return result;
}
