/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

const hasOwnProperty = Object.prototype.hasOwnProperty;

/** Assign to the object if a given property in update is undefined */
export function assignIfUndefined<T extends {}>(to: Partial<T>, full: T): T {
    for (const k of Object.keys(full)) {
        if (!hasOwnProperty.call(full, k)) continue;

        if (typeof (to as any)[k] === 'undefined') {
            (to as any)[k] = (full as any)[k];
        }
    }
    return to as T;
}

/** Create new object if any property in "update" changes in "source". */
export function shallowMerge2<T extends {}>(source: T, update: Partial<T>): T {
    // Adapted from LiteMol (https://github.com/dsehnal/LiteMol)
    let changed = false;
    for (const k of Object.keys(update)) {
        if (!hasOwnProperty.call(update, k)) continue;

        if ((update as any)[k] !== (source as any)[k]) {
            changed = true;
            break;
        }
    }

    if (!changed) return source;
    return Object.assign({}, source, update);
}

export function shallowEqual<T extends {}>(a: T, b: T) {
    if (!a) {
        if (!b) return true;
        return false;
    }
    if (!b) return false;

    const keys = Object.keys(a);
    if (Object.keys(b).length !== keys.length) return false;
    for (const k of keys) {
        if (!hasOwnProperty.call(a, k) || (a as any)[k] !== (b as any)[k]) return false;
    }

    return true;
}

export function shallowMerge<T extends {}>(source: T, ...rest: (Partial<T> | undefined)[]): T {
    return shallowMergeArray(source, rest);
}

export function shallowMergeArray<T extends {}>(source: T, rest: (Partial<T> | undefined)[]): T {
    // Adapted from LiteMol (https://github.com/dsehnal/LiteMol)
    let ret: any = source;

    for (let s = 0; s < rest.length; s++) {
        if (!rest[s]) continue;
        ret = shallowMerge2(source, rest[s] as T);
        if (ret !== source) {
            for (let i = s + 1; i < rest.length; i++) {
                ret = Object.assign(ret, rest[i]);
            }
            break;
        }
    }
    return ret;
}

/** Simple deep clone for number, boolean, string, null, undefined, object, array */
export function deepClone<T>(source: T): T {
    if (null === source || 'object' !== typeof source) return source;

    if (source instanceof Array) {
        const copy: any[] = [];
        for (let i = 0, len = source.length; i < len; i++) {
            copy[i] = deepClone(source[i]);
        }
        return copy as any as T;
    }

    // `instanceof Object` does not find `Object.create(null)`
    if (typeof source === 'object' && !('prototype' in source)) {
        const copy: { [k: string]: any } = {};
        for (const k in source) {
            if (hasOwnProperty.call(source, k)) copy[k] = deepClone(source[k]);
        }
        return copy as any as T;
    }

    throw new Error(`Can't clone, type "${typeof source}" unsupported`);
}

/** Return a new object with the same keys, where function `f` is applied to each value.
 * Equivalent to Pythonic `{k: f(v) for k, v in obj.items()}` */
export function mapObjectMap<T, S>(obj: { [k: string]: T }, f: (v: T) => S): { [k: string]: S } {
    const ret: any = {};
    for (const k of Object.keys(obj)) {
        ret[k] = f((obj as any)[k]);
    }
    return ret;
}

/** Return an object with keys being the elements of `array` and values computed by `getValue` function.
 * Equivalent to Pythonic `{k: getValue(k) for k in array}` */
export function mapArrayToObject<K extends keyof any, V>(array: readonly K[], getValue: (key: K) => V): Record<K, V> {
    const result = {} as Record<K, V>;
    for (const key of array) {
        result[key] = getValue(key);
    }
    return result;
}

export function objectForEach<T extends {}, V extends T[K], K extends keyof T & string>(o: T, f: (v: V, k: K) => void) {
    if (!o) return;
    for (const k of Object.keys(o)) {
        f((o as any)[k], k as K);
    }
}

/** Return an object with keys `keys` and their values same as in `obj`, i.e. `{ key: obj[key] for key in keys }` */
export function pickObjectKeys<T extends {}, K extends keyof T>(obj: T, keys: readonly K[]): Pick<T, K> {
    const result: Partial<Pick<T, K>> = {};
    for (const key of keys) {
        if (Object.hasOwn(obj, key)) {
            result[key] = obj[key];
        }
    }
    return result as Pick<T, K>;
}

/** Same as `pickObjectKeys` but allows loading values into a different key or skipping keys, i.e. return `{ key: obj[remapping[key] ?? key] for key in keys if remapping[key] !== null }`  */
export function pickObjectKeysWithRemapping<V>(obj: Record<string, V>, keys: string[], remapping: Record<string, string | null>): Record<string, V> {
    const result: Record<string, V> = {};
    for (const key of keys) {
        let srcKey = remapping[key];
        if (srcKey === null) continue; // Ignore key
        if (srcKey === undefined) srcKey = key; // Implicit key mapping
        if (Object.hasOwn(obj, srcKey)) {
            result[key] = obj[srcKey];
        }
    }
    return result;
}

export function objectOfArraysToArrayOfObjects<T extends object>(objectOfArrays: { [key in keyof T]: T[key][] }): T[] {
    let n: number | undefined = undefined;
    for (const key in objectOfArrays) {
        const length = objectOfArrays[key].length;
        if (n !== undefined && n !== length) throw new Error('FormatError: arrays must have the same length.');
        n = length;
    }
    if (n === undefined) return []; // empty input object

    const out: T[] = new Array(n);
    for (let i = 0; i < n; i++) {
        const item: Partial<T> = {};
        for (const key in objectOfArrays) {
            item[key] = objectOfArrays[key][i];
        }
        out[i] = item as T;
    }
    return out;
}

/** Return an object same as `obj` but without keys `keys` */
export function omitObjectKeys<T extends {}, K extends keyof T>(obj: T, omitKeys: readonly K[]): Omit<T, K> {
    const result: T = { ...obj };
    for (const key of omitKeys) {
        delete result[key];
    }
    return result as Omit<T, K>;
}

/** Create an object from keys and values (first key maps to first value etc.) */
export function objectFromKeysAndValues<K extends keyof any, V>(keys: K[], values: V[]): Record<K, V> {
    const obj: Partial<Record<K, V>> = {};
    for (let i = 0; i < keys.length; i++) {
        obj[keys[i]] = values[i];
    }
    return obj as Record<K, V>;
}

/** Decide if `obj` is a good old object (not array or null or other type). */
export function isPlainObject(obj: any): boolean {
    return typeof obj === 'object' && obj !== null && !Array.isArray(obj);
}

/** Like `Promise.all` but with objects instead of arrays */
export async function promiseAllObj<T extends {}>(promisesObj: { [key in keyof T]: Promise<T[key]> }): Promise<T> {
    const keys = Object.keys(promisesObj);
    const promises = Object.values(promisesObj);
    const results = await Promise.all(promises);
    return objectFromKeysAndValues(keys, results) as any;
}
