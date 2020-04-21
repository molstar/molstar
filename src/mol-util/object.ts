/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const hasOwnProperty = Object.prototype.hasOwnProperty;

/** Assign to the object if a given property in update is undefined */
export function assignIfUndefined<T>(to: Partial<T>, full: T): T {
    for (const k of Object.keys(full)) {
        if (!hasOwnProperty.call(full, k)) continue;

        if (typeof (to as any)[k] === 'undefined') {
            (to as any)[k] = (full as any)[k];
        }
    }
    return to as T;
}

/** Create new object if any property in "update" changes in "source". */
export function shallowMerge2<T>(source: T, update: Partial<T>): T {
    // Adapted from LiteMol (https://github.com/dsehnal/LiteMol)
    let changed = false;
    for (let k of Object.keys(update)) {
        if (!hasOwnProperty.call(update, k)) continue;

        if ((update as any)[k] !== (source as any)[k]) {
            changed = true;
            break;
        }
    }

    if (!changed) return source;
    return Object.assign({}, source, update);
}

export function shallowEqual<T>(a: T, b: T) {
    if (!a) {
        if (!b) return true;
        return false;
    }
    if (!b) return false;

    let keys = Object.keys(a);
    if (Object.keys(b).length !== keys.length) return false;
    for (let k of keys) {
        if (!hasOwnProperty.call(a, k) || (a as any)[k] !== (b as any)[k]) return false;
    }

    return true;
}

export function shallowMerge<T>(source: T, ...rest: (Partial<T> | undefined)[]): T {
    return shallowMergeArray(source, rest);
}

export function shallowMergeArray<T>(source: T, rest: (Partial<T> | undefined)[]): T {
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
        for (let k in source) {
            if (hasOwnProperty.call(source, k)) copy[k] = deepClone(source[k]);
        }
        return copy as any as T;
    }

    throw new Error(`Can't clone, type "${typeof source}" unsupported`);
}

export function mapObjectMap<T, S>(o: { [k: string]: T }, f: (v: T) => S): { [k: string]: S } {
    const ret: any = { };
    for (const k of Object.keys(o)) {
        ret[k] = f((o as any)[k]);
    }
    return ret;
}

export function objectForEach<T>(o: { [k: string]: T }, f: (v: T, k: string) => void) {
    if (!o) return;
    for (const k of Object.keys(o)) {
        f((o as any)[k], k);
    }
}