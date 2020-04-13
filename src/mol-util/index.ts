/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import BitFlags from './bit-flags';
import StringBuilder from './string-builder';
import UUID from './uuid';
import Mask from './mask';
import { Progress } from '../mol-task';

export * from './value-cell';
export { BitFlags, StringBuilder, UUID, Mask };

export const noop = function () { };

export function round(n: number, d: number) {
    let f = Math.pow(10, d);
    return Math.round(f * n) / f;
}

export function arrayEqual<T>(arr1: T[], arr2: T[]) {
    const length = arr1.length;
    if (length !== arr2.length) return false;
    for (let i = 0; i < length; i++) {
        if (arr1[i] !== arr2[i]) {
            return false;
        }
    }
    return true;
}

const hasOwnProperty = Object.prototype.hasOwnProperty;

export function deepEqual(a: any, b: any) {
    // from https://github.com/epoberezkin/fast-deep-equal MIT
    if (a === b) return true;

    const arrA = Array.isArray(a);
    const arrB = Array.isArray(b);

    if (arrA && arrB) {
        if (a.length !== b.length) return false;
        for (let i = 0; i < a.length; i++) {
            if (!deepEqual(a[i], b[i])) return false;
        }
        return true;
    }

    if (arrA !== arrB) return false;

    if (a && b && typeof a === 'object' && typeof b === 'object') {
        const keys = Object.keys(a);
        if (keys.length !== Object.keys(b).length) return false;

        const dateA = a instanceof Date;
        const dateB = b instanceof Date;
        if (dateA && dateB) return a.getTime() === b.getTime();
        if (dateA !== dateB) return false;

        const regexpA = a instanceof RegExp;
        const regexpB = b instanceof RegExp;
        if (regexpA && regexpB) return a.toString() === b.toString();
        if (regexpA !== regexpB) return false;

        for (let i = 0; i < keys.length; i++) {
            if (!hasOwnProperty.call(b, keys[i])) return false;
        }

        for (let i = 0; i < keys.length; i++) {
            if (!deepEqual(a[keys[i]], b[keys[i]])) return false;
        }

        return true;
    }

    return false;
}

export function shallowEqual(a: any, b: any) {
    if (a === b) return true;
    const arrA = Array.isArray(a);
    const arrB = Array.isArray(b);
    if (arrA && arrB) return shallowEqualArrays(a, b);
    if (arrA !== arrB) return false;
    if (a && b && typeof a === 'object' && typeof b === 'object') {
        return shallowEqualObjects(a, b);
    }
    return false;
}

export function shallowEqualObjects(a: {}, b: {}) {
    if (a === b) return true;
    if (!a || !b) return false;
    const keys = Object.keys(a);
    if (Object.keys(b).length !== keys.length) return false;
    for (const k of keys) {
        if (!hasOwnProperty.call(a, k) || (a as any)[k] !== (b as any)[k]) return false;
    }
    return true;
}

export default function shallowEqualArrays(a: any[], b: any[]) {
    if (a === b) return true;
    if (!a || !b) return false;
    if (a.length !== b.length) return false;
    for (let i = 0, il = a.length; i < il; ++i) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}

/** Returns `value` if not `undefined`, otherwise returns `defaultValue` */
export function defaults<T>(value: T | undefined, defaultValue: T): T {
    return value !== undefined ? value : defaultValue;
}

export function extend<S, T, U>(object: S, source: T, guard?: U): S & T & U {
    let v: any;

    let s = <any>source;
    let o = <any>object;
    let g = <any>guard;
    for (let k of Object.keys(source)) {
        v = s[k];
        if (v !== void 0) o[k] = v;
        else if (guard) o[k] = g[k];
    }

    if (guard) {
        for (let k of Object.keys(guard)) {
            v = o[k];
            if (v === void 0) o[k] = g[k];
        }
    }

    return <any>object;
}

export function shallowClone<T>(o: T): T {
    return extend({}, o) as T;
}

function _assign<T>(target: T): T {
    for (let s = 1; s < arguments.length; s++) {
        let from = arguments[s];
        for (let key of Object.keys(from)) {
            if (hasOwnProperty.call(from, key)) {
                (target as any)[key] = from[key];
            }
        }
    }
    return target;
}

export declare function _assignType<T>(o: T, ...from: any[]): T;
export const assign: (<T>(o: T, ...from: any[]) => T) = (Object as any).assign || _assign;

function _shallowMerge1<T>(source: T, update: T) {
    let changed = false;
    for (let k of Object.keys(update)) {
        if (!hasOwnProperty.call(update, k)) continue;

        if ((update as any)[k] !== (source as any)[k]) {
            changed = true;
            break;
        }
    }

    if (!changed) return source;
    return assign(shallowClone(source), update);
}

function _shallowMerge<T>(source: T) {
    let ret = source;

    for (let s = 1; s < arguments.length; s++) {
        if (!arguments[s]) continue;
        ret = _shallowMerge1(source, arguments[s]);
        if (ret !== source) {
            for (let i = s + 1; i < arguments.length; i++) {
                ret = assign(ret, arguments[i]);
            }
            break;
        }
    }
    return ret;
}

export const merge: (<T>(source: T, ...rest: Partial<T>[]) => T) = _shallowMerge;

function padTime(n: number) { return (n < 10 ? '0' : '') + n; }
export function formatTime(d: Date) {
    const h = d.getHours(), m = d.getMinutes(), s = d.getSeconds();
    return `${h}:${padTime(m)}:${padTime(s)}`;
}

export function formatProgress(p: Progress) {
    const tp = p.root.progress;
    if (tp.isIndeterminate) return tp.message;
    const x = (100 * tp.current / tp.max).toFixed(2);
    return `${tp.message} ${x}%`;
}