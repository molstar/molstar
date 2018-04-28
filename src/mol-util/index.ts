/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import BitFlags from './bit-flags'
import StringBuilder from './string-builder'
import UUID from './uuid'
import Mask from './mask'

export * from './value-cell'
export { BitFlags, StringBuilder, UUID, Mask }

export function arrayEqual<T>(arr1: T[], arr2: T[]) {
    const length = arr1.length
    if (length !== arr2.length) return false
    for (let i = 0; i < length; i++) {
        if (arr1[i] !== arr2[i]) {
            return false
        }
    }
    return true
}

export function deepEqual(a: any, b: any) {
    // from https://github.com/epoberezkin/fast-deep-equal MIT
    if (a === b) return true;

    const arrA = Array.isArray(a)
    const arrB = Array.isArray(b)

    if (arrA && arrB) {
        if (a.length !== b.length) return false
        for (let i = 0; i < a.length; i++) {
            if (!deepEqual(a[i], b[i])) return false
        }
        return true
    }

    if (arrA !== arrB) return false

    if (a && b && typeof a === 'object' && typeof b === 'object') {
        const keys = Object.keys(a)
        if (keys.length !== Object.keys(b).length) return false;

        const dateA = a instanceof Date
        const dateB = b instanceof Date
        if (dateA && dateB) return a.getTime() === b.getTime()
        if (dateA !== dateB) return false

        const regexpA = a instanceof RegExp
        const regexpB = b instanceof RegExp
        if (regexpA && regexpB) return a.toString() === b.toString()
        if (regexpA !== regexpB) return false

        for (let i = 0; i < keys.length; i++) {
            if (!Object.prototype.hasOwnProperty.call(b, keys[i])) return false
        }

        for (let i = 0; i < keys.length; i++) {
            if (!deepEqual(a[keys[i]], b[keys[i]])) return false
        }

        return true
    }

    return false
}

export function defaults(value: any, defaultValue: any) {
    return value !== undefined ? value : defaultValue
}