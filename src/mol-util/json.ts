/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { isPlainObject } from './object';


/** A JSON-serializable value */
export type Jsonable = string | number | boolean | null | Jsonable[] | { [key: string]: Jsonable | undefined }

/** Return a canonical string representation for a JSON-able object,
 * independent from object key order and undefined properties. */
export function canonicalJsonString(obj: Jsonable) {
    return JSON.stringify(obj, (key, value) => isPlainObject(value) ? sortObjectKeys(value) : value);
}

/** Return a pretty JSON representation for a JSON-able object,
 * (single line, but use space after comma). E.g. '{"name": "Bob", "favorite_numbers": [1, 2, 3]}' */
export function onelinerJsonString(obj: Jsonable) {
    return JSON.stringify(obj, undefined, '\t').replace(/,\n\t*/g, ', ').replace(/\n\t*/g, '');
}

/** Return a copy of object `obj` with alphabetically sorted keys and dropped keys whose value is undefined. */
function sortObjectKeys<T extends {}>(obj: T): T {
    const result = {} as T;
    for (const key of Object.keys(obj).sort() as (keyof T)[]) {
        const value = obj[key];
        if (value !== undefined) {
            result[key] = value;
        }
    }
    return result;
}
