/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

/** Convert a value to a JSON string on a single line (no newlines) */
export function onelinerJsonString(value: any): string {
    return JSON.stringify(value);
}

/** Convert a value to a canonical JSON string (deterministic ordering) */
export function canonicalJsonString(value: any): string {
    return JSON.stringify(value, (_key, val) => {
        if (val && typeof val === 'object' && !Array.isArray(val)) {
            // Sort object keys for canonical representation
            const sorted: any = {};
            Object.keys(val).sort().forEach(k => {
                sorted[k] = val[k];
            });
            return sorted;
        }
        return val;
    });
}

/** Type for JSON-serializable values */
export type Jsonable = string | number | boolean | null | undefined | Jsonable[] | { [key: string]: Jsonable };
