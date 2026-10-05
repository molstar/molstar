import { isPlainObject } from '@molstar/mvs-builder/object';

export type Jsonable = string | number | boolean | null | Jsonable[] | { [key: string]: Jsonable | undefined };

export function canonicalJsonString(obj: Jsonable) {
    return JSON.stringify(obj, (_key, value) => isPlainObject(value) ? sortObjectKeys(value) : value);
}

export function onelinerJsonString(obj: Jsonable) {
    return JSON.stringify(obj, undefined, '\t').replace(/,\n\t*/g, ', ').replace(/\n\t*/g, '');
}

function sortObjectKeys<T extends {}>(obj: T): T {
    const result = {} as T;
    for (const key of Object.keys(obj).sort() as (keyof T)[]) {
        const value = obj[key];
        if (value !== undefined) result[key] = value;
    }
    return result;
}
