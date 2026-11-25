/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

/** Whether `value` is a record made by an object literal (its prototype is `Object.prototype` or `null`). */
export function isPlainRecord(value: unknown): value is Record<string, unknown> {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === Object.prototype || prototype === null;
}

/**
 * The same own enumerable string keys with the same values (`Object.is`), counted with `for … in`: no key arrays.
 * Symbol keys are not compared.
 */
export function shallowEqualObjects(a: {}, b: {}): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  let count = 0;
  for (const key in a) {
    if (!hasOwnProperty.call(a, key)) continue;
    if (!hasOwnProperty.call(b, key) || !Object.is((a as any)[key], (b as any)[key])) return false;
    count++;
  }
  for (const key in b) {
    if (hasOwnProperty.call(b, key)) count--;
  }
  return count === 0;
}

/** The same length and items (`Object.is`). */
export function shallowEqualArrays(a: readonly any[], b: readonly any[]): boolean {
  if (a === b) return true;
  if (!a || !b || a.length !== b.length) return false;
  for (let i = 0, il = a.length; i < il; ++i) {
    if (!Object.is(a[i], b[i])) return false;
  }
  return true;
}

/**
 * Two objects with the same own keys and values (`Object.is`; two arrays item by item); two nullish values (`null`,
 * `undefined`) are equal.
 */
export function shallowEqual<T extends {}>(a: T, b: T): boolean {
  if (!a || !b) return !a && !b;
  if (Array.isArray(a) && Array.isArray(b)) return shallowEqualArrays(a, b);
  return shallowEqualObjects(a, b);
}

/**
 * The equality that decides "changed" for values held in states: `Object.is`; two arrays item by item; two plain
 * records (`isPlainRecord`) key by key; anything else (a class instance, a `Map`, a `Date`, a function) only to itself.
 */
export function sameShallow(a: unknown, b: unknown): boolean {
  if (Object.is(a, b)) return true;
  if (Array.isArray(a)) return Array.isArray(b) && shallowEqualArrays(a, b);
  return isPlainRecord(a) && isPlainRecord(b) && shallowEqualObjects(a, b);
}

/**
 * Every level of plain arrays and records (`isPlainRecord`) compared as `sameShallow` compares one; anything else by
 * identity (`Object.is`). For values built fresh from data (a list of items, a label's lines) that must not count as a
 * change when equal.
 */
export function sameData(a: unknown, b: unknown): boolean {
  if (Object.is(a, b)) return true;
  if (Array.isArray(a)) {
    if (!Array.isArray(b) || a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!sameData(a[i], b[i])) return false;
    return true;
  }
  if (!isPlainRecord(a) || !isPlainRecord(b)) return false;
  let count = 0;
  for (const key in a) {
    if (!hasOwnProperty.call(a, key)) continue;
    if (!hasOwnProperty.call(b, key) || !sameData(a[key], b[key])) return false;
    count++;
  }
  for (const key in b) {
    if (hasOwnProperty.call(b, key)) count--;
  }
  return count === 0;
}

/** Structural equality of plain data (arrays, objects, `Date`, `RegExp`); adapted from fast-deep-equal (MIT). */
export function deepEqual(a: any, b: any): boolean {
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

/**
 * Simple deep clone for number, boolean, string, null, undefined, object, array, typed array, Map and Set
 * (a Map's keys are kept, its values cloned). Other class instances are copied as plain objects of their own keys.
 */
export function deepClone<T>(source: T): T {
  if (null === source || 'object' !== typeof source) return source;

  // typed arrays keep their type (and DataView its bytes)
  if (ArrayBuffer.isView(source)) {
    if (source instanceof DataView)
      return new DataView(source.buffer.slice(source.byteOffset, source.byteOffset + source.byteLength)) as any as T;
    return (source as any).slice() as T;
  }
  if (source instanceof Map) {
    const copy = new Map();
    source.forEach((v, k) => copy.set(k, deepClone(v)));
    return copy as any as T;
  }
  if (source instanceof Set) {
    const copy = new Set();
    source.forEach((v) => copy.add(deepClone(v)));
    return copy as any as T;
  }

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
export function pickObjectKeysWithRemapping<V>(
  obj: Record<string, V>,
  keys: string[],
  remapping: Record<string, string | null>,
): Record<string, V> {
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

/** Return `true` if `value` is not `undefined` or `null`.
 * Prefer this over `value !== undefined`
 * (for maybe if we want to allow `null` in `AnnotationRow` in the future) */
export function isDefined<T>(value: T | undefined | null): value is T {
  return value !== undefined && value !== null;
}
