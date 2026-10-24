/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type Mutable<T> = { -readonly [P in keyof T]: T[P] };
export type PartialBy<T, K extends keyof T> = Omit<T, K> & Partial<Pick<T, K>>;
export type PickRequired<T, K extends keyof T> = T & Required<Pick<T, K>>;

export type TypedIntArray = Int8Array | Int16Array | Int32Array | Uint8Array | Uint16Array | Uint32Array;
export type TypedFloatArray = Float32Array | Float64Array;

export type TypedArray = TypedIntArray | TypedFloatArray;
export type NumberArray = TypedArray | number[];
export type UintArray = Uint8Array | Uint16Array | Uint32Array | number[];
export type ValueOf<T> = T[keyof T];
export type ArrayCtor<T> = { new (size: number): { [i: number]: T; length: number } };
/**
 * A `Uint8Array` that its reader may not write into: no index assignment, `set`, `fill`, `copyWithin`, `reverse` or
 * `sort`, and its `subarray` is read-only too. For bytes shared with others (a response's, held by the asset manager),
 * checked at the type level only: a `Uint8Array` is one, nothing is frozen, and `buffer` still reaches the memory.
 */
export interface ReadonlyUint8Array
  extends Omit<Uint8Array, number | 'set' | 'fill' | 'copyWithin' | 'reverse' | 'sort' | 'subarray'> {
  readonly [index: number]: number;
  subarray(begin?: number, end?: number): ReadonlyUint8Array;
}

export namespace ReadonlyUint8Array {
  /** Whether `x` is bytes (a `Uint8Array`): narrows `StringLike | ReadonlyUint8Array`, which `instanceof` does not. */
  export function is(x: unknown): x is ReadonlyUint8Array {
    return x instanceof Uint8Array;
  }
}

/** assignable ArrayLike version */
export type AssignableArrayLike<T> = { [i: number]: T; length: number };

export type NonNullableArray<T extends any[] | ReadonlyArray<any>> = T extends any[]
  ? NonNullable<T[0]>[]
  : ReadonlyArray<NonNullable<T[0]>>;

export function ObjectKeys<T extends object>(o: T) {
  return Object.keys(o) as (keyof T)[];
}
export interface FiniteArray<T, L extends number = number> extends ReadonlyArray<T> {
  length: L;
}

export function assertUnreachable(x: never): never {
  throw new Error('unreachable');
}

export function isPromiseLike<T = any>(x: any): x is Promise<T> {
  return typeof x?.then === 'function';
}
