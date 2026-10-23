/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ArrayCtor } from '@molstar/core/util/type-helpers';

export function arrayFind<T>(array: ArrayLike<T>, f: (v: T) => boolean): T | undefined {
  for (let i = 0, _i = array.length; i < _i; i++) {
    if (f(array[i])) return array[i];
  }
  return void 0;
}

export function iterableToArray<T>(it: IterableIterator<T>): T[] {
  if (Array.from) return Array.from(it);

  const ret = [];
  while (true) {
    const { done, value } = it.next();
    if (done) break;
    ret[ret.length] = value;
  }
  return ret;
}

/** Fills the array so that array[0] = start and array[array.length - 1] = end */
export function createRangeArray(start: number, end: number, ctor?: ArrayCtor<number>) {
  const len = end - start + 1;
  const array = ctor ? new ctor(len) : new Int32Array(len);
  for (let i = 0; i < len; i++) {
    array[i] = i + start;
  }
  return array;
}

/**
 * `array[indices[i]]` for each i, in an array of `array`'s kind. The common kinds have a loop each, so each keeps its
 * own type feedback: gathering columns of every kind (a sorted `atom_site`: floats, ints and strings) does not make one
 * loop's keyed loads and stores polymorphic for all of them (perf-audit-2 01/P10).
 */
export function arrayPickIndices<T>(
  array: ArrayLike<T>,
  indices: ArrayLike<number>,
): { [i: number]: T; length: number } {
  if (array instanceof Float32Array) return pickFloat32(array, indices) as any;
  if (array instanceof Float64Array) return pickFloat64(array, indices) as any;
  if (array instanceof Int32Array) return pickInt32(array, indices) as any;
  if (array instanceof Uint8Array) return pickUint8(array, indices) as any;
  if (array instanceof Uint16Array) return pickUint16(array, indices) as any;
  if (Array.isArray(array)) return pickArray(array, indices);
  return pickAny(array, indices);
}

function pickFloat32(array: Float32Array, indices: ArrayLike<number>) {
  const ret = new Float32Array(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

function pickFloat64(array: Float64Array, indices: ArrayLike<number>) {
  const ret = new Float64Array(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

function pickInt32(array: Int32Array, indices: ArrayLike<number>) {
  const ret = new Int32Array(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

function pickUint8(array: Uint8Array, indices: ArrayLike<number>) {
  const ret = new Uint8Array(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

function pickUint16(array: Uint16Array, indices: ArrayLike<number>) {
  const ret = new Uint16Array(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

function pickArray<T>(array: readonly T[], indices: ArrayLike<number>) {
  const ret = new Array<T>(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

function pickAny<T>(array: ArrayLike<T>, indices: ArrayLike<number>) {
  const ret = new (arrayGetCtor(array))(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) ret[i] = array[indices[i]];
  return ret;
}

export function arrayGetCtor<T>(data: ArrayLike<T>): ArrayCtor<T> {
  const ret = (data as any).constructor;
  if (!ret) throw new Error('data does not define a constructor and it should');
  return ret;
}
