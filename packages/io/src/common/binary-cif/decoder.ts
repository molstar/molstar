/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * From CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Encoding, type EncodedData } from './encoding.js';
import { IsNativeEndianLittle, flipByteOrder } from '../binary.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';

/**
 * Fixed point, delta, RLE, integer packing adopted from https://github.com/rcsb/mmtf-javascript/
 * by Alexander Rose <alexander.rose@weirdbyte.de>, MIT License, Copyright (c) 2016
 */

export function decode(data: EncodedData): any[] {
  return (decodeFused(data) as any[]) ?? decodeChain(data, true);
}

/** The decode one step at a time, an array per step: the reference for the fused decode. */
export function decodeSteps(data: EncodedData): any[] {
  return decodeChain(data, false);
}

function decodeChain(data: EncodedData, fuse: boolean): any[] {
  let current: any = data.data;
  for (let i = data.encoding.length - 1; i >= 0; i--) {
    current = decodeStep(current, data.encoding[i], fuse);
  }
  return current as any[];
}

function decodeStep(data: any, encoding: Encoding, fuse = false): any {
  switch (encoding.kind) {
    case 'ByteArray': {
      switch (encoding.type) {
        case Encoding.IntDataType.Uint8:
          return data;
        case Encoding.IntDataType.Int8:
          return int8(data);
        case Encoding.IntDataType.Int16:
          return int16(data);
        case Encoding.IntDataType.Uint16:
          return uint16(data);
        case Encoding.IntDataType.Int32:
          return int32(data);
        case Encoding.IntDataType.Uint32:
          return uint32(data);
        case Encoding.FloatDataType.Float32:
          return float32(data);
        case Encoding.FloatDataType.Float64:
          return float64(data);
        default:
          assertUnreachable(encoding.type);
      }
    }
    case 'FixedPoint':
      return fixedPoint(data, encoding);
    case 'IntervalQuantization':
      return intervalQuantization(data, encoding);
    case 'RunLength':
      return runLength(data, encoding);
    case 'Delta':
      return delta(data, encoding);
    case 'IntegerPacking':
      return integerPacking(data, encoding);
    case 'StringArray':
      return stringArray(data, encoding, fuse);
  }
}

function getIntArray(type: Encoding.IntDataType, size: number) {
  switch (type) {
    case Encoding.IntDataType.Int8:
      return new Int8Array(size);
    case Encoding.IntDataType.Int16:
      return new Int16Array(size);
    case Encoding.IntDataType.Int32:
      return new Int32Array(size);
    case Encoding.IntDataType.Uint8:
      return new Uint8Array(size);
    case Encoding.IntDataType.Uint16:
      return new Uint16Array(size);
    case Encoding.IntDataType.Uint32:
      return new Uint32Array(size);
    default:
      return new Int32Array(size);
  }
}

function getFloatArray(type: Encoding.FloatDataType, size: number) {
  switch (type) {
    case Encoding.FloatDataType.Float32:
      return new Float32Array(size);
    case Encoding.FloatDataType.Float64:
      return new Float64Array(size);
    default:
      return new Float64Array(size);
  }
}

function int8(data: Uint8Array) {
  return new Int8Array(data.buffer, data.byteOffset);
}

function view<T>(data: Uint8Array, byteSize: number, c: new (buffer: ArrayBufferLike) => T) {
  if (IsNativeEndianLittle) return new c(data.buffer);
  return new c(flipByteOrder(data, byteSize));
}

function int16(data: Uint8Array) {
  return view(data, 2, Int16Array);
}
function uint16(data: Uint8Array) {
  return view(data, 2, Uint16Array);
}
function int32(data: Uint8Array) {
  return view(data, 4, Int32Array);
}
function uint32(data: Uint8Array) {
  return view(data, 4, Uint32Array);
}
function float32(data: Uint8Array) {
  return view(data, 4, Float32Array);
}
function float64(data: Uint8Array) {
  return view(data, 8, Float64Array);
}

function fixedPoint(data: Int32Array, encoding: Encoding.FixedPoint) {
  const n = data.length;
  const output = getFloatArray(encoding.srcType, n);
  const f = 1 / encoding.factor;
  for (let i = 0; i < n; i++) {
    output[i] = f * data[i];
  }
  return output;
}

function intervalQuantization(data: Int32Array, encoding: Encoding.IntervalQuantization) {
  const n = data.length;
  const output = getFloatArray(encoding.srcType, n);
  const delta = (encoding.max - encoding.min) / (encoding.numSteps - 1);
  const min = encoding.min;
  for (let i = 0; i < n; i++) {
    output[i] = min + delta * data[i];
  }
  return output;
}

function runLength(data: Int32Array, encoding: Encoding.RunLength) {
  const output = getIntArray(encoding.srcType, encoding.srcSize);
  let dataOffset = 0;
  for (let i = 0, il = data.length; i < il; i += 2) {
    const value = data[i]; // value to be repeated
    const length = data[i + 1]; // number of repeats
    for (let j = 0; j < length; ++j) {
      output[dataOffset++] = value;
    }
  }
  return output;
}

function delta(data: Int8Array | Int16Array | Int32Array, encoding: Encoding.Delta) {
  const n = data.length;
  const output = getIntArray(encoding.srcType, n);
  if (!n) return data;
  output[0] = data[0] + (encoding.origin | 0);
  for (let i = 1; i < n; ++i) {
    output[i] = data[i] + output[i - 1];
  }
  return output;
}

function integerPackingSigned(data: Int8Array | Int16Array, encoding: Encoding.IntegerPacking) {
  const upperLimit = encoding.byteCount === 1 ? 0x7f : 0x7fff;
  const lowerLimit = -upperLimit - 1;
  const n = data.length;
  const output = new Int32Array(encoding.srcSize);
  let i = 0;
  let j = 0;
  while (i < n) {
    let value = 0,
      t = data[i];
    while (t === upperLimit || t === lowerLimit) {
      value += t;
      i++;
      t = data[i];
    }
    value += t;
    output[j] = value;
    i++;
    j++;
  }
  return output;
}

function integerPackingUnsigned(data: Uint8Array | Uint16Array, encoding: Encoding.IntegerPacking) {
  const upperLimit = encoding.byteCount === 1 ? 0xff : 0xffff;
  const n = data.length;
  const output = new Int32Array(encoding.srcSize);
  let i = 0;
  let j = 0;
  while (i < n) {
    let value = 0,
      t = data[i];
    while (t === upperLimit) {
      value += t;
      i++;
      t = data[i];
    }
    value += t;
    output[j] = value;
    i++;
    j++;
  }
  return output;
}

function integerPacking(data: Int8Array | Int16Array | Uint8Array | Uint16Array, encoding: Encoding.IntegerPacking) {
  if (data.length === encoding.srcSize) return data;
  return encoding.isUnsigned
    ? integerPackingUnsigned(data as any, encoding)
    : integerPackingSigned(data as any, encoding);
}

/** `''` and the strings of a StringArray: index `i` of the data is `strings[i + 1]` (-1 is `''`). */
function stringArrayStrings(encoding: Encoding.StringArray, fuse: boolean): string[] {
  const offsetData = { encoding: encoding.offsetEncoding, data: encoding.offsets };
  const offsets = fuse ? decode(offsetData) : decodeSteps(offsetData);

  const str = encoding.stringData;
  const strings = new Array(offsets.length);
  strings[0] = '';
  for (let i = 1, _i = offsets.length; i < _i; i++) {
    strings[i] = str.substring(offsets[i - 1], offsets[i]);
  }
  return strings;
}

/**
 * A StringArray column as a dictionary, without its array of strings: `values` is `''` and the distinct strings, and
 * `values[indices[i]]` is the decoded value `i` (the indices of the data plus one, in the smallest unsigned type that
 * holds them). `undefined` for any other encoding, or for an index out of the strings' range (whose decoded value is
 * `undefined`).
 */
export function decodeStringDictionary(
  data: EncodedData,
): { indices: Uint8Array | Uint16Array | Int32Array; values: string[] } | undefined {
  if (data.encoding.length !== 1 || data.encoding[0].kind !== 'StringArray') return;
  const encoding = data.encoding[0];
  const values = stringArrayStrings(encoding, true);
  const raw = decode({ encoding: encoding.dataEncoding, data: data.data }) as ArrayLike<number>;
  const n = raw.length;
  const max = values.length;
  const indices = max <= 0x100 ? new Uint8Array(n) : max <= 0x10000 ? new Uint16Array(n) : new Int32Array(n);
  for (let i = 0; i < n; i++) {
    const k = raw[i] + 1;
    if (!(k >= 0 && k < max)) return;
    indices[i] = k;
  }
  return { indices, values };
}

function stringArray(data: Uint8Array, encoding: Encoding.StringArray, fuse: boolean) {
  const strings = stringArrayStrings(encoding, fuse);

  // the indices decoded and mapped to their strings in one pass when the chain allows
  const fused = fuse ? decodeFused({ encoding: encoding.dataEncoding, data }, strings) : void 0;
  if (fused) return fused;

  const indices = decodeSteps({ encoding: encoding.dataEncoding, data });
  let offset = 0;
  const result = new Array(indices.length);
  for (let i = 0, _i = indices.length; i < _i; i++) {
    result[offset++] = strings[indices[i] + 1];
  }
  return result;
}

/**
 * One pass for the common integer chains, `[FixedPoint | IntervalQuantization] > [Delta] > [RunLength] >
 * [IntegerPacking] > ByteArray` with at least two steps above the byte array: coordinates, B-factors, ids, sequence ids
 * and the index and offset arrays of string columns. No intermediate array is made.
 *
 * Bit-identical to the step-by-step decode: values wrap as the stores into the Int32 intermediates do (`| 0`), and a
 * short or overlong input reads as zeros or is cut as the typed arrays of the steps do. Chains with non-Int32
 * intermediates or no rows return `undefined` and take the generic path. With `lookup` (the strings of a StringArray),
 * the output is `lookup[value + 1]` per row, as `stringArray` maps its indices.
 */
function decodeFused(data: EncodedData, lookup?: ArrayLike<string>): ArrayLike<any> | undefined {
  const encoding = data.encoding;
  let e = encoding.length - 1;
  const byteArray = encoding[e];
  if (e < 0 || byteArray.kind !== 'ByteArray' || !(byteArray.type >= 1 && byteArray.type <= 6)) return;
  e--;
  let packing: Encoding.IntegerPacking | undefined;
  let runLength: Encoding.RunLength | undefined;
  let delta: Encoding.Delta | undefined;
  let scale: Encoding.FixedPoint | Encoding.IntervalQuantization | undefined;
  let steps = 0;
  if (e >= 0 && encoding[e].kind === 'IntegerPacking') (packing = encoding[e--] as Encoding.IntegerPacking), steps++;
  if (e >= 0 && encoding[e].kind === 'RunLength') (runLength = encoding[e--] as Encoding.RunLength), steps++;
  if (e >= 0 && encoding[e].kind === 'Delta') (delta = encoding[e--] as Encoding.Delta), steps++;
  if (e >= 0 && (encoding[e].kind === 'FixedPoint' || encoding[e].kind === 'IntervalQuantization')) {
    scale = encoding[e--] as Encoding.FixedPoint | Encoding.IntervalQuantization;
    steps++;
  }
  // one step is already one pass, unless the output is mapped to strings
  if (e !== -1 || (!lookup && steps < 2)) return;
  if (runLength && runLength.srcType !== Encoding.IntDataType.Int32) return;
  if (delta && delta.srcType !== Encoding.IntDataType.Int32) return;

  const src = decodeStep(data.data, byteArray) as ArrayLike<number>;
  const n = src.length;
  // `integerPacking` returns its input when the lengths match
  const unpack = !!packing && n !== packing.srcSize;
  /** length of the integer stream: the IntegerPacking output or the byte array */
  const m = unpack ? packing!.srcSize : n;
  const count = runLength ? runLength.srcSize : m;
  if (!(count > 0)) return;

  const plan: FusedPlan = { src, n, unpack, m, count, runLength: !!runLength, delta, upper: 0, lower: 0 };
  if (packing) {
    plan.upper = packing.isUnsigned
      ? packing.byteCount === 1
        ? 0xff
        : 0xffff
      : packing.byteCount === 1
        ? 0x7f
        : 0x7fff;
    plan.lower = packing.isUnsigned ? plan.upper : -plan.upper - 1;
  }
  if (lookup) return scale ? void 0 : fusedStrings(plan, lookup);
  if (!scale) return fusedNumbers(plan, new Int32Array(count), 0, 0, 0);
  return scale.kind === 'FixedPoint'
    ? fusedNumbers(plan, getFloatArray(scale.srcType, count), 1, 1 / scale.factor, 0)
    : fusedNumbers(
        plan,
        getFloatArray(scale.srcType, count),
        2,
        (scale.max - scale.min) / (scale.numSteps - 1),
        scale.min,
      );
}

interface FusedPlan {
  src: ArrayLike<number>;
  /** length of `src` */
  n: number;
  /** IntegerPacking to undo */
  unpack: boolean;
  /** length of the integer stream: the IntegerPacking output or the byte array */
  m: number;
  /** length of the output */
  count: number;
  runLength: boolean;
  delta: Encoding.Delta | undefined;
  upper: number;
  lower: number;
}

// Two copies of one loop, numbers and strings: one function writing typed arrays and a JS array went megamorphic and
// was slower than the step decoder.

/** `kind` 0: the integers, 1: FixedPoint `factor * v`, 2: IntervalQuantization `min + factor * v` */
function fusedNumbers(
  plan: FusedPlan,
  out: Int32Array | Float32Array | Float64Array,
  kind: 0 | 1 | 2,
  factor: number,
  min: number,
) {
  const { src, n, unpack, m, count, runLength, upper, lower } = plan;
  const hasDelta = !!plan.delta;
  let acc = plan.delta ? plan.delta.origin | 0 : 0;
  // the integer stream, inlined below (closures over these counters were 2x slower than the step decoder): `i` is the
  // position in `src`; past the end of the input the IntegerPacking output holds zeros
  let i = 0;
  let j = 0;
  let v = 0;
  let t = 0;

  if (runLength) {
    for (let k = 0; k < m && j < count; ) {
      if (!unpack) v = src[i++];
      else if (i >= n) v = 0;
      else {
        for (v = 0, t = src[i]; t === upper || t === lower; t = src[++i]) v += t;
        v = i < n ? (v + src[i++]) | 0 : 0;
      }
      const value = v | 0;
      if (++k >= m) break;
      if (!unpack) v = src[i++];
      else if (i >= n) v = 0;
      else {
        for (v = 0, t = src[i]; t === upper || t === lower; t = src[++i]) v += t;
        v = i < n ? (v + src[i++]) | 0 : 0;
      }
      k++;
      for (let length = v; length > 0 && j < count; length--) {
        v = hasDelta ? (acc = (value + acc) | 0) : value;
        out[j++] = kind === 0 ? v : kind === 1 ? factor * v : min + factor * v;
      }
    }
    for (; j < count; j++) {
      v = hasDelta ? acc : 0;
      out[j] = kind === 0 ? v : kind === 1 ? factor * v : min + factor * v;
    }
  } else {
    for (; j < count; j++) {
      if (!unpack) v = src[i++];
      else if (i >= n) v = 0;
      else {
        for (v = 0, t = src[i]; t === upper || t === lower; t = src[++i]) v += t;
        v = i < n ? (v + src[i++]) | 0 : 0;
      }
      if (hasDelta) v = acc = (v + acc) | 0;
      out[j] = kind === 0 ? v : kind === 1 ? factor * v : min + factor * v;
    }
  }
  return out;
}

/** `stringArray`'s `strings[index + 1]` per row */
function fusedStrings(plan: FusedPlan, strings: ArrayLike<string>) {
  const { src, n, unpack, m, count, runLength, upper, lower } = plan;
  const out: string[] = new Array(count);
  const hasDelta = !!plan.delta;
  let acc = plan.delta ? plan.delta.origin | 0 : 0;
  // the integer stream, inlined below (closures over these counters were 2x slower than the step decoder): `i` is the
  // position in `src`; past the end of the input the IntegerPacking output holds zeros
  let i = 0;
  let j = 0;
  let v = 0;
  let t = 0;

  if (runLength) {
    for (let k = 0; k < m && j < count; ) {
      if (!unpack) v = src[i++];
      else if (i >= n) v = 0;
      else {
        for (v = 0, t = src[i]; t === upper || t === lower; t = src[++i]) v += t;
        v = i < n ? (v + src[i++]) | 0 : 0;
      }
      const value = v | 0;
      if (++k >= m) break;
      if (!unpack) v = src[i++];
      else if (i >= n) v = 0;
      else {
        for (v = 0, t = src[i]; t === upper || t === lower; t = src[++i]) v += t;
        v = i < n ? (v + src[i++]) | 0 : 0;
      }
      k++;
      for (let length = v; length > 0 && j < count; length--) {
        v = hasDelta ? (acc = (value + acc) | 0) : value;
        out[j++] = strings[v + 1];
      }
    }
    for (; j < count; j++) {
      v = hasDelta ? acc : 0;
      out[j] = strings[v + 1];
    }
  } else {
    for (; j < count; j++) {
      if (!unpack) v = src[i++];
      else if (i >= n) v = 0;
      else {
        for (v = 0, t = src[i]; t === upper || t === lower; t = src[++i]) v += t;
        v = i < n ? (v + src[i++]) | 0 : 0;
      }
      if (hasDelta) v = acc = (v + acc) | 0;
      out[j] = strings[v + 1];
    }
  }
  return out;
}
