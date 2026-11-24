/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GpuDataWriter } from './renderable/util.js';

/**
 * The typed arrays a segmented array is made of: vertex attributes (`Float32Array`, encoded normals `Int16Array`),
 * elements (`Uint32Array`) and data images (`Uint8Array`, `Float32Array`)
 */
export type SegmentedArrayType = Float32Array | Uint32Array | Int32Array | Uint8Array | Int16Array;

type ArrayCtor<T extends SegmentedArrayType> = { new (length: number): T; readonly BYTES_PER_ELEMENT: number };

/**
 * An array made of parts of other arrays, uploaded part by part instead of copied into one array first: the vertex
 * attributes, elements and theme layers of merged render objects are their members' arrays. A render item writes the
 * parts into its buffer or texture at their offsets, and on an update (the cell's version changed) the parts whose
 * version changed; a change of the parts' layout is a new segmented array
 */
export interface SegmentedArray<T extends SegmentedArrayType = SegmentedArrayType> {
  readonly kind: 'segmented-array';
  readonly ctor: ArrayCtor<T>;
  readonly BYTES_PER_ELEMENT: number;
  /** elements in all */
  readonly length: number;
  readonly byteLength: number;
  /** by offset, not overlapping */
  readonly parts: readonly SegmentedArray.Part<T>[];
  /**
   * the array this one extends (`SegmentedArray.extend`): its parts this one's first, the same objects, its length or
   * shorter; an uploader of it writes the parts added (and those changed) alone
   */
  readonly from?: SegmentedArray<T>;
}

export namespace SegmentedArray {
  export interface Part<T extends SegmentedArrayType = SegmentedArrayType> {
    /** the source: its first `count` elements */
    array: T;
    /** elements before the part */
    readonly offset: number;
    readonly count: number;
    /** added to each element (elements: the member's first vertex), 0 for none */
    readonly bias: number;
    /** advanced where the part's data changed */
    version: number;
    /**
     * what writes the part's data on the GPU (a member's `DataArray.gpu`): `array` zeros, the uploader writing the part
     * by it where it wrote its range or the words around it
     */
    gpu?: GpuDataWriter;
  }

  export function is(x: unknown): x is SegmentedArray {
    return !!x && (x as SegmentedArray).kind === 'segmented-array';
  }

  /** `length` at least the parts' end (more leaves zeros after them) */
  export function create<T extends SegmentedArrayType>(
    ctor: ArrayCtor<T>,
    parts: readonly Part<T>[],
    length = 0,
  ): SegmentedArray<T> {
    for (const p of parts) length = Math.max(length, p.offset + p.count);
    const { BYTES_PER_ELEMENT } = ctor;
    return { kind: 'segmented-array', ctor, BYTES_PER_ELEMENT, length, byteLength: length * BYTES_PER_ELEMENT, parts };
  }

  /** Extensions kept in a line (`from`): an uploader further behind writes all */
  const MaxExtensions = 16;

  /**
   * `a` with `parts` added after its own (by offset, within `length`): a new array extending `a` (`from`), so that what
   * uploaded `a` writes the added parts alone; `length` at least `a`'s (its by default: longer, room to grow, zeros
   * after the parts, its parts where they were)
   */
  export function extend<T extends SegmentedArrayType>(
    a: SegmentedArray<T>,
    parts: readonly Part<T>[],
    length = a.length,
  ): SegmentedArray<T> {
    if (length < a.length) throw new Error('segmented array: an extension shorter');
    let end = 0;
    for (const p of a.parts) end = Math.max(end, p.offset + p.count);
    for (const p of parts) {
      if (p.offset < end || p.offset + p.count > length)
        throw new Error('segmented array: a part not after the others');
      end = p.offset + p.count;
    }
    // the line cut behind its last extensions (the arrays before kept alive by it otherwise), not at `a`: an uploader of
    // `a` (one update behind, the usual) still writes the added parts alone, one further behind all
    let last: SegmentedArray<T> = a;
    for (let k = 1; k < MaxExtensions && last.from; ++k) last = last.from;
    if (last.from) (last as { from?: SegmentedArray<T> }).from = undefined;
    const byteLength = length * a.BYTES_PER_ELEMENT;
    return { ...a, length, byteLength, parts: [...a.parts, ...parts], from: a };
  }

  /** whether `a` is `b` or extends it (`from`) */
  function isExtensionOf(a: SegmentedArray, b: SegmentedArray) {
    for (let x: SegmentedArray | undefined = a; x; x = x.from) if (x === b) return true;
    return false;
  }

  /** The part's data, biased where it has a bias (into a new array: kept only for the write) */
  export function getPartData<T extends SegmentedArrayType>(part: Part<T>): T {
    const data = part.array.subarray(0, part.count) as T;
    if (part.bias === 0) return data;
    const biased = new (data.constructor as new (length: number) => T)(part.count);
    for (let i = 0, il = part.count; i < il; ++i) biased[i] = data[i] + part.bias;
    return biased;
  }

  /** Calls `write` for each part with its data and element offset */
  export function forEachPart<T extends SegmentedArrayType>(
    a: SegmentedArray<T>,
    write: (data: T, offset: number) => void,
  ) {
    for (const p of a.parts) if (p.count > 0) write(getPartData(p), p.offset);
  }

  /** The index of the first of `parts` ending after `index` (by offset, not overlapping: their ends ascend) */
  function findFirstPartEndingAfter(parts: readonly Part[], index: number) {
    let lo = 0;
    let hi = parts.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (parts[mid].offset + parts[mid].count <= index) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  /**
   * The elements `[start, end)` of `a` composed from its parts (zeros between them): for writers of whole words or rows,
   * which a part may share with its neighbors
   */
  export function read<T extends SegmentedArrayType>(a: SegmentedArray<T>, start: number, end: number): T {
    const out = new a.ctor(end - start);
    const { parts } = a;
    // from the first part ending after `start`, to the first starting at `end`: a read per part of many (a merged
    // object's members) stays linear
    for (let i = findFirstPartEndingAfter(parts, start), il = parts.length; i < il; ++i) {
      const p = parts[i];
      if (p.offset >= end) break;
      const s = Math.max(start, p.offset);
      const e = Math.min(end, p.offset + p.count);
      if (e <= s) continue;
      const data = p.array.subarray(s - p.offset, e - p.offset);
      if (p.bias === 0) out.set(data, s - start);
      else for (let k = 0, kl = e - s; k < kl; ++k) out[s - start + k] = data[k] + p.bias;
    }
    return out;
  }

  /** What an uploader has written of a segmented array: the array and its parts' versions (none: a part not written) */
  export type Uploaded = { array: SegmentedArray; versions: number[] };

  export function uploaded(a: SegmentedArray): Uploaded {
    return { array: a, versions: a.parts.map((p) => p.version) };
  }

  /**
   * The element ranges (`[start, end)` pairs, by offset) of the parts of `a` changed since `last` was uploaded (those
   * an extension added among them), and notes them uploaded in `last`; undefined (all) where `last` is of another
   * array
   */
  export function changedRanges(a: SegmentedArray, last: Uploaded | undefined): number[] | undefined {
    if (!last || !isExtensionOf(a, last.array)) return undefined;
    last.array = a;
    const ranges: number[] = [];
    const { parts } = a;
    for (let i = 0, il = parts.length; i < il; ++i) {
      const p = parts[i];
      if (p.version === last.versions[i]) continue;
      if (p.count > 0) ranges.push(p.offset, p.offset + p.count);
      last.versions[i] = p.version;
    }
    return ranges;
  }

  /**
   * Calls `write` for the parts of `a` changed since `last` was uploaded, and returns what is uploaded then; all parts
   * where `last` is of another array
   */
  export function writeChanged<T extends SegmentedArrayType>(
    a: SegmentedArray<T>,
    last: Uploaded | undefined,
    write: (data: T, offset: number) => void,
  ): Uploaded {
    if (!last || !isExtensionOf(a, last.array)) {
      forEachPart(a, write);
      return uploaded(a);
    }
    last.array = a;
    const { parts } = a;
    for (let i = 0, il = parts.length; i < il; ++i) {
      const p = parts[i];
      if (p.version === last.versions[i]) continue;
      if (p.count > 0) write(getPartData(p), p.offset);
      last.versions[i] = p.version;
    }
    return last;
  }

  /** The array in one piece (for consumers that need it whole, e.g. interleaving) */
  export function flatten<T extends SegmentedArrayType>(a: SegmentedArray<T>): T {
    return read(a, 0, a.length);
  }
}
