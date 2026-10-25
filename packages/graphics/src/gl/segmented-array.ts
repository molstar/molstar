/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The typed arrays a segmented array is made of: vertex attributes (`Float32Array`), elements (`Uint32Array`) and data
 * images (`Uint8Array`, `Float32Array`)
 */
export type SegmentedArrayType = Float32Array | Uint32Array | Uint8Array;

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

  /**
   * The elements `[start, end)` of `a` composed from its parts (zeros between them): for writers of whole words or rows,
   * which a part may share with its neighbors
   */
  export function read<T extends SegmentedArrayType>(a: SegmentedArray<T>, start: number, end: number): T {
    const out = new a.ctor(end - start);
    const { parts } = a;
    for (let i = 0, il = parts.length; i < il; ++i) {
      const p = parts[i];
      const s = Math.max(start, p.offset);
      const e = Math.min(end, p.offset + p.count);
      if (e <= s) continue;
      const data = p.array.subarray(s - p.offset, e - p.offset);
      if (p.bias === 0) out.set(data, s - start);
      else for (let k = 0, kl = e - s; k < kl; ++k) out[s - start + k] = data[k] + p.bias;
    }
    return out;
  }

  /** What an uploader has written of a segmented array: the array and its parts' versions */
  export type Uploaded = { array: SegmentedArray; versions: number[] };

  export function uploaded(a: SegmentedArray): Uploaded {
    return { array: a, versions: a.parts.map((p) => p.version) };
  }

  /**
   * The element ranges (`[start, end)` pairs, by offset) of the parts of `a` changed since `last` was uploaded, and
   * notes them uploaded in `last`; undefined (all) where `last` is of another array
   */
  export function changedRanges(a: SegmentedArray, last: Uploaded | undefined): number[] | undefined {
    if (!last || last.array !== a) return undefined;
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
    if (!last || last.array !== a) {
      forEachPart(a, write);
      return uploaded(a);
    }
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
