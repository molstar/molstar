/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export type Comparer<T = any> = (data: T, i: number, j: number) => number;
export type Swapper<T = any> = (data: T, i: number, j: number) => void;

type Ctx = { cmp: Comparer; swap: Swapper; parts: number[]; data: any };

export function arrayLess(arr: ArrayLike<number>, i: number, j: number) {
  return arr[i] - arr[j];
}

export function arraySwap(arr: ArrayLike<any>, i: number, j: number) {
  const temp = arr[i];
  (arr as any[])[i] = arr[j];
  (arr as any[])[j] = temp;
}

function medianPivotIndex(data: any, cmp: Comparer, l: number, r: number) {
  const m = (l + r) >> 1;
  if (cmp(data, l, r) > 0) return cmp(data, l, m) > 0 ? (cmp(data, m, r) > 0 ? m : r) : l;
  else return cmp(data, r, m) > 0 ? (cmp(data, m, l) > 0 ? m : l) : r;
}

function partition(ctx: Ctx, l: number, r: number) {
  const { cmp, swap, data, parts } = ctx;
  let equals = l + 1,
    tail = r;

  // move the median to the 1st spot
  swap(data, l, medianPivotIndex(data, cmp, l, r));

  while (cmp(data, tail, l) > 0) {
    --tail;
  }
  for (let i = l + 1; i <= tail; i++) {
    const c = cmp(data, i, l);
    if (c > 0) {
      swap(data, i, tail);
      --tail;
      while (cmp(data, tail, l) > 0) {
        --tail;
      }
      i--;
    } else if (c === 0) {
      swap(data, i, equals);
      equals++;
    }
  }

  // move the medians to the correct spots
  for (let i = l; i < equals; i++) {
    swap(data, i, l + tail - i);
  }
  parts[0] = tail - equals + l + 1;
  parts[1] = tail;
}

function insertionSort({ data, cmp, swap }: Ctx, start: number, end: number) {
  for (let i = start + 1; i <= end; i++) {
    let j = i - 1;
    while (j >= start && cmp(data, j, j + 1) > 0) {
      swap(data, j, j + 1);
      j = j - 1;
    }
  }
}

function quickSort(ctx: Ctx, low: number, high: number) {
  const { parts } = ctx;
  while (low < high) {
    if (high - low < 16) {
      insertionSort(ctx, low, high);
      return;
    }

    partition(ctx, low, high);
    const li = parts[0],
      ri = parts[1];

    if (li - low < high - ri) {
      quickSort(ctx, low, li - 1);
      low = ri + 1;
    } else {
      quickSort(ctx, ri + 1, high);
      high = li - 1;
    }
  }
}

function partitionArrayAsc(data: number[], parts: number[], l: number, r: number) {
  let equals = l + 1,
    tail = r;

  // move the median to the 1st spot
  arraySwap(data, l, medianPivotIndex(data, arrayLess, l, r));
  const pivot = data[l];

  while (data[tail] > pivot) {
    --tail;
  }
  for (let i = l + 1; i <= tail; i++) {
    const v = data[i];
    if (v > pivot) {
      arraySwap(data, i, tail);
      --tail;
      while (data[tail] > pivot) {
        --tail;
      }
      i--;
    } else if (v === pivot) {
      arraySwap(data, i, equals);
      ++equals;
    }
  }

  // move all medians to the correct spots
  for (let i = l; i < equals; i++) {
    arraySwap(data, i, l + tail - i);
  }
  parts[0] = tail - equals + l + 1;
  parts[1] = tail;
}

function insertionSortArrayAsc(data: number[], start: number, end: number) {
  for (let i = start + 1; i <= end; i++) {
    const key = data[i];
    let j = i - 1;
    while (j >= start && data[j] > key) {
      data[j + 1] = data[j];
      j = j - 1;
    }
    data[j + 1] = key;
  }
}

function quickSortArrayAsc(data: number[], parts: number[], low: number, high: number) {
  while (low < high) {
    if (high - low < 16) {
      insertionSortArrayAsc(data, low, high);
      return;
    }

    partitionArrayAsc(data, parts, low, high);
    const li = parts[0],
      ri = parts[1];

    if (li - low < high - ri) {
      quickSortArrayAsc(data, parts, low, li - 1);
      low = ri + 1;
    } else {
      quickSortArrayAsc(data, parts, ri + 1, high);
      high = li - 1;
    }
  }
}

export function sortArray(data: ArrayLike<number>, cmp: Comparer<ArrayLike<number>> = arrayLess): ArrayLike<number> {
  return sortArrayRange(data, 0, data.length, cmp);
}

export function sortArrayRange(
  data: ArrayLike<number>,
  start: number,
  end: number,
  cmp: Comparer<ArrayLike<number>> = arrayLess,
): ArrayLike<number> {
  if (cmp === arrayLess) quickSortArrayAsc(data as any, [0, 0], start, end - 1);
  else quickSort({ data, cmp, swap: arraySwap, parts: [0, 0] }, start, end - 1);
  return data;
}

export function sort<T>(data: T, start: number, end: number, cmp: Comparer<T>, swap: Swapper<T>): T {
  const ctx: Ctx = { data, cmp, swap, parts: [0, 0] };
  quickSort(ctx, start, end - 1);
  return data;
}

const _radixCounts = new Uint32Array(257);

/**
 * In-place LSD radix sort (8-bit digits) for non-negative integers, ascending.
 * `maxValue` is an inclusive upper bound of the values and determines the number
 * of passes; it is computed from the data when omitted.
 *
 * Typically outperforms comparison sorts, at the cost of a scratch
 * array allocation of the same size as the input.
 */
export function radixSort<T extends Int32Array | Uint32Array>(data: T, maxValue?: number): T {
  const n = data.length;
  if (n < 2) return data;

  let max = maxValue;
  if (max === undefined) {
    max = 0;
    for (let i = 0; i < n; i++) {
      if (data[i] > max) max = data[i];
    }
  }

  const counts = _radixCounts;
  let a: Int32Array | Uint32Array = data;
  let b: Int32Array | Uint32Array = new (data.constructor as new (n: number) => T)(n);
  for (let shift = 0; shift < 32 && max >>> shift > 0; shift += 8) {
    counts.fill(0);
    for (let i = 0; i < n; i++) counts[((a[i] >>> shift) & 0xff) + 1]++;
    for (let i = 0; i < 256; i++) counts[i + 1] += counts[i];
    for (let i = 0; i < n; i++) b[counts[(a[i] >>> shift) & 0xff]++] = a[i];
    const t = a;
    a = b;
    b = t;
  }
  // odd number of passes ends in the scratch array
  if (a !== data) data.set(a);
  return data;
}

/** from this many keys of more than two digits `radixSortIndices` sorts by bucket (`radixSortIndicesByBucket`) */
const RadixBucketMin = 1 << 16;
const _radixBounds = new Uint32Array(257);

/**
 * The indices `0..n-1` of `keys` (non-negative integers) ordered by their key, ascending, equal keys by index (stable):
 * an LSD radix sort (8-bit digits) of the indices, the keys read in place. `maxValue` is an inclusive upper bound of the
 * keys and determines the number of passes; it is computed from the keys when omitted. Two index arrays of `n`. Many
 * keys of more than two digits by bucket (`radixSortIndicesByBucket`, the same order)
 */
export function radixSortIndices(keys: Int32Array | Uint32Array, maxValue?: number): Int32Array {
  const n = keys.length;
  let max = maxValue;
  if (max === undefined) {
    max = 0;
    for (let i = 0; i < n; i++) {
      if (keys[i] > max) max = keys[i];
    }
  }
  if (n >= RadixBucketMin && max >>> 16 > 0) return radixSortIndicesByBucket(keys, max);

  let a = new Int32Array(n);
  for (let i = 0; i < n; i++) a[i] = i;
  if (n < 2) return a;
  let b = new Int32Array(n);
  const counts = _radixCounts;
  for (let shift = 0; shift < 32 && max >>> shift > 0; shift += 8) {
    counts.fill(0);
    for (let i = 0; i < n; i++) counts[((keys[a[i]] >>> shift) & 0xff) + 1]++;
    for (let i = 0; i < 256; i++) counts[i + 1] += counts[i];
    for (let i = 0; i < n; i++) {
      const k = a[i];
      b[counts[(keys[k] >>> shift) & 0xff]++] = k;
    }
    const t = a;
    a = b;
    b = t;
  }
  return a;
}

/**
 * `radixSortIndices` of many keys (`max` their bound, above two digits): the indices into buckets by the keys' top
 * digit (a pass reading them in order), then each bucket by LSD with its keys alongside, the bucket's arrays in cache;
 * a key read out of order once, not a pass each (31.5 M Morton codes of 27 bits: 3.8 to 1.0 s in Node). Stable, as the LSD
 * sort: a bucket's indices ascending, its passes stable. Scratch of the largest bucket's size
 */
function radixSortIndicesByBucket(keys: Int32Array | Uint32Array, max: number): Int32Array {
  const n = keys.length;
  const top = 32 - Math.clz32(max) - 8;
  const counts = _radixCounts;
  const bounds = _radixBounds;
  counts.fill(0);
  for (let i = 0; i < n; i++) counts[(keys[i] >>> top) + 1]++;
  for (let i = 0; i < 256; i++) counts[i + 1] += counts[i];
  bounds.set(counts);
  let largest = 0;
  for (let d = 0; d < 256; d++) largest = Math.max(largest, bounds[d + 1] - bounds[d]);
  const byBucket = new Int32Array(n);
  for (let i = 0; i < n; i++) byBucket[counts[keys[i] >>> top]++] = i;

  const out = new Int32Array(n);
  const kA = new Uint32Array(largest);
  const kB = new Uint32Array(largest);
  const iA = new Int32Array(largest);
  const iB = new Int32Array(largest);
  const mask = 2 ** top - 1;
  for (let d = 0; d < 256; d++) {
    const s = bounds[d];
    const m = bounds[d + 1] - s;
    if (m < 2) {
      if (m === 1) out[s] = byBucket[s];
      continue;
    }
    let ka = kA;
    let kb = kB;
    let ia = iA;
    let ib = iB;
    for (let j = 0; j < m; j++) {
      const k = byBucket[s + j];
      ia[j] = k;
      ka[j] = keys[k] & mask;
    }
    for (let shift = 0; shift < top; shift += 8) {
      counts.fill(0);
      for (let j = 0; j < m; j++) counts[((ka[j] >>> shift) & 0xff) + 1]++;
      for (let j = 0; j < 256; j++) counts[j + 1] += counts[j];
      for (let j = 0; j < m; j++) {
        const p = counts[(ka[j] >>> shift) & 0xff]++;
        kb[p] = ka[j];
        ib[p] = ia[j];
      }
      const tk = ka;
      ka = kb;
      kb = tk;
      const ti = ia;
      ia = ib;
      ib = ti;
    }
    out.set(ia.subarray(0, m), s);
  }
  return out;
}

let _countingCounts = new Uint32Array(0);

/**
 * The indices `0..n-1` of the first `n` of `keys` (integers from 0 to `maxValue`) ordered by key, ascending, equal keys
 * by index (stable), written into `out` from `offset`: a counting sort, one pass reading the keys in order and a count
 * per value, for keys of a small range (to 2^16 or so: `radixSortIndices` beyond). Faster than `radixSortIndices`'s two
 * passes for 16-bit keys (1 M keys about 2x, 5 M 4x in Node: no reads by index); scratch of `maxValue + 2` counts, kept
 */
export function countingSortIndices(
  keys: ArrayLike<number>,
  n: number,
  maxValue: number,
  out: Uint32Array | Int32Array,
  offset = 0,
): void {
  if (_countingCounts.length < maxValue + 2) _countingCounts = new Uint32Array(maxValue + 2);
  const counts = _countingCounts;
  counts.fill(0, 0, maxValue + 2);
  for (let i = 0; i < n; i++) counts[keys[i] + 1]++;
  for (let k = 0; k <= maxValue; k++) counts[k + 1] += counts[k];
  for (let i = 0; i < n; i++) out[offset + counts[keys[i]]++] = i;
}
