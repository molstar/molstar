/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Ranges of a fixed capacity: first fit from a list of free ranges sorted by offset, merged with their neighbors when
 * freed. Sizes and offsets are multiples of `alignment`, in units of the caller (bytes, elements)
 */
export interface RangeAllocator {
  readonly capacity: number;
  /** units in use */
  readonly used: number;
  /** the offset of a range of `size` units, -1 where none fits */
  allocate(size: number): number;
  free(offset: number, size: number): void;
}

export function createRangeAllocator(capacity: number, alignment: number): RangeAllocator {
  // free ranges by offset, as pairs in one array: [offset, size, offset, size, ...]
  const free: number[] = capacity > 0 ? [0, capacity] : [];
  let used = 0;
  const align = (n: number) => Math.ceil(n / alignment) * alignment;
  return {
    capacity,
    get used() {
      return used;
    },
    allocate(size) {
      const n = align(Math.max(size, 1));
      for (let i = 0, il = free.length; i < il; i += 2) {
        if (free[i + 1] < n) continue;
        const offset = free[i];
        if (free[i + 1] === n) free.splice(i, 2);
        else {
          free[i] += n;
          free[i + 1] -= n;
        }
        used += n;
        return offset;
      }
      return -1;
    },
    free(offset, size) {
      const n = align(Math.max(size, 1));
      used -= n;
      // the first free range after it
      let i = 0;
      while (i < free.length && free[i] < offset) i += 2;
      const mergesBefore = i > 0 && free[i - 2] + free[i - 1] === offset;
      const mergesAfter = i < free.length && offset + n === free[i];
      if (mergesBefore && mergesAfter) {
        free[i - 1] += n + free[i + 1];
        free.splice(i, 2);
      } else if (mergesBefore) {
        free[i - 1] += n;
      } else if (mergesAfter) {
        free[i] = offset;
        free[i + 1] += n;
      } else {
        free.splice(i, 0, offset, n);
      }
    },
  };
}
