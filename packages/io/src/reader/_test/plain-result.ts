/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * A parser's result as plain data to compare with `toEqual`: a column as the array of its values, functions dropped.
 * Two results are equal when their values are, whatever backs their columns (a string or bytes).
 */
export function plainResult(x: unknown): unknown {
  if (Array.isArray(x)) return x.map(plainResult);
  if (x === null || typeof x !== 'object' || ArrayBuffer.isView(x)) return x;
  const o = x as Record<string, unknown>;
  if (typeof o.value === 'function' && typeof o.rowCount === 'number') {
    const value = o.value as (row: number) => unknown;
    return Array.from({ length: o.rowCount }, (_, i) => value(i));
  }
  const ret: Record<string, unknown> = {};
  for (const k of Object.keys(o)) {
    if (typeof o[k] !== 'function') ret[k] = plainResult(o[k]);
  }
  return ret;
}
