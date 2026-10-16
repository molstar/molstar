/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

function isPlainObject(value: object): value is Record<string, unknown> {
  return Object.getPrototypeOf(value) === Object.prototype;
}

export function assertStrictJsonData(value: unknown, depth = 0): void {
  if (depth > 32) throw new Error('data is too deeply nested');
  if (value === null || typeof value === 'string' || typeof value === 'boolean') return;
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) throw new Error('data contains a non-finite number');
    return;
  }
  if (Array.isArray(value)) {
    for (const item of value) assertStrictJsonData(item, depth + 1);
    return;
  }
  if (typeof value === 'object' && value !== null && isPlainObject(value)) {
    for (const [key, item] of Object.entries(value)) {
      if (item === undefined) throw new Error(`data contains undefined at ${key}`);
      assertStrictJsonData(item, depth + 1);
    }
    return;
  }
  throw new Error('data must contain only plain JSON-compatible values');
}

export function cloneData<T>(value: T, depth = 0): T {
  if (depth > 64) throw new Error('data is too deeply nested');
  if (value === null || typeof value === 'string' || typeof value === 'boolean') return value;
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) throw new Error('data contains a non-finite number');
    return value;
  }
  if (value === undefined) return value;
  if (Array.isArray(value)) return value.map((item) => cloneData(item, depth + 1)) as T;
  if (typeof value === 'object' && isPlainObject(value)) {
    return Object.fromEntries(
      Object.entries(value)
        .filter(([, item]) => item !== undefined)
        .map(([key, item]) => [key, cloneData(item, depth + 1)]),
    ) as T;
  }
  throw new Error('data must contain only plain JSON-compatible values');
}
