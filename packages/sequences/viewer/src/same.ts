/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The equality `derived` (DOM) and `useDerived` (React) decide "changed" by: no serialization, no allocation. A function,
// a Set, a Map or a class instance is equal only to itself, so a changed handler is a change.

const isRecord = (v: unknown): v is Record<string, unknown> => {
  if (typeof v !== "object" || v === null || Array.isArray(v)) return false;
  const proto = Object.getPrototypeOf(v);
  return proto === Object.prototype || proto === null;
};

/** One level: the same value, or arrays of the same length or plain records of the same keys whose members are the same values. */
export function sameShallow(a: unknown, b: unknown): boolean {
  if (Object.is(a, b)) return true;
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!Object.is(a[i], b[i])) return false;
    return true;
  }
  if (!isRecord(a) || !isRecord(b)) return false;
  let n = 0;
  for (const k in a) {
    if (!Object.hasOwn(b, k) || !Object.is(a[k], b[k])) return false;
    n++;
  }
  for (const _ in b) n--;
  return n === 0;
}

/** Every level of plain arrays and records (a list of items, a label's lines); anything else by identity. */
export function sameData(a: unknown, b: unknown): boolean {
  if (Object.is(a, b)) return true;
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!sameData(a[i], b[i])) return false;
    return true;
  }
  if (!isRecord(a) || !isRecord(b)) return false;
  let n = 0;
  for (const k in a) {
    if (!Object.hasOwn(b, k) || !sameData(a[k], b[k])) return false;
    n++;
  }
  for (const _ in b) n--;
  return n === 0;
}
