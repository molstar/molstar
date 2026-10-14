/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Equality of document values without building their text: two values are equal when their canonical JSON is
// (`canonical`, the patch's), so an edit compares as the stored document would. Identity first; a typed array
// equals a list of the same numbers (a Float32Array's as their shortest float32 decimals); null, a missing value and
// a number that is not finite (JSON's null) are equal; key order does not matter.
import { float32Number } from '@molstar/seq-core/model';

type List = ArrayLike<unknown>;

/** JSON writes each as null. */
const isNull = (x: unknown): boolean => x == null || (typeof x === 'number' && !Number.isFinite(x));

/** Whether `a` and `b` have equal canonical JSON, compared in place: no text, no copies. */
export function equal(a: unknown, b: unknown): boolean {
  if (a === b) return true;
  const na = isNull(a);
  const nb = isNull(b);
  if (na || nb) return na && nb;
  if (typeof a !== 'object' || typeof b !== 'object') return false;
  const va = ArrayBuffer.isView(a);
  const vb = ArrayBuffer.isView(b);
  const la = va || Array.isArray(a);
  const lb = vb || Array.isArray(b);
  if (la || lb) return la && lb && lists(a as List, b as List, va, vb);
  const ka = Object.keys(a as object);
  if (ka.length !== Object.keys(b as object).length) return false;
  for (const k of ka) {
    if (!Object.hasOwn(b as object, k)) return false;
    if (!equal((a as Record<string, unknown>)[k], (b as Record<string, unknown>)[k])) return false;
  }
  return true;
}

function lists(a: List, b: List, va: boolean, vb: boolean): boolean {
  const n = a.length;
  if (n !== b.length) return false;
  const fa = a instanceof Float32Array;
  const fb = b instanceof Float32Array;
  if (va && vb && fa === fb) {
    // Two typed arrays of one precision: equal numbers, or both written as null.
    const x = a as ArrayLike<number>;
    const y = b as ArrayLike<number>;
    for (let i = 0; i < n; i++) {
      const p = x[i];
      const q = y[i];
      if (p !== q && (Number.isFinite(p) || Number.isFinite(q))) return false;
    }
    return true;
  }
  for (let i = 0; i < n; i++) if (!equal(at(a, i, fa), at(b, i, fb))) return false;
  return true;
}

/** An element as its JSON reads it: a float32 as its shortest decimal. */
const at = (l: List, i: number, f32: boolean): unknown => {
  const x = l[i];
  return f32 && Number.isFinite(x as number) ? float32Number(x as number) : x;
};
