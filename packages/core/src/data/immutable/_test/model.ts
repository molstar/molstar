/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Shared fixtures of the model-based tests (README §8): a seeded generator and key pools with −0 and NaN.

import { PCG } from '../../util/hash-functions.js';

export class Rng {
  private pcg: PCG;
  constructor(seed: number) {
    this.pcg = new PCG(seed);
  }
  int(n: number) {
    return this.pcg.int() % n;
  }
  pick<T>(xs: readonly T[]): T {
    return xs[this.int(xs.length)];
  }
  shuffle<T>(xs: T[]): T[] {
    for (let i = xs.length - 1; i > 0; i--) {
      const j = this.int(i + 1);
      const t = xs[i];
      xs[i] = xs[j];
      xs[j] = t;
    }
    return xs;
  }
}

const Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_';

/** A 22-character id, as Mol*'s refs (`UUID.create22`). */
export function id22(rng: Rng) {
  let s = '';
  for (let i = 0; i < 22; i++) s += Alphabet[rng.int(64)];
  return s;
}

/** Keys of every kind the default hash covers, −0 and NaN included. */
export function keyPool(rng: Rng, n: number): unknown[] {
  const keys: unknown[] = [0, -0, NaN, 1, -1, 0.5, 1e300, null, undefined, true, false, '', 'a', 10n, Symbol('s')];
  const objects = [{}, {}, [], () => {}];
  keys.push(...objects);
  for (let i = 0; keys.length < n; i++) {
    switch (i % 3) {
      case 0:
        keys.push(i);
        break;
      case 1:
        keys.push(id22(rng));
        break;
      default:
        keys.push(`k${i}`);
    }
  }
  return keys;
}

/** Sorted by a stable string of each key (entries as `[k, v]`), −0 as 0: unordered contents compared. */
export function sorted<T>(xs: T[]): T[] {
  return xs.slice().sort((a, b) => (keyString(a) < keyString(b) ? -1 : keyString(a) > keyString(b) ? 1 : 0));
}

const objectIds = new WeakMap<object, number>();
let nextObjectId = 0;
function keyString(x: unknown): string {
  if (Array.isArray(x) && x.length === 2) return `${keyString(x[0])}=${keyString(x[1])}`;
  if (typeof x === 'object' && x !== null) {
    if (!objectIds.has(x)) objectIds.set(x, nextObjectId++);
    return `o${objectIds.get(x)}`;
  }
  if (typeof x === 'function') {
    if (!objectIds.has(x)) objectIds.set(x, nextObjectId++);
    return `f${objectIds.get(x)}`;
  }
  if (typeof x === 'symbol') return `s${x.toString()}`;
  // a JS `Map` keeps −0 as 0; a collection keeps the key first set
  if (Object.is(x, -0)) return 'number:0';
  return `${typeof x}:${String(x)}`;
}
