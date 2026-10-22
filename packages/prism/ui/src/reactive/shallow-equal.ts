/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const hasOwn = Object.prototype.hasOwnProperty;
const isEnumerable = Object.prototype.propertyIsEnumerable;

/**
 * Own enumerable keys, strings and symbols, equal by `Object.is`. The default compare of every record state
 * (`sameShallow`): string keys by `for…in` and a count, no key arrays made (perf audit 2, 11/P5: 925 → 120 ns at 12 keys).
 */
export function shallowEqual<T extends object>(left: T, right: T): boolean {
  if (left === right) {
    return true;
  }
  let count = 0;
  for (const key in left) {
    if (!hasOwn.call(left, key)) continue;
    if (!Object.is(left[key], right[key]) || !hasOwn.call(right, key)) return false;
    count++;
  }
  for (const key in right) {
    if (hasOwn.call(right, key)) count--;
  }
  return count === 0 && sameSymbols(left, right);
}

function sameSymbols(left: object, right: object): boolean {
  const leftSymbols = Object.getOwnPropertySymbols(left);
  const rightSymbols = Object.getOwnPropertySymbols(right);
  if (leftSymbols.length === 0 && rightSymbols.length === 0) return true;
  let count = 0;
  for (const key of leftSymbols) {
    if (!isEnumerable.call(left, key)) continue;
    if (!isEnumerable.call(right, key) || !Object.is(Reflect.get(left, key), Reflect.get(right, key))) return false;
    count++;
  }
  for (const key of rightSymbols) {
    if (isEnumerable.call(right, key)) count--;
  }
  return count === 0;
}
