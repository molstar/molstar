/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function shallowEqual<T extends object>(left: T, right: T): boolean {
  if (left === right) {
    return true;
  }
  const leftKeys = Reflect.ownKeys(left).filter((key) => Object.prototype.propertyIsEnumerable.call(left, key));
  const rightKeys = Reflect.ownKeys(right).filter((key) => Object.prototype.propertyIsEnumerable.call(right, key));
  if (leftKeys.length !== rightKeys.length) {
    return false;
  }
  for (const key of leftKeys) {
    if (!Object.hasOwn(right, key)) {
      return false;
    }
    if (!Object.is(Reflect.get(left, key), Reflect.get(right, key))) {
      return false;
    }
  }
  return true;
}
