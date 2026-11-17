/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BitFlags } from './bit-flags.js';
import { StringBuilder } from './string-builder.js';
import { UUID } from './uuid.js';
import { Mask } from './mask.js';
import { shallowEqualArrays, shallowEqualObjects } from './object.js';
import type { Progress } from '@molstar/core/task';

export * from './value-cell.js';
export { BitFlags, StringBuilder, UUID, Mask };

export const noop = function () {};

export function round(n: number, d: number) {
  const f = Math.pow(10, d);
  return Math.round(f * n) / f;
}

// One implementation each, in `object.ts` and `array.ts`; the names here are kept.
export { arrayEqual } from './array.js';
export { deepEqual, shallowEqualArrays, shallowEqualObjects } from './object.js';

/** Two arrays item by item, two objects by their own keys (`Object.is`); `null` and `undefined` differ. */
export function shallowEqual(a: any, b: any): boolean {
  if (Object.is(a, b)) return true;
  const arrA = Array.isArray(a);
  if (arrA !== Array.isArray(b)) return false;
  if (arrA) return shallowEqualArrays(a, b);
  if (a && b && typeof a === 'object' && typeof b === 'object') return shallowEqualObjects(a, b);
  return false;
}

const hasOwnProperty = Object.prototype.hasOwnProperty;

/** Returns `value` if not `undefined`, otherwise returns `defaultValue` */
export function defaults<T>(value: T | undefined, defaultValue: T): T {
  return value !== undefined ? value : defaultValue;
}

export function extend<S extends {}, T extends {}, U extends {}>(object: S, source: T, guard?: U): S & T & U {
  let v: any;

  const s = <any>source;
  const o = <any>object;
  const g = <any>guard;
  for (const k of Object.keys(source)) {
    v = s[k];
    if (v !== void 0) o[k] = v;
    else if (guard) o[k] = g[k];
  }

  if (guard) {
    for (const k of Object.keys(guard)) {
      v = o[k];
      if (v === void 0) o[k] = g[k];
    }
  }

  return <any>object;
}

export function shallowClone<T extends {}>(o: T): T {
  return extend({}, o) as T;
}

function _assign<T>(target: T): T {
  for (let s = 1; s < arguments.length; s++) {
    const from = arguments[s];
    for (const key of Object.keys(from)) {
      if (hasOwnProperty.call(from, key)) {
        (target as any)[key] = from[key];
      }
    }
  }
  return target;
}

export declare function _assignType<T>(o: T, ...from: any[]): T;
export const assign: <T>(o: T, ...from: any[]) => T = (Object as any).assign || _assign;

function _shallowMerge1<T extends {}>(source: T, update: T) {
  let changed = false;
  for (const k of Object.keys(update)) {
    if (!hasOwnProperty.call(update, k)) continue;

    if ((update as any)[k] !== (source as any)[k]) {
      changed = true;
      break;
    }
  }

  if (!changed) return source;
  return assign(shallowClone(source), update);
}

function _shallowMerge<T>(source: T) {
  let ret = source;

  for (let s = 1; s < arguments.length; s++) {
    if (!arguments[s]) continue;
    ret = _shallowMerge1(source, arguments[s]);
    if (ret !== source) {
      for (let i = s + 1; i < arguments.length; i++) {
        ret = assign(ret, arguments[i]);
      }
      break;
    }
  }
  return ret;
}

export const merge: <T>(source: T, ...rest: Partial<T>[]) => T = _shallowMerge;

function padTime(n: number) {
  return (n < 10 ? '0' : '') + n;
}
export function formatTime(d: Date) {
  const h = d.getHours(),
    m = d.getMinutes(),
    s = d.getSeconds();
  return `${h}:${padTime(m)}:${padTime(s)}`;
}

export function formatProgress(p: Progress) {
  const tp = p.root.progress;
  if (tp.isIndeterminate) return tp.message;
  const x = ((100 * tp.current) / tp.max).toFixed(2);
  return `${tp.message} ${x}%`;
}

export function formatBytes(count: number) {
  const units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];
  const i = Math.floor(Math.log(count) / Math.log(1024));
  return `${(count / Math.pow(1024, i)).toFixed(2)} ${units[i]}`;
}
