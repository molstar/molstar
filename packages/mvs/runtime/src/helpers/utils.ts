/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from '@molstar/core/state';
import { range } from '@molstar/core/util/array';
import { Color } from '@molstar/core/util/color';
import { SplitColor as CoreSplitColor } from '@molstar/core/util/color/split-color';
import { decodeColor as _decodeColor } from '@molstar/core/util/color/utils';
import { isDefined } from '@molstar/core/util/object';
import { ColorT } from '@molstar/mvs-builder/tree/mvs/param-types';

/** Represents either the result or the reason of failure of an operation that might have failed */
export type Maybe<T> = { ok: true; value: T } | { ok: false; error: any };

/** Try to await a promise and return an object with its result (if resolved) or with the error (if rejected) */
export async function safePromise<T>(promise: T): Promise<Maybe<Awaited<T>>> {
  try {
    const value = await promise;
    return { ok: true, value };
  } catch (error) {
    return { ok: false, error };
  }
}

/** Return filtered array containing all original elements except `undefined` or `null`. */
export function filterDefined<T>(elements: (T | undefined | null)[]): T[] {
  return elements.filter((x) => x !== undefined && x !== null) as T[];
}

/** Create an 16-hex-character hash for a given input string, e.g. 'spanish inquisition' -> '7f9ac4be544330be'*/
export function stringHash(input: string): string {
  // = stringHash32(input) + stringHash32(reversed input), in one pass over the UTF-16 code units
  let h = 0;
  let r = 0;
  for (let i = 0, l = input.length; i < l; i++) {
    h = ((h << 5) - h + input.charCodeAt(i)) | 0;
    r = ((r << 5) - r + input.charCodeAt(l - 1 - i)) | 0;
  }
  return (h >>> 0).toString(16).padStart(8, '0') + (r >>> 0).toString(16).padStart(8, '0');
}

/** Return type of elements in a set */
export type ElementOfSet<S> = S extends Set<infer T> ? T : never;

/** Convert `colorString` (either X11 color name like 'magenta' or hex code like '#ff00ff') to Color.
 * If `colorString` contains a "split color" (e.g. 'red/white'), process only the first part.
 * Return `undefined` if `colorString` cannot be converted. */
export function decodeColor(colorString: string | number | undefined | null): Color | undefined {
  if (typeof colorString === 'number') {
    return Color(colorString);
  }
  if (colorString == null) return undefined;
  const slashIdx = colorString.indexOf('/'); // do not crash on split colors like 'red/white'
  return _decodeColor(slashIdx >= 0 ? colorString.substring(0, slashIdx) : colorString);
}

/** Utility functions for handling "split color" strings with primary and secondary color, e.g. 'red/white', '#ff0000/#ffffff'.
 * Regular color strings are treated as the same color repeated (e.g. 'red' == 'red/red'). */
export const SplitColor = {
  ...CoreSplitColor,
  /** Convert to hex code string, '#ff0000/#ffffff', '#ff0000' */
  toHexStyle(primaryColor: Color, secondaryColor?: Color): ColorT {
    return CoreSplitColor.toHexStyle(primaryColor, secondaryColor) as ColorT;
  },
};

export function collectMVSReferences<T extends StateObject.Ctor>(
  type: T[],
  dependencies: Record<string, StateObject>,
): Record<string, StateObject.From<T>['data']> {
  const ret: any = {};

  for (const key of Object.keys(dependencies)) {
    const o = dependencies[key];
    let okType = false;
    for (const t of type) {
      if (t.is(o)) {
        okType = true;
        break;
      }
    }
    if (!okType || !o.tags) continue;
    for (const tag of o.tags) {
      if (tag.startsWith('mvs-ref:')) {
        ret[tag.substring(8)] = o.data;
        break;
      }
    }
  }

  return ret;
}

export function getMVSReferenceObject<T extends StateObject.Ctor>(
  type: T[],
  dependencies: Record<string, StateObject> | undefined,
  ref: string,
): StateObject | undefined {
  if (!dependencies) return undefined;

  for (const key of Object.keys(dependencies)) {
    const o = dependencies[key];
    let okType = false;
    for (const t of type) {
      if (t.is(o)) {
        okType = true;
        break;
      }
    }
    if (!okType || !o.tags) continue;
    for (const tag of o.tags) {
      if (tag.startsWith('mvs-ref:')) {
        if (tag.substring(8) === ref) return o;
      }
    }
  }
}

/** Data structure for an array divided into contiguous groups */
export interface GroupedArray<T> {
  /** Number of groups */
  count: number;
  /** Get size of i-th group as `offsets[i+1]-offsets[i]`.
   * Get j-th element in i-th group as `grouped[offsets[i]+j]` */
  offsets: number[];
  /** Get j-th element in i-th group as `grouped[offsets[i]+j]` */
  grouped: T[];
}

export const GroupedArray = {
  getGroup<T>(groupedArray: GroupedArray<T>, iGroup: number): T[] {
    return groupedArray.grouped.slice(groupedArray.offsets[iGroup], groupedArray.offsets[iGroup + 1]);
  },
  /** Return element indices grouped by `group_by(element, index)`. Elements with `group_by(element, index)===undefined` are treated as separate groups. */
  groupIndices<T>(
    elements: readonly T[],
    group_by: (element: T, index: number) => string | undefined,
  ): GroupedArray<number> {
    let counter = 0;
    const groupMap = new Map<string, number>();
    const groups: number[] = [];
    for (let i = 0; i < elements.length; i++) {
      const groupId = group_by(elements[i], i);
      if (!isDefined(groupId)) {
        groups.push(counter++);
      } else {
        const groupIndex = groupMap.get(groupId);
        if (groupIndex === undefined) {
          groupMap.set(groupId, counter);
          groups.push(counter);
          counter++;
        } else {
          groups.push(groupIndex);
        }
      }
    }
    const elementIndices = range(elements.length).sort((i, j) => groups[i] - groups[j]);
    const offsets: number[] = [];
    for (let i = 0; i < elements.length; i++) {
      if (i === 0 || groups[elementIndices[i]] !== groups[elementIndices[i - 1]]) offsets.push(i);
    }
    offsets.push(elementIndices.length);
    return { count: offsets.length - 1, offsets, grouped: elementIndices };
  },
};
