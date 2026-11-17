/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  elements,
  elementsEqual,
  hashElements,
  intersectElements,
  subtractElements,
  unionElements,
} from '../elements.js';
import {
  EmptyLocation,
  hashLocation,
  intersectLocations,
  isEmptyLocation,
  location,
  locationsEqual,
  normalizeRanges,
  overlapsRanges,
  positionCount,
  rangesLocation,
  subtractLocations,
  unionLocations,
} from '../location.js';

const ranges = (l: { parts: readonly { system: string; ranges: Int32Array }[] }) =>
  l.parts.map((p) => [p.system, [...p.ranges]]);

describe('query: Location', () => {
  it('holds sorted, merged, half-open ranges per system, one part each, empty ones dropped', () => {
    const l = location({
      parts: [
        { system: '/b', ranges: [10, 20, 0, 5, 5, 8, 15, 30, 40, 40] },
        { system: '/a', ranges: [3, 4] },
        { system: '/b', ranges: [50, 60] },
        { system: '/c', ranges: [7, 7] },
      ],
    });
    expect(ranges(l)).toEqual([
      ['/a', [3, 4]],
      ['/b', [0, 8, 10, 30, 50, 60]],
    ]);
    expect(l.parts[1].ranges).toBeInstanceOf(Int32Array);
    expect(positionCount(l)).toBe(1 + 8 + 20 + 10);
    expect(location({ parts: [{ system: '/c', ranges: [7, 7] }] })).toBe(EmptyLocation);
    expect(() => normalizeRanges([1])).toThrow();
    expect(() => normalizeRanges([0.5, 2])).toThrow();
  });

  it('is a value: frozen, equal and hashed alike however it was built', () => {
    const a = location({
      parts: [
        { system: '/s', ranges: [5, 9, 0, 2] },
        { system: '/t', ranges: [1, 2] },
      ],
    });
    const b = location({
      parts: [
        { system: '/t', ranges: [1, 2] },
        { system: '/s', ranges: [0, 1, 1, 2, 5, 9] },
      ],
    });
    expect(Object.isFrozen(a) && Object.isFrozen(a.parts) && Object.isFrozen(a.parts[0])).toBe(true);
    expect(locationsEqual(a, b)).toBe(true);
    expect(hashLocation(a)).toBe(hashLocation(b));
    const c = rangesLocation('/s', [0, 2, 5, 10]);
    expect(locationsEqual(a, c)).toBe(false);
    expect(hashLocation(c)).not.toBe(hashLocation(a));
  });

  it('keeps rects sorted and each once; a rect with an empty side is none', () => {
    const r1 = { a: { system: '/x', ranges: [10, 25] }, b: { system: '/x', ranges: [40, 55] } };
    const r0 = { a: { system: '/x', ranges: [0, 5] }, b: { system: '/y', ranges: [0, 5] } };
    const l = location({ rects: [r1, r0, r1, { a: r0.a, b: { system: '/y', ranges: [] } }] });
    expect(l.rects.map((r) => [r.a.system, [...r.a.ranges], r.b.system, [...r.b.ranges]])).toEqual([
      ['/x', [0, 5], '/y', [0, 5]],
      ['/x', [10, 25], '/x', [40, 55]],
    ]);
    expect(locationsEqual(l, location({ rects: [r0, r1] }))).toBe(true);
    expect(hashLocation(l)).toBe(hashLocation(location({ rects: [r0, r1] })));
  });

  it('unions, intersects and subtracts per system', () => {
    const a = location({
      parts: [
        { system: '/s', ranges: [0, 10, 20, 30] },
        { system: '/t', ranges: [0, 5] },
      ],
    });
    const b = location({
      parts: [
        { system: '/s', ranges: [5, 25] },
        { system: '/u', ranges: [1, 2] },
      ],
    });
    expect(ranges(unionLocations(a, b))).toEqual([
      ['/s', [0, 30]],
      ['/t', [0, 5]],
      ['/u', [1, 2]],
    ]);
    expect(ranges(intersectLocations(a, b))).toEqual([['/s', [5, 10, 20, 25]]]);
    expect(ranges(subtractLocations(a, b))).toEqual([
      ['/s', [0, 5, 25, 30]],
      ['/t', [0, 5]],
    ]);
    expect(isEmptyLocation(subtractLocations(a, a))).toBe(true);
  });

  it('intersects rects over the same two systems side by side', () => {
    const rect = (a: number[], b: number[]) => ({ a: { system: '/x', ranges: a }, b: { system: '/y', ranges: b } });
    const p = location({ rects: [rect([0, 10], [0, 10])] });
    const q = location({ rects: [rect([5, 20], [8, 9]), rect([30, 40], [0, 10])] });
    const i = intersectLocations(p, q);
    expect(i.rects.map((r) => [[...r.a.ranges], [...r.b.ranges]])).toEqual([
      [
        [5, 10],
        [8, 9],
      ],
    ]);
    expect(unionLocations(p, q).rects.length).toBe(3);
  });

  it('tells overlap of an item with ranges; a zero-length item where its point is inside', () => {
    const r = normalizeRanges([10, 20, 30, 40]);
    expect([overlapsRanges(r, 0, 10), overlapsRanges(r, 0, 11), overlapsRanges(r, 20, 30)]).toEqual([
      false,
      true,
      false,
    ]);
    expect([overlapsRanges(r, 35, 50), overlapsRanges(r, 10, 10), overlapsRanges(r, 20, 20)]).toEqual([
      true,
      true,
      false,
    ]);
  });
});

describe('query: Elements', () => {
  const k = (owner: string, key: string | number) => ({ owner, key });

  it('are sorted and each once; a number key never equals its string', () => {
    const e = elements([k('/b', 'x'), k('/a', 'y'), k('/a', 2), k('/a', '2'), k('/b', 'x'), k('/a', 10)]);
    expect(e.keys).toEqual([k('/a', 2), k('/a', 10), k('/a', '2'), k('/a', 'y'), k('/b', 'x')]);
    expect(Object.isFrozen(e.keys)).toBe(true);
  });

  it('compare and hash as values, and combine as sets', () => {
    const a = elements([k('/a', 1), k('/a', 'x')]);
    const b = elements([k('/a', 'x'), k('/a', 1)]);
    expect(elementsEqual(a, b) && hashElements(a) === hashElements(b)).toBe(true);
    expect(hashElements(elements([k('/a', 1)]))).not.toBe(hashElements(elements([k('/a', '1')])));
    const c = elements([k('/a', 'x'), k('/b', 'z')]);
    expect(unionElements(a, c).keys).toEqual([k('/a', 1), k('/a', 'x'), k('/b', 'z')]);
    expect(intersectElements(a, c).keys).toEqual([k('/a', 'x')]);
    expect(subtractElements(a, c).keys).toEqual([k('/a', 1)]);
  });
});
