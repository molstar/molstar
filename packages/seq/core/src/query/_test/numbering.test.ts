/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { compareLabels, hasLabelSuffix, isNoInsertionCode, readLabel } from '../labels.js';
import {
  findNumbering,
  labelIndex,
  nearestLabels,
  numberingIds,
  sequentialNumbering,
  snapLabel,
} from '../numbering.js';
import { Systems } from './fixture.js';

describe('query: labels', () => {
  it('order: integer, then no suffix, then letters (shorter first, case-blind), then .n (EC-2)', () => {
    const sorted = ['101', '100.2', '100AA', '100z', '100', '100A', '99', '100.1', '-1', 'H1', 'h0'].sort(
      compareLabels,
    );
    expect(sorted).toEqual(['-1', '99', '100', '100A', '100z', '100AA', '100.1', '100.2', '101', 'h0', 'H1']);
    // total: case breaks the last tie
    expect(Math.sign(compareLabels('100a', '100A'))).toBe(-Math.sign(compareLabels('100A', '100a')));
    expect(compareLabels('100a', '100A')).not.toBe(0);
  });

  it('reads integers and words; an integer outside int32 is refused, never wrapped (EC-4)', () => {
    expect(readLabel(94)).toEqual({ label: '94' });
    expect(readLabel('94A')).toEqual({ label: '94A' });
    expect(readLabel(-5)).toEqual({ label: '-5' });
    expect(readLabel(2 ** 31 - 1)).toEqual({ label: '2147483647' });
    expect(readLabel(2 ** 31)).toEqual({ error: 'int32' });
    expect(readLabel('99999999999A')).toEqual({ error: 'int32' });
    for (const bad of [1.5, '', ' 9', true, null, [1]]) expect(readLabel(bad)).toEqual({ error: 'syntax' });
  });

  it('knows the spellings of no insertion code, and which labels have a suffix (EC-5)', () => {
    expect(['', '?', '.', 'A'].map(isNoInsertionCode)).toEqual([true, true, true, false]);
    expect(['94', '94A', '100.1', 'H1'].map(hasLabelSuffix)).toEqual([false, true, true, false]);
  });
});

describe('query: numberings', () => {
  it('index is 0-based and ordinal 1-based, over the system only (EC-1)', () => {
    const idx = sequentialNumbering('index', 5, 0);
    const ord = sequentialNumbering('ordinal', 5, 1);
    expect([idx.labelAt(0), ord.labelAt(0), ord.labelAt(5)]).toEqual(['0', '1', undefined]);
    expect([idx.indexOf('4'), idx.indexOf('5'), ord.indexOf('1'), ord.indexOf('0')]).toEqual([4, -1, 0, -1]);
    // one spelling per position
    expect([ord.indexOf('01'), ord.indexOf('+1'), idx.indexOf('-0'), ord.indexOf('1A')]).toEqual([-1, -1, -1, -1]);
  });

  it('a repeated label names its first position (EC-4)', () => {
    const n = findNumbering(Systems.p53, 'dup')!;
    expect(labelIndex(n, Systems.p53, '2')).toBe(1);
  });

  it("a unit's own numbering beats its molecule's of the same id; else the molecule's, shifted and clipped (EC-8)", () => {
    const unit = Systems.unit;
    const auth = findNumbering(unit, 'auth')!;
    expect(auth.shift).toBe(0);
    expect(labelIndex(auth, unit, '1')).toBe(0);
    const kabat = findNumbering(unit, 'kabat')!;
    expect(kabat.shift).toBe(-2);
    // the molecule's 82 is at 2, the unit's first position
    expect(labelIndex(kabat, unit, '82')).toBe(0);
    // the molecule's 80 is before the unit, 88 past it
    expect([labelIndex(kabat, unit, '80'), labelIndex(kabat, unit, '88')]).toEqual([-1, -1]);
    expect(findNumbering(unit, 'nope')).toBeUndefined();
    expect(numberingIds(unit).sort()).toEqual(['auth', 'imgt', 'index', 'kabat', 'ordinal']);
  });

  it('snaps a missing range end inward by the order of labels (EC-6)', () => {
    const ab = Systems.ab;
    const kabat = findNumbering(ab, 'kabat')!;
    // 84 is a hole: a start snaps to 85, an end to 83
    expect([snapLabel(kabat, ab, '84', 'start'), snapLabel(kabat, ab, '84', 'end')]).toEqual([8, 6]);
    // 82D is past 82C and before 83
    expect([snapLabel(kabat, ab, '82D', 'start'), snapLabel(kabat, ab, '82D', 'end')]).toEqual([6, 5]);
    expect([snapLabel(kabat, ab, '1', 'end'), snapLabel(kabat, ab, '1000', 'start')]).toEqual([-1, -1]);
    expect(nearestLabels(kabat, ab, '84', 4)).toEqual(['82C', '83', '85', '86']);
  });
});
