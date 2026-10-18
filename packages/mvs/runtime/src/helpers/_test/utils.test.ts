/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { range } from '@molstar/core/util/array';
import { MVSAnnotationRow } from '../schemas.js';
import { hashString } from '@molstar/core/data/util';
import { GroupedArray, stringHash } from '../utils.js';

describe('GroupedArray', () => {
  it('GroupedArray.groupIndices', async () => {
    const rows = [
      { label: 'A' },
      { label: 'B', group_id: 1 },
      { label: 'C', group_id: 'x' },
      { label: 'D', group_id: 1 },
      { label: 'E' },
      { label: 'F' },
      { label: 'G', group_id: 'x' },
      { label: 'H', group_id: 'x' },
    ] as any as MVSAnnotationRow[];
    const g = GroupedArray.groupIndices(rows, (row) => row.group_id);
    const groupedIndices = range(g.count).map((i) => g.grouped.slice(g.offsets[i], g.offsets[i + 1]));
    const groupedRows = groupedIndices.map((group) => group.map((j) => rows[j]));
    expect(groupedRows).toEqual([
      [{ label: 'A' }],
      [
        { label: 'B', group_id: 1 },
        { label: 'D', group_id: 1 },
      ],
      [
        { label: 'C', group_id: 'x' },
        { label: 'G', group_id: 'x' },
        { label: 'H', group_id: 'x' },
      ],
      [{ label: 'E' }],
      [{ label: 'F' }],
    ]);
  });
});

describe('stringHash', () => {
  // the earlier implementation: refs made from it are in saved sessions, so the output must not change
  const hash32 = (s: string) => (hashString(s) >>> 0).toString(16).padStart(8, '0');
  const reference = (s: string) => hash32(s) + hash32(s.split('').reverse().join(''));

  it('is the hash of the string and of its reverse', () => {
    for (const s of [
      '',
      'a',
      'spanish inquisition',
      '-=root=-:ms-plugin.download',
      'ünïcödé 😀 surrogates',
      'x'.repeat(1000),
    ]) {
      expect(stringHash(s)).toBe(reference(s));
    }
    expect(stringHash('spanish inquisition')).toBe('7f9ac4be544330be');
  });
});
