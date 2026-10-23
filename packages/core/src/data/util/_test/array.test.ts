/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { arrayPickIndices } from '../array.js';

describe('arrayPickIndices', () => {
  const indices = Int32Array.of(3, 0, 2, 2);
  const cases: [string, ArrayLike<any>][] = [
    ['Float32Array', Float32Array.of(0.5, 1.5, 2.5, 3.5)],
    ['Float64Array', Float64Array.of(0.1, 1.1, 2.1, 3.1)],
    ['Int32Array', Int32Array.of(-1, 10, 20, 30)],
    ['Uint8Array', Uint8Array.of(1, 2, 3, 255)],
    ['Int16Array', Int16Array.of(-5, 5, -50, 50)],
    ['Array', ['a', 'b', 'c', 'd']],
  ];

  it.each(cases)('gathers a %s into an array of its kind', (name, array) => {
    const picked = arrayPickIndices(array, indices);
    expect(picked.constructor).toBe(array.constructor);
    expect(Array.from(picked as ArrayLike<any>)).toEqual(Array.from(indices, (i) => array[i]));
  });

  it('backs Column.view of an array column of every kind (perf-audit-2 01/P10)', () => {
    for (const [, array] of cases) {
      const schema = typeof array[0] === 'string' ? Column.Schema.str : Column.Schema.float;
      const view = Column.view(Column.ofArray({ array, schema } as any), indices, false);
      expect(view.__array!.constructor).toBe(array.constructor);
      expect(Array.from({ length: indices.length }, (_, i) => view.value(i))).toEqual(
        Array.from(indices, (i) => array[i]),
      );
    }
  });
});
