/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { mapsEqual } from '../map.js';

describe('mapsEqual', () => {
  it('the same keys with the same values, whatever the order; a missing key, another value or size is not', () => {
    const a = new Map([
      ['x', 1],
      ['y', 2],
    ]);
    expect(mapsEqual(a, a)).toBe(true);
    expect(
      mapsEqual(
        a,
        new Map([
          ['y', 2],
          ['x', 1],
        ]),
      ),
    ).toBe(true);
    expect(mapsEqual(a, new Map([['x', 1]]))).toBe(false);
    expect(
      mapsEqual(
        a,
        new Map([
          ['x', 1],
          ['y', 3],
        ]),
      ),
    ).toBe(false);
    expect(
      mapsEqual(
        new Map<string, number | undefined>([['x', undefined]]),
        new Map<string, number | undefined>([['z', undefined]]),
      ),
    ).toBe(false);
    expect(mapsEqual(new Map([['n', Number.NaN]]), new Map([['n', Number.NaN]]))).toBe(true);
  });
});
