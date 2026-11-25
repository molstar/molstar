/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { ElementSymbol } from '../../../../types.js';
import { AtomicNumbers, AtomNumber, DefaultVdwRadius, VdwRadius, VdwRadiusOfNumber } from '../measures.js';

describe('VdwRadiusOfNumber', () => {
  it('is VdwRadius of the symbol for every element and for unknown ones', () => {
    const symbols = [...Object.keys(AtomicNumbers), 'X', '', 'Q', 'UNK'] as ElementSymbol[];
    for (const s of symbols) {
      expect(VdwRadiusOfNumber(AtomNumber(s))).toBe(VdwRadius(s));
    }
    expect(VdwRadiusOfNumber(0)).toBe(DefaultVdwRadius);
  });

  it('gives elements without a radius (110 to 118) the default one', () => {
    expect(VdwRadius('OG' as ElementSymbol)).toBe(DefaultVdwRadius);
    expect(VdwRadiusOfNumber(118)).toBe(DefaultVdwRadius);
  });
});
