/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hydrophobicity scale lives with the residue schemes in core; this package re-exports it where it always was, the
// same table (one source of the numbers).

import {
  DefaultResidueHydrophobicity as CoreDefault,
  ResidueHydrophobicity as CoreScale,
} from '@molstar/core/util/color/residue-schemes';
import { describe, expect, it } from 'vitest';
import { DefaultResidueHydrophobicity, ResidueHydrophobicity } from '../types.js';

describe('the hydrophobicity scale', () => {
  it('is core’s table, re-exported here unchanged', () => {
    expect(ResidueHydrophobicity).toBe(CoreScale);
    expect(DefaultResidueHydrophobicity).toBe(CoreDefault);
    // Wimley–White: arginine's interface scale, a known value
    expect(ResidueHydrophobicity.ARG[0]).toBeCloseTo(0.81, 2);
    expect(Object.keys(ResidueHydrophobicity).length).toBeGreaterThanOrEqual(20);
  });
});
