/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The picker's text for the BLOSUM62 scheme says what the scheme does (_spec/README.md "Themes", Column schemes): it
// compares each letter with the reference row's, it does not score against the column's consensus.

import { describe, expect, it } from 'vitest';
import { schemeProvenance } from '../provenance.js';

describe('scheme provenance', () => {
  it('blosum62 says it compares each letter with the reference row', () => {
    const text = schemeProvenance('blosum62') ?? '';
    expect(text).toContain('reference row');
    expect(text).not.toContain('consensus');
  });
});
