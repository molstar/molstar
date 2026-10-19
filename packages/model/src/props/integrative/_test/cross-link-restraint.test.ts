/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Structure } from '@molstar/model/model/structure';
import { CrossLinkRestraint } from '../cross-link-restraint/property.js';

describe('cross-link restraints', () => {
  it('isLoci recognizes its own loci only', () => {
    const loci = CrossLinkRestraint.Loci(Structure.Empty, { pairs: [] } as any, [0]);
    expect(CrossLinkRestraint.isLoci(loci)).toBe(true);
    expect(CrossLinkRestraint.isLoci({ ...loci, tag: 'interactions' })).toBe(false);
  });
});
