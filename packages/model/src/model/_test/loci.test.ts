/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Loci } from '../loci.js';

describe('Loci.GranularityOptions', () => {
  it('lists only the instance granularities under With Symmetry (perf-audit-2 04/B6)', () => {
    // an option is [value, label] or [value, label, group]
    const groups = new Map(Loci.GranularityOptions.map((o) => [o[0], (o as string[])[2]]));
    for (const k of ['residue', 'chain', 'entity', 'model', 'operator', 'element', 'structure'] as const) {
      expect(groups.get(k), k).toBeUndefined();
    }
    for (const k of ['residueInstances', 'chainInstances', 'elementInstances'] as const) {
      expect(groups.get(k), k).toBe('With Symmetry');
    }
  });
});
