/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The translation over a model with an overlay that holds no node (_spec/README.md "Structures"): the empty overlay a
// viewer's store gives before its first frame knows no path of its own, so a pair of the model's systems is the
// model's path, and both ways read as over the model alone (EC-48).

import { locationsEqual, rangesLocation } from '@molstar/seqx-core/query';
import { EmptyOverlay, resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { fromUnits, toUnits } from '../units.js';
import { Complex, EntStart, p53Spec, Up } from './p53.js';

describe('an empty overlay', () => {
  const model = resolveSpec(p53Spec());
  const withEmpty = { model, overlay: EmptyOverlay };

  it('carries a P04637 hover onto the four units as the model alone does (EC-48)', () => {
    const loc = rangesLocation(Up, [EntStart + 10, EntStart + 13]);
    const alone = toUnits(model, loc).target;
    const over = toUnits(withEmpty, loc).target;
    expect(over.parts).toHaveLength(4);
    expect(locationsEqual(over, alone)).toBe(true);
  });

  it('places a unit’s residues back on P04637’s section paths as the model alone does', () => {
    const t = rangesLocation(`${Complex}/B`, [20, 30]);
    expect(locationsEqual(fromUnits(withEmpty, t).location, fromUnits(model, t).location)).toBe(true);
  });
});
