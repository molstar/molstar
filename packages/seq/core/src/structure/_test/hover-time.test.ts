/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure translation's representatives (_spec/README.md "Structures", Performance): one 2OCJ hover carried
// onto its four units a frame and a selection on every one of 32,768 copies (quick, bounded in the heavy lane: linear
// in the parts and the units reached, one search a part and molecule), and a hover on an assembly of `maxUnits` units
// gated to whole units (EC-167, heavy: the scale).

import { positionCount, rangesLocation } from '@molstar/seqx-core/query';
import { DefaultResolveConfig, resolveSpec } from '@molstar/seqx-core/resolve';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { DefaultSyncLimits, toUnits } from '../units.js';
import { EntStart, Observed, p53Spec, Up } from './p53.js';

const hover = (k: number) => rangesLocation(Up, [EntStart + 20 + (k % 150), EntStart + 21 + (k % 150)]);

describe('a hover a frame', () => {
  const model = resolveSpec(p53Spec());

  itWithin('one 2OCJ hover carried onto its four units (EC-48)', (within) => {
    toUnits(model, hover(0), undefined, { limits: DefaultSyncLimits.hover }); // warm: the paths found once
    // measured: 100 hovers in 1.6 ms (an M-series laptop); 20 ms leaves room for a parallel heavy run
    const lit = within(20, () => {
      let n = 0;
      for (let k = 0; k < 100; k++)
        n += toUnits(model, hover(k), undefined, { limits: DefaultSyncLimits.hover }).target.parts.length;
      return n;
    });
    expect(lit).toBe(400);
  });
});

describe('a selection on every copy', () => {
  const copies = 32_768;

  itWithin(
    'a selection on 32,768 copies: past the gate, whole units; ungated, a run on each (EC-167)',
    (within) => {
      // made in the body: a try again searches again
      const model = resolveSpec(p53Spec({ copies, auth: false, zinc: false }));
      const sel = rangesLocation(Up, [EntStart + 10, EntStart + 60]);
      // measured: about 70 ms gated (the first: its searches) and 50 ms exact (an M-series laptop); the gate stops counting once over
      const gated = within(400, () => toUnits(model, sel));
      expect([gated.gated, gated.target.parts.length]).toEqual(['units', copies]);
      const exact = within(400, () => toUnits(model, sel, undefined, { limits: null }));
      expect(positionCount(exact.target)).toBe(copies * 50);
      expect(positionCount(gated.target)).toBe(copies * (Observed[1] - Observed[0]));
    },
    30_000,
  );
});

describe('an assembly of maxUnits units', { tags: ['heavy'] }, () => {
  // half the bound in copies of the domain, half in their zincs
  const copies = DefaultResolveConfig.maxUnits / 2;

  itWithin(
    'a hover on every copy, gated to whole units (EC-167)',
    (within) => {
      // made in the body: a try again searches again, and no other lane builds it
      const model = resolveSpec(p53Spec({ copies, auth: false }));
      // measured: about 60 ms each (an M-series laptop: one search a part and molecule, then each unit's step)
      const first = within(1_000, () => toUnits(model, hover(0), undefined, { limits: DefaultSyncLimits.hover }));
      expect([first.gated, first.target.parts.length]).toEqual(['units', copies]);
      // the paths kept: a hover a frame after the first
      const next = within(500, () => toUnits(model, hover(1), undefined, { limits: DefaultSyncLimits.hover }));
      expect(next.target.parts.length).toBe(copies);
    },
    60_000,
  );
});
