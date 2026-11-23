/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A hover carried onto a structure's units (_spec/README.md "Structures", Performance): bytes a 2OCJ hover once its
// paths are kept, on one copy and on four, so the cost a lit unit reads as their difference. Not zero: the target is a
// new Location a hover (a host keeps it); bounded here so it does not grow. By hand, never in CI:
// `pnpm -C packages/seq/core test:perf`.

import { rangesLocation } from '@molstar/seqx-core/query';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, it } from 'vitest';
import { DefaultSyncLimits, toUnits } from '../units.js';
import { EntStart, p53Spec, Up } from './p53.js';

const hovers = Array.from({ length: 8 }, (_, k) => rangesLocation(Up, [EntStart + 20 + k, EntStart + 21 + k]));
const stepOf = (copies: number) => {
  const model = resolveSpec(p53Spec({ copies }));
  return (i: number) => void toUnits(model, hovers[i & 7], undefined, { limits: DefaultSyncLimits.hover });
};

describe('a hover onto the units', () => {
  it('2OCJ: bytes a hover on one copy and on four (a lit unit their difference over three)', ({ task }) => {
    const one = measureBytesPerStep(stepOf(1)).bytesPerStep;
    const four = stepOf(4);
    const { r } = quickCase(task, 'toUnits hover onto 2OCJ', four);
    const perUnit = (r.bytesPerStep - one) / 3;
    console.info(
      `toUnits, a 2OCJ hover: ${r.bytesPerStep.toFixed(0)} bytes on four units, ${one.toFixed(0)} on one; ${perUnit.toFixed(0)} a lit unit (${r.note})`,
    );
    // measured: 9,240 bytes on four units, 4,030 on one, 1,740 a lit unit (before the cheap cuts: 29,720, 9,380, 6,780)
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(10_500);
  });
});
