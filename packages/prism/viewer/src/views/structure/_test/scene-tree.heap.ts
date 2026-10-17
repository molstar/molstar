// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The scene tree's hot paths keep nothing (`test:heap`): what stays after a collection over n and 4n steps is flat;
// a step that keeps what it made is caught.

import { describe, expect, it } from 'vitest';
import { describeScaling, isRetainedFlat, retainedScaling } from '../../../_test/perf.js';
import { jump, pick, scrollStep } from './scene-tree-cases.js';

const warm = (step: (i: number) => void, n: number) => {
  for (let i = 0; i < n; i++) step(i);
};

describe('scene tree heap guards', () => {
  for (const [name, make, n] of [
    ['a pick', pick, 2_000],
    ['a scroll step', scrollStep, 2_000],
    ['a jump', jump, 300],
  ] as const) {
    it(name, () => {
      const { step } = make();
      warm(step, n);
      const r = retainedScaling(n, step);
      expect(isRetainedFlat(r, n), describeScaling(r, n)).toBe(true);
    });
  }

  it('a jump that keeps what it made is caught (positive control)', () => {
    const { step } = jump();
    warm(step, 300);
    const kept: unknown[] = [];
    expect(
      isRetainedFlat(
        retainedScaling(300, (i) => {
          step(i);
          kept.push(new Uint8Array(1 << 10));
        }),
        300,
      ),
    ).toBe(false);
  });
});
