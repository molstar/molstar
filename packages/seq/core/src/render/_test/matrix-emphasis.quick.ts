/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The matrix layer's emphasis redraws (_spec/README.md "Render", The matrix layer, Emphasis): the selection of its
// states, own and projected, and the hover of a member, own hit and projected, allocate nothing once warm. By hand,
// never in CI: `pnpm -C packages/seq/core test:perf`.

import { rangesLocation } from '@molstar/seqx-core/query';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import type { Emphasis, FrameInput } from '../index.js';
import { settle } from './helpers.js';
import { matrixSpec, States } from './matrix-fixtures.js';

const Input: FrameInput = { size: { width: 800, height: 320 }, dpr: 2, settled: true };

const Cases: [string, Emphasis][] = [
  ['own selection', { selection: rangesLocation(States, [1, 3, 5, 6]) }],
  ['projected selection', { projected: { selection: rangesLocation(States, [1, 3, 5, 6]) } }],
  ['a hover on another layer’s member', { hover: { layer: 'o', section: 'x', system: States, position: 3 } }],
  ['a projected hover', { projected: { hover: rangesLocation(States, [2, 4]), hoverSection: 'x' } }],
];

describe('the matrix layer’s emphasis redraws with nothing allocated once warm', () => {
  for (const [name, emphasis] of Cases)
    it(name, async ({ task }) => {
      const s = await settle(matrixSpec());
      const at = { ...Input, emphasis };
      for (let i = 0; i < 4; i++) await s.frame(at);
      const { r } = quickCase(task, `matrix ${name}`, () => void s.renderer.frame(at), {});
      console.info(`matrix ${name}: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    });
});
