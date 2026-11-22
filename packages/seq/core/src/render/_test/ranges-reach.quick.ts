/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's reachable read (`noteReachable`, _spec/README.md "Render", Track heights) over a genome (83 Mb): where the
// fit level holds no items, or the read is past the layer's bounds (said), the verdict is kept for its fit and a warm
// layout allocates nothing. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import type { LayoutContext, StoreReader } from '../layer-type.js';
import type { LayerSource } from '../layer-data.js';
import type { RangesPrep } from '../ranges-lay.js';
import { noteReachable } from '../ranges-reach.js';

function warm(levels: number[], fits: boolean) {
  const coverage = { levels };
  const pieces = [0, 83_000_000];
  const issues: unknown[] = [];
  const ctx = {
    axis: { pieces, dpr: 2 },
    box: { width: 1000 },
    key: '/views/v/s/t:0',
    settled: true,
    report: (i: unknown) => issues.push(i),
  } as unknown as LayoutContext;
  const data = {
    coverage: () => coverage,
    layerLoad: () => ({ fits }),
    layerRefused: () => undefined,
  } as unknown as StoreReader;
  const src = { source: '/data/s', system: 'chr', ratio: 1 } as unknown as LayerSource;
  const prep = { arrow: 0 } as unknown as RangesPrep;
  return { run: () => noteReachable(prep, src, ctx, data), issues };
}

describe('a source’s reachable read over a genome: a kept verdict, nothing allocated warm', () => {
  for (const [what, levels, fits] of [
    ['a coarse fit level (a density)', [1, 1024], true],
    ['a read past the layer’s bounds (said)', [1], false],
  ] as const)
    it(`${what}: 0 bytes warm`, async ({ task }) => {
      const w = warm([...levels], fits);
      w.run();
      const { r } = quickCase(task, `reach verdict: ${what}`, () => w.run(), { control: false, why: 'a pooled issue' });
      console.info(`reach verdict, ${what}: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      if (!fits) expect(w.issues.length).toBeGreaterThan(0);
    });
});
