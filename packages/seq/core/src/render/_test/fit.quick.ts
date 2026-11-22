/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's level fitted to the store's bounds over blocks (_spec/README.md "Render", Data: the level rules): a layer
// laid out a block at a time (a grid's block pairs, a matrix over an assembly's blocks) keeps each block's fit, so a
// frame over 25 blocks fits nothing again and allocates nothing once warm. By hand, never in CI:
// `pnpm -C packages/seq/core test:perf`.

import type { Coverage } from '@molstar/seqx-core/data';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { fittedLevel, type LayerSource } from '../layer-data.js';
import type { StoreReader } from '../layer-type.js';

const Coverage: Coverage = { systems: [], shapes: [], levels: [1, 4, 16, 64] };

describe('a layer’s fit over 25 blocks', () => {
  it('each block pair’s fit kept: a frame over 25 × 25 blocks, 0 bytes warm', ({ task }) => {
    let asked = 0;
    const data: StoreReader = {
      coverage: () => Coverage,
      keysOf: () => [],
      get: () => undefined,
      settings: { tileBins: 1024, maxTilesPerWant: 64, gridTileBins: 256 },
      layerLoad: () => {
        asked++;
        return { fits: true };
      },
    };
    const srcs: LayerSource[] = Array.from({ length: 25 }, (_, i) => ({
      source: '/data/hic',
      system: `chr${i + 1}`,
      dataSystem: `/data/g/(sequence :name chr${i + 1})`,
      shape: 'grid',
      ratio: 1,
    }));
    const ys = srcs.map((s) => ({ system: s.system, start: 0, end: 4000 }));
    const regions = [0, 4000];
    const ctx = { key: '/views/v/plot/plot/contacts' };
    const frame = () => {
      for (const s of srcs) for (const y of ys) fittedLevel(ctx, data, s, regions, 4, y);
    };
    frame();
    const once = asked;
    frame();
    // the second frame asks the store nothing: every block pair's fit kept
    expect(asked).toBe(once);
    const { r } = quickCase(task, 'fit over 25 × 25 blocks', frame);
    console.info(`fit over 25 × 25 blocks: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
