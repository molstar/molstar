/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import type { CustomProperty } from '../../../common/custom-property.js';
import type { Volume } from '@molstar/model/model/volume/volume';
import type { Grid } from '@molstar/model/model/volume/grid';
import { calculateBasicStreamlines, type BasicStreamlineCalculationProps } from '../basic.js';
import { Orders, grid, hash, viaAccessors } from '../../../../model/volume/_test/grid-test-util.js';

async function streamlines(g: Grid, props: BasicStreamlineCalculationProps) {
  const ctx = { runtime: SyncRuntimeContext } as CustomProperty.Context;
  const lines = await calculateBasicStreamlines(ctx, { grid: g } as Volume, props);
  const points: number[] = [];
  for (const line of lines) for (const p of line) points.push(p[0], p[1], p[2]);
  return { lines: lines.length, points: points.length / 3, hash: hash(points) };
}

async function sameForEveryOrder(props: BasicStreamlineCalculationProps) {
  const reference = await streamlines(grid([0, 1, 2]), props);
  expect(reference.lines).toBeGreaterThan(0);
  for (const order of Orders) {
    expect(await streamlines(grid(order), props), `order ${order}`).toEqual(reference);
    expect(await streamlines(viaAccessors(grid(order)), props), `order ${order}, accessors`).toEqual(reference);
  }
  return reference;
}

describe('basic streamlines', () => {
  it('are the same for every axis order', async () => {
    expect(await sameForEveryOrder({ seedDensity: 6, stepSize: 0.35, minSeparation: 0 })).toMatchInlineSnapshot(`
      {
        "hash": "b82e366c",
        "lines": 441,
        "points": 10467,
      }
    `);
  });

  it('are the same for every axis order, with a minimum separation', async () => {
    expect(await sameForEveryOrder({ seedDensity: 6, stepSize: 0.35, minSeparation: 1.5 })).toMatchInlineSnapshot(`
      {
        "hash": "c4ea8817",
        "lines": 14,
        "points": 182,
      }
    `);
  });
});
