/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of prepare.heap.ts (`test:perf`): LY2's reads, segment windows, `featureAt`, `ofPartAt` and a device
// column's (prepare-cases.ts, the guard's cases), each with its control.

import { beforeAll, describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { LENGTH, prepareReads } from './prepare-cases.js';
import { itemsIn } from './prepare-fixtures.js';

let c: Awaited<ReturnType<typeof prepareReads>>;

beforeAll(async () => {
  c = await prepareReads();
}, 120_000);

describe("LY2's reads, quick", () => {
  it('segment windows and the count in view', ({ task }) => {
    expect(itemsIn(c.genes, 0, 0, LENGTH, -1, 1e9)).toBeGreaterThan(100);
    quickCase(task, 'prepare: windows', c.windows);
  });

  it("the hit's featureAt and positionKey", ({ task }) => {
    quickCase(task, 'prepare: featureAt', c.features);
  });

  it('ofPartAt stepping and jumping', ({ task }) => {
    expect(c.of.count).toBeGreaterThan(10_000);
    quickCase(task, 'prepare: ofPartAt', c.ofParts);
  });

  it("a device column's reads", ({ task }) => {
    quickCase(task, 'prepare: columns', c.columns);
  });
});
