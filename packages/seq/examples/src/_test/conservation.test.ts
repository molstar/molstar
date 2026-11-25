/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultOpenLimits, type ValuesBins, type ValuesPositions } from '@molstar/seqx-core/data';
import { BigWigRuntime } from '@molstar/seqx-core/runtimes/bigwig';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets } from '../node-assets.js';
import { GeneLength } from '../p53/data.js';
import { ConservationFile, ConservationFormat, ConservationUncovered, ConservationZooms } from '../p53/conservation.js';
import { GeneOffset, RegionLength, RegionSystem } from '../p53/region.js';

// p53's conservation bigWig, read as a view reads it: by region and zoom level through the `bigwig` runtime, from the
// fixtures folder by the examples' Node adapter. It covers the gene's region (its system and length), and the gene's
// stretch is served at level 1 and at both zoom levels.

const signal = new AbortController().signal;

async function open() {
  const { reader } = nodeFixtureAssets();
  return BigWigRuntime.open(
    { address: '/data/conservation', format: ConservationFormat, uri: ConservationFile },
    { assets: reader, signal, limits: DefaultOpenLimits },
  );
}

describe('p53’s conservation over the gene’s region', () => {
  it('covers the region’s system with its length, at level 1 and its zoom levels', async () => {
    const s = await open();
    expect(s.coverage.systems).toEqual([{ name: RegionSystem, length: RegionLength }]);
    expect(s.coverage.levels).toEqual([1, ...ConservationZooms]);
    s.dispose();
  });

  it('the gene’s stretch: a value a base in [0, 1); bins whose means are in their bounds; the end without data', async () => {
    const s = await open();
    const region = {
      system: RegionSystem,
      start: GeneOffset,
      end: GeneOffset + GeneLength,
      shapes: ['values' as const],
    };
    const a = await s.read({ ...region, level: 1 }, signal);
    if (!a.ok) throw new Error(a.issue.message);
    const v = a.tables.values as ValuesPositions;
    expect(v.value.length).toBe(GeneLength);
    expect(v.value.every((x) => x >= 0 && x < 1)).toBe(true);
    for (const level of ConservationZooms) {
      const z = await s.read({ ...region, level }, signal);
      if (!z.ok) throw new Error(z.issue.message);
      const b = z.tables.values as ValuesBins;
      expect(b.start % level).toBe(0);
      for (let i = 0; i < b.count.length; i++) {
        expect(b.count[i]).toBeGreaterThan(0);
        expect(b.mean[i]).toBeGreaterThanOrEqual(b.min[i]);
        expect(b.mean[i]).toBeLessThanOrEqual(b.max[i] + 1e-6);
      }
    }
    const end = await s.read(
      {
        system: RegionSystem,
        start: RegionLength - ConservationUncovered,
        end: RegionLength,
        level: 1,
        shapes: ['values'],
      },
      signal,
    );
    if (!end.ok) throw new Error(end.issue.message);
    expect((end.tables.values as ValuesPositions).value.every(Number.isNaN)).toBe(true);
    s.dispose();
  });
});
