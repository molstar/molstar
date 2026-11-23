/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A whole genome's layout at the fit (_spec/README.md "Render", Blocks): 25 chromosomes × 3
// layers, and 3,400 sequences (25 chromosomes and 3,375 contigs, GRCh38 with its alternates and decoys) in runs, each a
// frame laid out again (its width changed by a pixel) under its bound. Each bound is the heavy lane's (`itWithin`, best
// of three); the quick lane runs them unbounded.

import type { Node } from '@molstar/core/repr/tree';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { itWithin } from '@molstar/test-harness/timing';
import { beforeAll, describe, expect } from 'vitest';
import type { BlockPlacement, RingPlacement } from '../index.js';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';

beforeAll(async () => {
  await loadGenome();
});

const Frames = 20;

/** The median ms of a layout at the fit, its width changed by a pixel each frame (so each frame lays it out again). */
async function fitTime(spec: Node, height: number) {
  const s = await settle(spec);
  const input = (k: number) => ({ size: { width: 1200 + (k & 1), height }, dpr: 2 });
  await s.frame(input(0));
  await s.frame(input(1));
  const time = () => {
    const t: number[] = [];
    for (let k = 0; k < Frames; k++) {
      const t0 = performance.now();
      const f = s.renderer.frame(input(k));
      t.push(performance.now() - t0);
      if (!f.rebuilt) throw new Error('a frame at another width lays the genome out again');
    }
    t.sort((a, b) => a - b);
    return t[t.length >> 1];
  };
  return { s, time, last: () => s.renderer.frame(input(0)) };
}

describe('a whole genome at the fit', () => {
  itWithin('25 chromosomes × 3 layers: a layout and a draw', async (within) => {
    const { time, last } = await fitTime(genomeSpec(toyGenome(25, 0), { tracks: 3 }), 400);
    time();
    const median = within(4 * Frames, time);
    expect((last().layout.sections[0].placement as BlockPlacement).strips.n).toBe(25);
    console.info(`whole genome, 25 blocks × 3 layers: median ${median.toFixed(3)} ms`);
  });

  itWithin('3,400 sequences in runs: a layout and a draw', async (within) => {
    const { time, last } = await fitTime(genomeSpec(toyGenome(25, 3375)), 300);
    time();
    const median = within(6 * Frames, time);
    const p = last().layout.sections[0].placement as BlockPlacement;
    expect(p.runN).toBeGreaterThan(0);
    // the joint rule's read set: the blocks a device pixel wide or more (a strip each), the rest the runs' bands
    expect(p.strips.n).toBeLessThan(400);
    console.info(
      `whole genome, 3,400 blocks in runs: median ${median.toFixed(3)} ms; ${p.strips.n} blocks read, ${p.runN} runs`,
    );
  });

  itWithin('3,400 sequences on a ring, in runs: a layout and a draw', async (within) => {
    const { time, last } = await fitTime(genomeSpec(toyGenome(25, 3375), { layout: 'rings' }), 1200);
    time();
    const median = within(12 * Frames, time);
    expect((last().layout.sections[0].placement as RingPlacement).runN).toBeGreaterThan(0);
    console.info(`whole genome ring, 3,400 blocks in runs: median ${median.toFixed(3)} ms`);
  });
});
