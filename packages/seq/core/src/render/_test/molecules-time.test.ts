/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A complex's frame target: 200 copies of a chain merged into one
// block (EC-67: one block's cost, not 200), its glycan, ATP and an ion as chips, the modified residues marked and the
// links as arcs; a pan of the merged block and a zoom about it; the 200 copies split, pairs within the chain's molecule
// on each copy (EC-48, EC-65), a pan along them. The bound is the heavy lane's (`itWithin`, best of
// three); the quick lane runs it unbounded, a representative of the frame target.

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import type { BlockPlacement, FrameInput, PairsLay } from '../index.js';
import { settle } from './helpers.js';
import { moleculesSpec } from './molecules-fixtures.js';

const Size = { width: 1440, height: 400 };
const Frames = 40;

describe('a complex’s frame time', () => {
  itWithin('200 merged copies and a glycan chip: a pan and a zoom frame (EC-67)', async (within) => {
    const s = await settle(moleculesSpec({ chains: 200 }));
    const inputs: FrameInput[] = [];
    for (let k = 0; k < 8; k++) {
      const start = 10 + k * 3;
      inputs.push({ size: Size, dpr: 2, viewport: { start, end: start + 60 + (k % 2) * 20 } });
    }
    for (const i of inputs) await s.frame(i);
    const f = s.renderer.frame(inputs[0]);
    const p = f.layout.sections[0].placement as BlockPlacement;
    // one block for the 200 copies: its strip and the chips' after it
    expect(f.layout.sections[0].section.blocks![0].members.length).toBe(200);
    expect(p.strips.n).toBeLessThanOrEqual(4);
    const time = () => {
      const t: number[] = [];
      for (let k = 0; k < Frames; k++) {
        const t0 = performance.now();
        s.renderer.frame(inputs[(k + 1) % inputs.length]);
        t.push(performance.now() - t0);
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    time();
    const median = within(8 * Frames, time);
    console.info(`200 copies and a glycan chip, a pan or zoom: median ${median.toFixed(3)} ms`);
  });

  itWithin('200 split copies, pairs within the chain on each: a pan frame (EC-48, EC-65)', async (within) => {
    const s = await settle(moleculesSpec({ chains: 200, copies: 'each', within: true }));
    const inputs: FrameInput[] = [];
    for (let k = 0; k < 8; k++) {
      const start = 4000 + k * 37;
      inputs.push({ size: Size, dpr: 2, viewport: { start, end: start + 300 } });
    }
    for (const i of inputs) await s.frame(i);
    const f = s.renderer.frame(inputs[0]);
    // three pairs on each of the 200 copies, the ends placed once a section
    const lay = f.layout.sections[0].sectionLays.find((l) => (l.lay as PairsLay).key.includes('/within/'))!
      .lay as PairsLay;
    expect(lay.prep!.units.get(f.layout.sections[0].section)!.n).toBe(600);
    expect(lay.n).toBeGreaterThan(0);
    const time = () => {
      const t: number[] = [];
      for (let k = 0; k < Frames; k++) {
        const t0 = performance.now();
        s.renderer.frame(inputs[(k + 1) % inputs.length]);
        t.push(performance.now() - t0);
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    time();
    const median = within(8 * Frames, time);
    console.info(`200 split copies, pairs on each, a pan: median ${median.toFixed(3)} ms`);
  });
});
