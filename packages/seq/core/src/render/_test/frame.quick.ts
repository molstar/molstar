/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The allocation targets (performance.md, Memory; _spec/README.md "Render", Draw lists): a warm draw (every pass of a
// laid-out view into its lists) and a warm hover (the renderer's frame with an equal key: the check and the emphasis
// redraw) allocate nothing; a pan or zoom frame (a layout and a draw) allocates at most a bound per scenario, a size
// that follows the screen and the view, not the data's size. `quickCase` (test-harness) reads the bytes allocated a
// step exactly; the inputs of a walk are made before it starts, so a step's bytes are the renderer's. By hand, never in
// CI: `pnpm -C packages/seq/core test:perf` (`--reporter=verbose` prints the bytes a step).

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import { projectOntoView } from '@molstar/seqx-core/view';
import { drawView } from '../index.js';
import { hoverInput, hovers, prepare, reusedHover, Scenarios, STEPS, walk } from './scenarios.js';

// a draw or a hover costs microseconds: the harness's own warm-up, to the code it settles in
const Options = {};
// a frame of the 3M bars costs a millisecond: its warm-up is capped by steps, so a case takes under a minute; a layout
// allocates kilobytes by design, and the control's 88 bytes are lost in that
const LayoutOptions = {
  warmSteps: 500,
  steps: 200,
  windows: 3,
  maxWarmSteps: 4000,
  control: false,
  why: 'a layout allocates kilobytes: they hide the control',
} as const;

/**
 * Bytes a pan or zoom frame may allocate, per scenario: the larger of its two walks' measurements (2026-11-18, M-series,
 * Node 26; the windows of a walk read up to a third apart) with half again as headroom. What a layout makes once (the
 * wants, the regions, the tiles' keys, the key's text) and a little an item in view (the 3M bars show 2,000 ranges
 * items at 100,000 positions: their keys and lane maps); never a buffer a position or a bin (44 MB a frame before).
 */
const LayoutBound: Record<string, number> = {
  // measured 25.7 KB
  'p53 tracks, letters': 40_000,
  'p53 tracks, whole protein': 40_000,
  // 40.5 KB
  'p53 on its gene, letters': 64_000,
  'p53 on its gene, whole protein': 64_000,
  // 27.9 KB
  '3M letters (120 positions)': 48_000,
  // 37.3 KB
  '3M cells (800 positions)': 64_000,
  // 765 KB (pan) at the default, a pixel's level; 730 KB at two pixels'
  '3M bars (100,000 positions)': 1_200_000,
  '3M bars (100,000 positions), level of two pixels': 1_200_000,
  // 21 KB (pan: level 64 throughout), 740 KB (zoom, a window 1.0 MB: below 160,000 positions the ranges are items)
  '3M bars (200,000 positions)': 1_500_000,
};

describe('warm draws and hovers allocate nothing; a layout a bounded size', () => {
  for (const s of Scenarios) {
    for (const kind of ['pan', 'zoom'] as const)
      it(`${s.name}: ${kind} (layout and draw)`, async ({ task }) => {
        const settled = await prepare(s);
        const inputs = walk(s, kind);
        const { r } = quickCase(
          task,
          `${s.name} ${kind}`,
          (i) => void settled.renderer.frame(inputs[i % STEPS]),
          LayoutOptions,
        );
        console.info(`${s.name} ${kind}: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
        expect(r.bytesPerStep, `${s.name} ${kind}: ${r.note}`).toBeLessThan(LayoutBound[s.name]);
      });

    it(`${s.name}: draw`, async ({ task }) => {
      const settled = await prepare(s);
      const frame = settled.renderer.frame(hoverInput(s));
      const { layout, lists } = frame;
      const { r } = quickCase(task, `${s.name} draw`, () => drawView(layout, lists), Options);
      console.info(`${s.name} draw: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, `${s.name} draw: ${r.note}`).toBeLessThanOrEqual(inject('perfThreshold'));
    });

    it(`${s.name}: hover`, async ({ task }) => {
      const settled = await prepare(s);
      const input = hoverInput(s);
      const { hits } = hovers(settled, input);
      const hover = reusedHover(hits);
      const at = { ...input, emphasis: hover.emphasis };
      const step = (i: number) => {
        hover.set(i);
        settled.renderer.frame(at);
      };
      const { r } = quickCase(task, `${s.name} hover`, step, Options);
      console.info(`${s.name} hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, `${s.name} hover: ${r.note}`).toBeLessThanOrEqual(inject('perfThreshold'));
    });

    it(`${s.name}: selection and cursor`, async ({ task }) => {
      const settled = await prepare(s);
      const input = hoverInput(s);
      const axis = settled.renderer.frame(input).layout.sections[0].axis!;
      const vis = axis.visible([]);
      const a = vis[0];
      const w = Math.max(1, Math.floor((vis[1] - vis[0]) / 8));
      // the viewer's records: a selection of a few ranges in view, the cursor walking across it
      const selections = [0, 1, 2, 3].map((k) => rangesLocation(axis.system, [a + k * w, a + (k + 1) * w]));
      const cursor = { system: axis.system, position: a };
      const emphasis: { selection?: Location; cursor?: typeof cursor } = { selection: selections[0], cursor };
      const at = { ...input, emphasis };
      const step = (i: number) => {
        emphasis.selection = selections[i & 3];
        cursor.position = a + (i % (4 * w));
        settled.renderer.frame(at);
      };
      const { r } = quickCase(task, `${s.name} selection and cursor`, step, Options);
      console.info(`${s.name} selection and cursor: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, `${s.name} selection and cursor: ${r.note}`).toBeLessThanOrEqual(inject('perfThreshold'));
    });

    it(`${s.name}: projected hover and selection, highlights`, async ({ task }) => {
      const settled = await prepare(s);
      const input = hoverInput(s);
      const axis = settled.renderer.frame(input).layout.sections[0].axis!;
      const vis = axis.visible([]);
      const a = vis[0];
      const w = Math.max(1, Math.floor((vis[1] - vis[0]) / 8));
      // the viewer's records (viewer.ts, `_updateEmphasis`): projections made when what they carry changes, cached;
      // a redraw only sets the fields to them
      const own = [0, 1, 2, 3].map((k) => rangesLocation(axis.system, [a + k * w, a + (k + 1) * w]));
      const projected = own.map((l) => projectOntoView(settled.model, settled.view, l));
      const hovered = own.map((l) => projectOntoView(settled.model, settled.view, l, { own: true, close: true }));
      const highlights = [{ location: hovered[2], color: 0x7aa3a6ff }];
      const p: { hover?: Location; selection?: Location } = {};
      const emphasis: { selection?: Location; projected: typeof p; highlights: typeof highlights } = {
        projected: p,
        highlights,
      };
      const at = { ...input, emphasis };
      const step = (i: number) => {
        emphasis.selection = own[i & 3];
        p.selection = projected[i & 3];
        p.hover = hovered[(i + 1) & 3];
        settled.renderer.frame(at);
      };
      const { r } = quickCase(task, `${s.name} projected`, step, Options);
      console.info(`${s.name} projected: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, `${s.name} projected: ${r.note}`).toBeLessThanOrEqual(inject('perfThreshold'));
    });
  }
});
