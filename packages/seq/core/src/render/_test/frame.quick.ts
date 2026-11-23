/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The allocation targets (_spec/README.md "Render", Draw lists): a warm draw (every pass of a
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
} as const;

/**
 * Bytes a pan or zoom frame may allocate, per scenario: the larger of its two walks' measurements (2026-11-18, M-series,
 * Node 26; the windows of a walk read up to a third apart) with half again as headroom. What a layout makes once (the
 * wants, the regions, the tiles' keys, the key's text) and a little an item in view (the 3M bars show 2,000 ranges
 * items at 100,000 positions: their keys and lane maps); never a buffer a position or a bin (44 MB a frame before).
 */
const LayoutBound: Record<string, number> = {
  // 44.5 KB (zoom; pan 42.3 KB, 2026-11-19): three tracks, the ruler's ticks and labels, the run and its note
  'holey rulers': 64_000,
  // measured 25.7 KB
  'p53 tracks, letters': 40_000,
  // 43.6 KB (pan; zoom 36.7 KB, 2026-11-20): one block for 200 copies, three chips, the marks, the links' arcs
  'complex, 200 copies': 70_000,
  'p53 tracks, whole protein': 40_000,
  // 65.0 KB (seven tracks: the ruler's labels, secondary structure's and the variants' marks and sorts, the lines' paths)
  'p53 small layers, whole protein': 100_000,
  // 58.5 KB
  'p53 small layers, letters': 90_000,
  // 40.5 KB
  'p53 on its gene, letters': 64_000,
  'p53 on its gene, whole protein': 64_000,
  // 64.8 KB (zoom; pan 58.4 KB): the lines in view laid out again, their transforms and marks
  'p53 text': 100_000,
  // 161 KB (2026-11-19): the codon path's parts made a letter (`projectRanges` clears its buffer, TODO)
  'p53 on its gene, coding': 240_000,
  // 272 KB (313 KB run first): the breaks and pieces in view, 1,000 of each
  'reversed, 150,000 pieces': 470_000,
  // 121.6 KB (pan; zoom 110.4 KB, 2026-11-21): the transcripts in view, their lanes by item, boxes and hairlines
  'gene models (60 kb)': 182_000,
  // scrolled band: the lanes in view only
  'gene models, scrolled (60 kb)': 190_000,
  'gene models, scrolled far (60 kb)': 190_000,
  // the lanes rail dragged: as the scrolled band, the rail's window elsewhere
  'gene models, rail dragged (60 kb)': 190_000,
  // 27.9 KB
  '3M letters (120 positions)': 48_000,
  // 37.3 KB
  '3M cells (800 positions)': 64_000,
  // 114 KB (pan; zoom 94 KB, 2026-11-19): 160 items in view, their marks, arrowheads and lane maps; their labels too
  // long for the marks, skipped unmeasured
  '3M items with arrowheads (8,000 positions)': 180_000,
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

    it(`${s.name}: hover with its projection (the view's mark where a layer drew none)`, async ({ task }) => {
      // the hover as the viewer draws it (viewer-emphasis.ts): its hit, and its projection with its own section, so a
      // hit its layer draws no mark for (a position, a density bin) takes the view's (Emphasis, The hover contract)
      const settled = await prepare(s);
      const input = hoverInput(s);
      const layout = settled.renderer.frame(input).layout;
      const { hits } = hovers(settled, input);
      const hover = reusedHover(hits);
      const shown = hits.map((h) =>
        projectOntoView(
          settled.model,
          settled.view,
          h.extent ?? rangesLocation(h.system, h.position === undefined ? [] : [h.position, h.position + 1]),
          { own: true, close: true },
        ),
      );
      const sections = hits.map((h) => layout.sections.find((x) => x.id === h.section)?.section.address);
      const p: { hover?: Location; hoverSection?: string } = {};
      const at = { ...input, emphasis: { ...hover.emphasis, projected: p } };
      const step = (i: number) => {
        // the hit record is the emphasis's (`reusedHover`: one object, its fields set)
        hover.set(i);
        p.hover = shown[i % shown.length];
        p.hoverSection = sections[i % sections.length];
        settled.renderer.frame(at);
      };
      const { r } = quickCase(task, `${s.name} own hover`, step, Options);
      console.info(`${s.name} own hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, `${s.name} own hover: ${r.note}`).toBeLessThanOrEqual(inject('perfThreshold'));
    });

    it(`${s.name}: overlapping highlight channels`, async ({ task }) => {
      const settled = await prepare(s);
      const input = hoverInput(s);
      const axis = settled.renderer.frame(input).layout.sections[0].axis!;
      const vis = axis.visible([]);
      const a = vis[0];
      const w = Math.max(2, Math.floor((vis[1] - vis[0]) / 8));
      // four channels overlapping in turn, in each style: stripes, merged runs and the styles' bars and outlines
      const over = (k: number) => rangesLocation(axis.system, [a + k * w, a + (k + 3) * w]);
      const made = (k: number) => [
        { location: over(k), color: 0x7aa3a6ff },
        { location: over(k + 1), color: 0xd86c17ff, style: 'band' as const },
        { location: over(k + 2), color: 0x2a9d8fff, style: 'outline' as const },
        { location: over(k + 3), color: 0xe9c46aff },
      ];
      // two sets the step alternates between (made before the walk: a changed list is a redrawn emphasis)
      const sets = [made(0), made(1)];
      const emphasis = { highlights: sets[0] };
      const at = { ...input, emphasis };
      const step = (i: number) => {
        emphasis.highlights = sets[i & 1];
        settled.renderer.frame(at);
      };
      const { r } = quickCase(task, `${s.name} overlapping highlights`, step, Options);
      console.info(`${s.name} overlapping highlights: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, `${s.name} overlapping highlights: ${r.note}`).toBeLessThanOrEqual(
        inject('perfThreshold'),
      );
    });
  }
});
