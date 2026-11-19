/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows' allocation targets (the design's "Performance and size"): a warm draw at letters and at the fit (image
// tiles) and a warm hover over the cells allocate nothing; the selection's rows and the cursor on a cell nothing; a
// layout a bounded size a pan frame, the same at 100,000 and 1,000,000 lazy rows. By hand, never in CI:
// `pnpm -C packages/seq/core test:perf`.

import { rangesLocation } from '@molstar/seqx-core/query';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, type Hit, hitView, type RowsLay } from '../index.js';
import { settle } from './helpers.js';
import { Columns, rowsSpec, Section, Synthetic, Track } from './rows-fixtures.js';
import { reusedHover } from './scenarios.js';

const Size = { width: 1440, height: 800 };
const at = (bandTop: number, rowScale: number, start = 0, end = Columns): FrameInput => ({
  size: Size,
  dpr: 2,
  settled: true,
  viewports: { [Section]: { start, end, bands: { [Track]: { bandTop, rowScale } } } },
});

const LayoutOptions = {
  warmSteps: 300,
  steps: 120,
  windows: 3,
  maxWarmSteps: 2000,
  control: false,
  why: 'a layout allocates kilobytes: they hide the control',
} as const;

/** Bytes a pan frame of the rows may allocate: what a layout makes once and its labels and tiles in view. */
const LayoutBound = 120_000;

function rowsLay(f: {
  layout: { sections: readonly { lays: readonly (readonly { lay: unknown; type: { name: string } }[])[] }[] };
}) {
  for (const s of f.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
}

describe('the rows draw and hover with nothing allocated once warm; a layout a bounded size', () => {
  for (const [name, input] of [
    ['letters', at(20 * 4000, 20, 100, 160)],
    ['the fit (image tiles)', at(0, 0.06)],
  ] as const)
    it(`10,000 × 300 at ${name}: draw`, async ({ task }) => {
      const s = await settle(rowsSpec(10_000, true), [Synthetic]);
      const f = await s.frame(input);
      const { r } = quickCase(task, `rows draw at ${name}`, () => drawView(f.layout, f.lists), {});
      console.info(`rows draw at ${name}: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    });

  it('10,000 × 300 at letters: hover over the cells, the selection’s rows and the cursor on a cell', async ({
    task,
  }) => {
    const s = await settle(rowsSpec(10_000, true), [Synthetic]);
    const input = at(20 * 4000, 20, 100, 160);
    const f = await s.frame(input);
    const lay = rowsLay(f);
    const hits: Hit[] = [];
    for (let i = lay.first; i < lay.last; i += 3)
      for (let j = 100; j < 160; j += 7) {
        if (!lay.axis.cell(j)) continue;
        const h = hitView(f.layout, (lay.axis.sx0 + lay.axis.sx1) / 2, lay.rowY(i) + 5);
        if (h) hits.push(h);
      }
    const hover = reusedHover(hits);
    const hoverAt = { ...input, emphasis: hover.emphasis };
    const h = quickCase(
      task,
      'rows hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(hoverAt);
      },
      {},
    );
    console.info(`rows hover: ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const rows = lay.rowsSystem;
    const selections = [0, 1, 2, 3].map((k) => rangesLocation(rows, [lay.first + k, lay.first + k + 3]));
    const cursor = { system: lay.system, position: 100, ySystem: rows, yPosition: lay.first };
    const emphasis: { selection?: (typeof selections)[number]; cursor?: typeof cursor } = { cursor };
    const selAt = { ...input, emphasis };
    const c = quickCase(
      task,
      'rows selection and cursor',
      (i) => {
        emphasis.selection = selections[i & 3];
        cursor.position = 100 + (i % 50);
        s.renderer.frame(selAt);
      },
      {},
    );
    console.info(`rows selection and cursor: ${c.r.bytesPerStep.toFixed(0)} bytes a frame (${c.r.note})`);
    expect(c.r.bytesPerStep, c.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  for (const n of [100_000, 1_000_000])
    it(`${n.toLocaleString('en-US')} lazy rows: a pan frame at letters, bounded`, async ({ task }) => {
      const s = await settle(rowsSpec(n, false), [Synthetic]);
      const inputs = [0, 1, 2, 3].map((k) => at(20 * (n / 2) + 7 * k, 20, 100, 160));
      for (const i of inputs) await s.frame(i);
      const { r } = quickCase(task, `rows pan ${n}`, (i) => void s.renderer.frame(inputs[i & 3]), LayoutOptions);
      console.info(`rows pan, ${n} lazy rows: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThan(LayoutBound);
    });
});
