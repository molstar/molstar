/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of input.heap.ts (`test:perf`): IN2's pointer-rate paths (input-cases.ts, the guard's cases), each
// with its control: moves within an element at every spot, element changes, items, a wrapped section and matrix bands,
// the label churn (a step is one hover change: eight moves, its label's strings by design), drags, a press, drag and
// release (a step one cycle), the wheel; and the public streams with no subscriber: element changes and a storm of
// context-menu requests unobserved, and element changes flat again once a hover subscriber unsubscribed.

import { describe, expect, it } from 'vitest';
import { case3, case10 } from '../../core/_test/format/cases.js';
import { measureBytesPerStep, quickCase } from './heap.js';
import {
  CHURNS,
  case3Spots,
  churnCase,
  cycleDoc,
  cycling,
  DRAGS,
  dragDoc,
  elementChanges,
  familySpots,
  featureItems,
  hoverStream,
  matrixBands,
  menuStorm,
  railSpots,
  rowsSweep,
  type Spot,
  scrollGroupDrag,
  startDrag,
  wheelCases,
  within,
  wrappedMoves,
} from './input-cases.js';

const DOCS = { case3, case10 };
/** Every spot's moves within its element, a case each. */
const each = (task: Parameters<typeof quickCase>[0], spots: Spot[]) => {
  for (const s of spots) quickCase(task, `input: within, ${s.name}`, within(s));
};

describe('moves and element changes, quick', () => {
  it('moves within an element over every region of a family with a tree, labels on', async ({ task }) => {
    const { spots } = await familySpots();
    expect(spots.length).toBeGreaterThan(5);
    each(task, spots);
  });

  it("over case 3 and its navigator's body, edge and outside", async ({ task }) => {
    each(task, await case3Spots());
  });

  it('over the rulerRail', async ({ task }) => {
    const spots = await railSpots();
    expect(spots.some((s) => s.name.includes('region 7'))).toBe(true);
    each(task, spots);
  });

  it('an element change every 8 moves, labels off', async ({ task }) => {
    const { step, changes } = await elementChanges();
    // A higher cap: even under the quick flags these changes allocate 8 to 20 bytes a step until about the 18,000th
    // (2026-10-15, Node 24, M4); the warm-up runs on until it reads quiet.
    quickCase(task, 'input: element changes', step, { warmSteps: 2_000, maxWarmSteps: 60_000 });
    expect(changes()).toBeGreaterThan(1_000);
  });

  it("a points track's items and split rows, a wrapped section, matrix bands", async ({ task }) => {
    const items = await featureItems();
    quickCase(task, 'input: items and split rows', items.step);
    expect(items.seen()).toBeGreaterThan(0);
    const wrapped = await wrappedMoves();
    expect(wrapped.lines).toBeGreaterThan(5);
    quickCase(task, 'input: wrapped', wrapped.step);
    const bands = await matrixBands();
    quickCase(task, 'input: matrix bands', bands.step);
    expect(bands.states()).toBeGreaterThan(1_000);
  });

  it.for(CHURNS)('the label churn over a $name: a label per hover change', async (c, { task }) => {
    const churn = await churnCase(DOCS[c.doc](), c.section, c.track, c.row, 20);
    // A step is one change: eight moves over one cell, the next step's over the next.
    const step = (k: number) => {
      for (let j = 0; j < 8; j++) churn.moved(8 * k + j);
    };
    churn.reset();
    // Allocates by design (a label's strings), so never quiet: a short cap (alloc.ts).
    quickCase(task, `input: label churn, ${c.name}`, step, { warmSteps: 500, steps: 250, maxWarmSteps: 4_000 });
    const { computed, changes } = churn.counts();
    expect(changes).toBeGreaterThan(0);
    expect(computed).toBe(changes);
  });
});

describe('drags, quick', () => {
  it("dragPan, a select drag, the navigator's body and edge, the scrollbar's thumb", async ({ task }) => {
    const t = await dragDoc();
    for (const d of DRAGS) {
      const { step, release } = startDrag(t, d);
      quickCase(task, `input: drag, ${d.name}`, step);
      expect([d.name, release()]).toEqual([d.name, true]);
    }
  });

  it('dragPan in a scroll group of three', async ({ task }) => {
    const { step, release, scrollX } = await scrollGroupDrag();
    quickCase(task, 'input: drag, scroll group', step);
    release();
    const [one, two, three] = scrollX();
    expect([two, three]).toEqual([one, one]);
  });

  it.for([
    { name: 'dragPan', alt: false },
    { name: 'select', alt: true },
  ])('a full press, drag and release: $name, a step a cycle', async ({ name, alt }, { task }) => {
    const t = await cycleDoc();
    const cycle = cycling(t);
    // A cycle allocates by design (the drag's records), so never quiet: a short cap (alloc.ts).
    quickCase(task, `input: press-drag-release, ${name}`, (k) => cycle(alt, k), {
      warmSteps: 1_000,
      steps: 250,
      maxWarmSteps: 4_000,
    });
    expect(t.ctx.ui.drag.value).toBeNull();
  });
});

describe('the wheel, quick', () => {
  it("Ctrl zoom, a horizontal pan, Shift with only deltaX and only deltaY, a box section's rows", async ({ task }) => {
    const { cases, version } = await wheelCases();
    const v0 = version();
    for (const [name, step] of cases) quickCase(task, `input: wheel, ${name}`, step);
    expect(version()).toBeGreaterThan(v0 + 10_000);
  });

  it('a rows sweep under a resting pointer, the hover hit again each flush', async ({ task }) => {
    const { t, step, changes } = await rowsSweep();
    quickCase(task, 'input: wheel, rows sweep', step);
    expect(changes()).toBeGreaterThan(1_000);
    expect(t.ctx.hover.live.current.owner).toBe('pointer');
  });
});

describe('the public streams, unobserved, quick', () => {
  it("a storm of the browser's context-menu requests, no menu subscriber", async ({ task }) => {
    const { step, asked, observed } = await menuStorm();
    expect(observed()).toBe(false);
    quickCase(task, 'input: streams, context menu unobserved', step);
    expect(asked()).toBeGreaterThan(1_000);
  });

  it('element changes unobserved; with a hover subscriber a payload each; once it unsubscribed, flat again', async ({
    task,
  }) => {
    const { step, subscribe, observed } = await hoverStream();
    quickCase(task, 'input: streams, hover unobserved', step);
    const sub = subscribe();
    // Not vacuous: with the subscriber the changes allocate (a payload each, by design).
    const observedRun = measureBytesPerStep(step);
    expect(observedRun.bytesPerStep, observedRun.note).toBeGreaterThan(4);
    sub.unsubscribe();
    expect(observed()).toBe(false);
    quickCase(task, 'input: streams, hover unsubscribed', step);
  });
});
