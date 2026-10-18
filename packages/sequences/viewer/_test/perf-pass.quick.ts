/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The performance pass's quick cases (seqspec:viewer#performance; perf-pass-cases.ts), each with its
// positive control: quick-only, as a new hot path's first guard is (packages/sequences/AGENTS.md, Tests). The input
// storms (the unobserved context menu, the streams with no subscriber and after one unsubscribed) are input.quick.ts's;
// the strip lanes wrapped-lanes.quick.ts's.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import {
  centeredText,
  chips,
  FRAME,
  TEXT_TRACKS,
  textTracks,
  unobservedLinear,
  unobservedText,
} from './perf-pass-cases.js';

describe('unobserved positions (unobservedGray), quick', () => {
  it("letters' cells zoomed in on a linear section over many observed ranges, panned", async ({ task }) => {
    const d = await unobservedLinear();
    quickCase(task, 'perf-pass: unobserved linear cells', d.step, FRAME);
    expect(d.be.rectItems.n).toBeGreaterThan(50);
  });

  it.for([true, false])(
    'wrapped text over underline and strip lanes, scrolled (wrappedColors %s)',
    async (tinted, { task }) => {
      const d = await unobservedText(tinted);
      quickCase(task, `perf-pass: unobserved text${tinted ? ', tinted' : ''}`, d.step, FRAME);
      expect(d.be.rectItems.n).toBeGreaterThan(50);
    },
  );
});

describe('a centered reading column (wrappedMaxWidth), quick', () => {
  it('a select drag over its lines, the pending range drawn each move', async ({ task }) => {
    const d = await centeredText('select');
    quickCase(task, 'perf-pass: centered text select', d.step);
    const s = d.state();
    expect([s.owner, s.active]).toEqual(['select', true]);
    expect(s.ends).toBeGreaterThan(1_000);
  });

  it('hits on its lines', async ({ task }) => {
    const d = await centeredText('hit');
    quickCase(task, 'perf-pass: centered text hits', d.step);
    expect(d.state().hits).toBeGreaterThan(1_000);
  });
});

describe('chips (short units in the text), quick', () => {
  it('the pointer over chips, an overlay frame each', async ({ task }) => {
    const d = await chips('hover');
    quickCase(task, 'perf-pass: chips hover', d.step);
    expect(d.hovered()).toBeGreaterThan(1_000);
  });

  it('rows of chips scrolled', async ({ task }) => {
    const d = await chips('scroll');
    quickCase(task, 'perf-pass: chips scroll', d.step, FRAME);
    expect(d.be.rectItems.n).toBeGreaterThan(0);
  });
});

describe('tracks drawn in full in the text (wrappedTracks: full), quick', () => {
  it.for(TEXT_TRACKS)('%s, scrolled', async (kind, { task }) => {
    const d = await textTracks(kind);
    quickCase(task, `perf-pass: text tracks ${kind}`, d.step, FRAME);
    expect(d.be.opCount).toBeGreaterThan(0);
  });
});
