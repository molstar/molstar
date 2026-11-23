/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text layout's targets: a reflow frame (line breaking and line heights, then the
// lines in view laid out and drawn) under 8 ms for the p53 text figure and a 35,000-residue sequence, flat in the
// sequence's length (35,000 against 3,500 at most 1.2×); a hover (the placement's `at`, the layer's hit, the emphasis
// redraw) under 2 ms. Each bound is the heavy lane's (`itWithin`, best of three); the quick lane runs them unbounded.

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { type FrameInput, hitView } from '../index.js';
import { settle } from './helpers.js';
import { p53TextSpec, seeded } from './fixtures.js';

/** A protein of `n` residues as text, with a domain a 500 residues under its letters. */
function longTextSpec(n: number): Node {
  const r = seeded(5);
  const alphabet = 'ACDEFGHIKLMNPQRSTVWY';
  let letters = '';
  for (let i = 0; i < n; i++) letters += alphabet[Math.floor(r() * alphabet.length)];
  const starts: number[] = [];
  const ends: number[] = [];
  const ids: string[] = [];
  for (let s = 0; s + 200 < n; s += 500) {
    starts.push(s);
    ends.push(s + 200);
    ids.push(`d${s}`);
  }
  const b = createSeqSpec({ title: 'long text' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters });
  const d = b.data.ranges('d', {
    target: p,
    rows: { start: Int32Array.from(starts), end: Int32Array.from(ends), id: categoryColumn(ids) },
  });
  b.views.view('v', { layout: 'text' }, (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('domains').layer('ranges', { type: 'ranges', data: d });
  });
  return b.build();
}

const Widths = [1440, 1200, 1000, 1320];
const Frames = 40;

/** The median time of a reflow frame: the width changes every frame (each lays out and draws again). */
async function reflow(spec: Node, start: number): Promise<{ median: number; time: () => number }> {
  const s = await settle(spec);
  const inputs: FrameInput[] = Widths.map((width) => ({
    size: { width, height: 800 },
    dpr: 2,
    viewport: { start, end: start },
  }));
  for (const i of inputs) await s.frame(i);
  const time = () => {
    const t: number[] = [];
    for (let k = 0; k < Frames; k++) {
      const t0 = performance.now();
      const f = s.renderer.frame(inputs[k % inputs.length]);
      t.push(performance.now() - t0);
      if (!f.rebuilt) throw new Error('a reflow frame lays out again');
    }
    t.sort((a, b) => a - b);
    return t[t.length >> 1];
  };
  time();
  return { median: time(), time };
}

describe('the text layout’s time', () => {
  itWithin('a reflow frame of the p53 text figure', async (within) => {
    const r = await reflow(p53TextSpec(), 100);
    // the median of a walk of reflow frames under 8 ms (the walk's total under 8 ms a frame)
    const median = within(8 * Frames, () => r.time());
    console.info(`p53 text reflow: median ${median.toFixed(3)} ms`);
    expect(r.median).toBeGreaterThan(0);
  });

  itWithin('a reflow frame of a 35,000-residue text', async (within) => {
    const r = await reflow(longTextSpec(35_000), 20_000);
    const median = within(8 * Frames, () => r.time());
    console.info(`35,000-residue text reflow: median ${median.toFixed(3)} ms`);
    expect(r.median).toBeGreaterThan(0);
  });

  it('a reflow is flat in the sequence’s length: 35,000 residues against 3,500 at most 1.2×', {
    tags: ['heavy'],
  }, async () => {
    const small = await reflow(longTextSpec(3_500), 2_000);
    const large = await reflow(longTextSpec(35_000), 20_000);
    let best = Infinity;
    for (let i = 0; i < 3; i++) best = Math.min(best, large.time() / small.time());
    console.info(
      `text reflow: 3,500 ${small.median.toFixed(3)} ms, 35,000 ${large.median.toFixed(3)} ms, ratio ${best.toFixed(2)}`,
    );
    expect(best).toBeLessThan(1.2);
  });

  itWithin('a hover over the text: the placement, the hit, the emphasis redraw', async (within) => {
    const s = await settle(p53TextSpec());
    const input: FrameInput = { size: { width: 1440, height: 800 }, dpr: 2 };
    const f = await s.frame(input);
    const layout = f.layout;
    const p = layout.sections[0].placement;
    const points: { x: number; y: number }[] = [];
    for (let i = 0; i < p.strips.n; i++) {
      const a = p.strips.axes[i];
      for (let u = a.u0; u < a.u1; u += 7)
        points.push({ x: a.xOfUnit(u + 0.5), y: p.strips.top[i] + p.strips.bandY[i * 2] + 3 });
    }
    const emphasis = points.map((q) => ({ ...input, emphasis: { hover: hitView(layout, q.x, q.y) } }));
    const hover = () => {
      const t: number[] = [];
      for (let k = 0; k < 200; k++) {
        const q = points[k % points.length];
        const t0 = performance.now();
        hitView(layout, q.x, q.y);
        s.renderer.frame(emphasis[k % emphasis.length]);
        t.push(performance.now() - t0);
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    hover();
    const median = within(2 * 200, hover);
    console.info(`text hover: median ${median.toFixed(3)} ms`);
    expect(hitView(layout, points[3].x, points[3].y)?.position).toBe(21);
  });
});
