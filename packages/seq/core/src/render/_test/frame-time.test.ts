/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

// The frame-time targets (performance.md, Time), measured through `Renderer.frame` at 1440 CSS px and dpr 2, warm, over
// 200 steps: a pan or zoom frame under 8 ms (median), never over 16 ms (its maximum), a hover under 2 ms (the
// hit-test and the emphasis redraw), the SVG export of p53 as tracks under 200 ms. Each case prints its median, 95th
// percentile and maximum. Run: `pnpm test:heavy packages/seq/core/src/render`.

import { beforeAll, describe, expect, it } from 'vitest';
import { hitView, renderSvg } from '../index.js';
import { p53TracksSpec } from './fixtures.js';
import { hoverInput, hovers, prepare, Scenarios, stats, STEPS, time, walk } from './scenarios.js';
import type { Settled } from './helpers.js';

const PAN_MS = 8;
const FRAME_MAX_MS = 16;
const HOVER_MS = 2;

const report = (what: string, t: number[]) => {
  const s = stats(t);
  console.info(`${what}: median ${s.median.toFixed(3)} ms, p95 ${s.p95.toFixed(3)} ms, max ${s.max.toFixed(3)} ms`);
  return s;
};

for (const s of Scenarios) {
  describe(s.name, () => {
    let settled: Settled;
    beforeAll(async () => {
      settled = await prepare(s);
    }, 60_000);

    for (const kind of ['pan', 'zoom'] as const) {
      it(`${kind} frame: median under ${PAN_MS} ms, never over ${FRAME_MAX_MS} ms`, { tags: ['heavy'] }, () => {
        const inputs = walk(s, kind);
        const t = time((i) => void settled.renderer.frame(inputs[i % STEPS]));
        const r = report(`${s.name} ${kind}`, t);
        expect(r.median).toBeLessThan(PAN_MS);
        expect(r.p95).toBeLessThan(FRAME_MAX_MS);
        expect(r.max).toBeLessThan(FRAME_MAX_MS);
      });
    }

    it(`hover: hit-test and emphasis redraw under ${HOVER_MS} ms (median)`, { tags: ['heavy'] }, () => {
      const input = hoverInput(s);
      const { points, emphasis } = hovers(settled, input);
      expect(points.length).toBeGreaterThan(0);
      const layout = settled.renderer.frame(input).layout;
      const withEmphasis = emphasis.map((e) => ({ ...input, emphasis: e }));
      const frame = report(
        `${s.name} hover frame`,
        time((i) => void settled.renderer.frame(withEmphasis[i % STEPS])),
      );
      const hit = report(
        `${s.name} hit-test`,
        time((i) => void hitView(layout, points[i % points.length].x, points[i % points.length].y)),
      );
      const both = report(
        `${s.name} hover (hit-test + frame)`,
        time((i) => {
          const p = points[i % points.length];
          hitView(layout, p.x, p.y);
          settled.renderer.frame(withEmphasis[i % STEPS]);
        }),
      );
      // the hover pass draws something for most of the walk (a position hit over a pending or bar layer draws nothing)
      let drawn = 0;
      for (const e of withEmphasis) {
        const f = settled.renderer.frame(e);
        expect(f.rebuilt).toBe(false);
        if (f.lists.passes.hover.size > 0) drawn++;
      }
      console.info(`${s.name} hover: ${drawn} of ${STEPS} steps draw a hover mark`);
      // measured: over bars (below a device pixel a position) no hover mark is drawn at all, so those steps time the pass's reset alone
      if (!s.name.includes('bars')) expect(drawn).toBeGreaterThan(0);
      expect(both.median).toBeLessThan(HOVER_MS);
      expect(frame.median + hit.median).toBeLessThan(HOVER_MS + 0.5);
    });
  });
}

describe('SVG export', () => {
  it('p53 as tracks, 1440 px at dpr 2: under 200 ms (median over 50 runs, warm)', { tags: ['heavy'] }, async () => {
    const run = () => renderSvg(p53TracksSpec(), { size: { width: 1440, height: 160 }, dpr: 2 });
    for (let i = 0; i < 20; i++) await run();
    const t: number[] = [];
    for (let i = 0; i < 50; i++) {
      const t0 = performance.now();
      const r = await run();
      t.push(performance.now() - t0);
      expect(r.settled).toBe(true);
    }
    const r = report('renderSvg p53 tracks', t);
    expect(r.median).toBeLessThan(200);
  });

  it('p53 as tracks: the first (cold) export, reported', { tags: ['heavy'] }, async () => {
    const t0 = performance.now();
    await renderSvg(p53TracksSpec(), { size: { width: 1440, height: 160 }, dpr: 2 });
    console.info(`renderSvg p53 tracks, first call in this process: ${(performance.now() - t0).toFixed(1)} ms`);
  });
});
