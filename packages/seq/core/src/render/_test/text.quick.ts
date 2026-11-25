/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grouped axis's placement (the text layout: groups of ten with a gap after each, _spec/README.md "Render", Lanes and
// Hits): `at` reads the page point and the strip-frame x through arrays, so placing a point on a letter or in a group's
// gap allocates nothing. The whole text `hitView` is measured and reported, not bounded: its other hit paths still
// allocate (the package's TODO). By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { hitView, PlacePoint } from '../index.js';
import { p53TextSpec } from './fixtures.js';
import { settle } from './helpers.js';

describe('a grouped axis places a point with nothing allocated', () => {
  it('the text layout: points on letters and in the groups’ gaps', async ({ task }) => {
    const s = await settle(p53TextSpec());
    const f = await s.frame({ size: { width: 700, height: 900 }, dpr: 1.5, settled: true });
    const sec = f.layout.sections[0];
    const p = sec.placement;
    const a = p.strips.axes[0];
    expect(a.gsize).toBeGreaterThan(0);
    // fractional page points over the first lines' letters band: on letters and in the gaps between groups
    const points: number[] = [];
    const pt = new PlacePoint();
    let gaps = 0;
    for (let i = 0; i < Math.min(4, p.strips.n); i++) {
      const ax = p.strips.axes[i];
      // the strip's first band (its letters), in page coordinates
      const y = p.strips.dy[i] + p.strips.bandY[i * p.strips.bands] + 1.37;
      for (let k = 0; k < 40; k++) {
        const x = p.strips.dx[i] + ax.x0 + ((k + 0.31) * ax.width) / 40;
        points.push(x, y);
        expect(p.at([x, y], pt)).toBe(true);
        if (pt.kind === 'gap') gaps++;
      }
    }
    expect(gaps).toBeGreaterThan(0);
    const flat = Float64Array.from(points);
    const m = flat.length / 2;
    const PA = new Float64Array(2);
    const probe = quickCase(
      task,
      'text placement at grouped',
      (i) => {
        PA[0] = flat[2 * (i % m)];
        PA[1] = flat[2 * (i % m) + 1];
        p.at(PA, pt);
      },
      {},
    );
    console.info(`text at (grouped): ${probe.r.bytesPerStep.toFixed(0)} bytes (${probe.r.note})`);
    expect(probe.r.bytesPerStep, probe.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const hit = quickCase(
      task,
      'text hitView grouped',
      (i) => void hitView(f.layout, flat[2 * (i % m)], flat[2 * (i % m) + 1]),
      {},
    );
    console.info(`text hitView (grouped): ${hit.r.bytesPerStep.toFixed(0)} bytes a hit (${hit.r.note})`);
  });
});
