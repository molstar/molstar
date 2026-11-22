/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover's allocation targets (_spec/README.md "Render", A liftover): a warm draw of the figure (the old genes
// lifted onto the new assembly through a map read by region, the bridge's ribbons) and a warm hover over the lifted
// genes allocate nothing. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { quickCase } from '@molstar/test-harness/quick';
import { beforeAll, describe, expect, inject, it } from 'vitest';
import { loadBridges } from '../bridge-load.js';
import { drawView, type FrameInput, type Hit, hitView } from '../index.js';
import { settle } from './helpers.js';
import { liftSpec, NewSection } from './lift-fixtures.js';
import { reusedHover } from './scenarios.js';

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
});

const input: FrameInput = { size: { width: 1440, height: 600 }, dpr: 2, settled: true };

describe('a liftover draws and hovers with nothing allocated once warm', () => {
  it('the lifted genes: a warm frame, a hover over them', async ({ task }) => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const f = await s.frame(input);
    const d = quickCase(task, 'liftover draw', () => drawView(f.layout, f.lists), {});
    console.info(`liftover draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    // the lifted genes' hits on the new section
    const sec = f.layout.sections.find((x) => x.section.address === NewSection)!;
    const hits: Hit[] = [];
    const keys = new Set<string>();
    for (let y = sec.y; y < sec.y + sec.height; y += 2)
      for (let x = 0; x < 1440; x += 3) {
        const h = hitView(f.layout, x, y);
        if (h?.element && h.layer.endsWith('/new/genes/genes') && !keys.has(`${h.system}|${h.element.key}`)) {
          keys.add(`${h.system}|${h.element.key}`);
          hits.push({ ...h });
        }
      }
    expect(hits.length).toBeGreaterThan(4);
    const hover = reusedHover(hits);
    const hoverAt = { ...input, emphasis: hover.emphasis };
    const h = quickCase(
      task,
      'liftover hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(hoverAt);
      },
      {},
    );
    console.info(`liftover hover: ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
