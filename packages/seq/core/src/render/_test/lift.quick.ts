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
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { beforeAll, describe, expect, inject, it } from 'vitest';
import { loadBridges } from '../bridge-load.js';
import { drawView, type FrameInput, type Hit, hitView } from '../index.js';
import { settle } from './helpers.js';
import { liftSpec, NewSection, OldSection } from './lift-fixtures.js';
import { reusedHover } from './scenarios.js';
import { regionProjector } from '../lift-project.js';

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

  it('a hover carried across the map and its card’s places, both ways', async ({ task }) => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    await s.frame(input);
    const p = regionProjector(s.store);
    const old = '/data/old/(sequence :name chr1)';
    const neu = '/data/new/(sequence :name chr2)';
    const other = '/data/new/(sequence :name chr1)';
    const gene = [72_000, 78_000];
    const image = [62_000, 68_000];
    // the caller's buffer, cut back to a kept prefix each step: V8 drops the backing store of an array set to length 0
    // (and trims one cut far below its capacity), so emptying it would measure its regrowth (152 B), not the projector
    const Kept = 64;
    const out: number[] = new Array(Kept).fill(0);
    const step = (i: number) => {
      out.length = Kept;
      if (i & 1) {
        p.project(s.model, neu, image, old, out);
        p.project(s.model, neu, image, other, out);
        p.places!(s.model, s.view, neu, image[0], image[1]);
      } else {
        p.project(s.model, old, gene, neu, out);
        p.places!(s.model, s.view, old, gene[0], gene[1]);
      }
    };
    step(0);
    step(1);
    expect(out.slice(Kept)).toEqual([72_000, 78_000, 62_000, 68_000]);
    const r = quickCase(task, 'liftover projection', step, {});
    console.info(`liftover projection: ${r.r.bytesPerStep.toFixed(0)} bytes a hover (${r.r.note})`);
    expect(r.r.bytesPerStep, r.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});

/**
 * A pan of a section lays out that section and the bridge again (counted in `lift.test.ts`): about a megabyte a frame
 * here, the bridge's ribbons and the section's strips (2026-11-22: 1.12 MB with a ruler alone on the new section). The
 * bound is the same genes read directly: the old section panned alike (1.50 MB). A lifted layer keeps its rows by the
 * runs' tiles and the window's level, so its pan costs no more than the plain layer's (1.30 MB); rebuilding its runs,
 * maps and far reads each frame would cost more.
 */
const PanSteps = { warmSteps: 200, maxWarmSteps: 2_000, steps: 200 };

async function panner(section: string) {
  const s = await settle(liftSpec(), [SeqRuntime]);
  const at = (x: number) => ({ ...input, viewports: { [section]: { start: 1_000 + x, end: 120_000 + x } } });
  const frames = [at(0), at(37)];
  for (const f of frames) await s.frame(f);
  return (i: number) => void s.renderer.frame(frames[i & 1]);
}

describe('a liftover laid out again', () => {
  it('a pan of the new section costs no more than the genes’ own section panned alike', async ({ task }) => {
    const lifted = await panner(NewSection);
    const plain = await panner(OldSection);
    const reads = [lifted, plain, lifted, plain].map((s) => measureBytesPerStep(s, PanSteps).bytesPerStep);
    const own = Math.min(reads[1], reads[3]);
    // a frame of a megabyte: `leak`'s 88 bytes are under its noise, so no control
    const { r } = quickCase(task, 'liftover pan', lifted, { ...PanSteps, control: false, why: 'a megabyte a frame' });
    console.info(
      `liftover pan: ${r.bytesPerStep.toFixed(0)} bytes a frame, the plain genes ${own.toFixed(0)} (${r.note})`,
    );
    expect(Math.min(reads[0], reads[2], r.bytesPerStep)).toBeLessThanOrEqual(own);
  });
});
