/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gridCellHover } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openHic } from '../openings.js';
import { advance, drawn, exported, imageColours, openRegion } from './region-twin.js';

// Features/Region reads/Hi-C, headless on a fake clock (the story's composition): a map by region, its tiles arriving
// coarse to fine, a pending one drawn from the coarser held, and over the genome at the fit read through the
// whole-assembly level (EC-169, EC-223).

let viewer: SeqViewer | undefined;
// the genome chunk landed first: its import is real I/O, never on the fake clock
beforeAll(async () => {
  await loadGenome();
});
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

interface GridLay {
  readonly images: { readonly length: number };
  readonly pend: { readonly length: number };
  readonly binLevel: number;
  readonly tiles: { readonly length: number };
  readonly system: string;
  readonly ySystem: string;
  readonly region?: unknown;
}

/** The plot's layer (the grid section's one strip). */
const plot = (v: SeqViewer) => v.frame!.layout.sections[1].tracks[0].layers[0].lay as unknown as GridLay;
const box = (v: SeqViewer) => JSON.stringify(v.frame!.layout.sections[1].tracks[0].box);

describe('Features/Region reads/Hi-C, headless', () => {
  it('by region (EC-169, EC-223): pending, then the coarse tile, then sharp; the plot’s box fixed throughout', async () => {
    const v = (viewer = await openRegion('hic-region'));
    await advance(v, 50);
    const box0 = box(v);
    openHic(v);
    // asked, nothing held: a pending rect, never an empty plot
    await advance(v, 16);
    expect(plot(v).images.length).toBe(0);
    expect(plot(v).pend.length).toBeGreaterThan(0);
    expect(box(v)).toBe(box0);
    // the coarse level (a third of the delay) lands: its tile drawn, blocky
    await advance(v, 100);
    expect(plot(v).images.length).toBe(1);
    expect(plot(v).pend.length).toBe(0);
    const coarse = plot(v).binLevel;
    expect(coarse).toBeGreaterThan(1);
    // zoomed in the window asks the finest level: its tile pending, the coarse one held under it, never empty
    await advance(v, 100);
    expect(plot(v).binLevel).toBe(1);
    expect(plot(v).images.length).toBe(1);
    expect(plot(v).tiles.length).toBe(0);
    expect(drawn(v, 'image')).toBe(1);
    // the held tile in the colours of the domain held, never those it was first drawn in (none held yet: saturated)
    expect(imageColours(v)[0]).toBeGreaterThan(5);
    expect(box(v)).toBe(box0);
    // the fine tiles land: sharp
    await advance(v, 2000);
    expect(plot(v).tiles.length).toBeGreaterThan(0);
    expect(plot(v).pend.length).toBe(0);
    expect(imageColours(v).every((n) => n > 5)).toBe(true);
    expect(box(v)).toBe(box0);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
      // an image is runs of rects in the SVG
      expect(svg.svg).toContain('<rect');
    }
    expect((await exported(v, true)).svg).not.toBe((await exported(v, false)).svg);
  });

  it('a cell’s hover names its bin on each axis and its value; the lower triangle is the upper one’s mirror', async () => {
    const v = (viewer = await openRegion('hic-region'));
    openHic(v);
    await advance(v, 3000);
    const lay = plot(v);
    const name = (s: string) => s;
    // positions 120 and 220 of chr1 either way round: the same value
    const above = gridCellHover(lay as never, 120, 220, name)!;
    const below = gridCellHover(lay as never, 220, 120, name)!;
    const value = (c: typeof above) => c.rows.find((r) => r[0] === 'Value')![1];
    expect(value(above)).toMatch(/^\d/);
    expect(value(below)).toBe(value(above));
    expect(above.rows.find((r) => r[0] === 'x')![1]).toMatch(/^\/data\/chr1:\d+/);
    expect(above.rows.find((r) => r[0] === 'y')![1]).toMatch(/^\/data\/chr1:\d+/);
  });

  it('over the genome at the fit (EC-223, EC-191): a block pair a cell of chromosomes, each read by block pair, the lower ones mirrored; in both themes', async () => {
    const v = (viewer = await openRegion('hic-genome'));
    await advance(v, 20_000);
    // six chromosomes: 36 block pairs, each a strip with its tiles drawn (the 15 below the diagonal mirrored)
    const lays = v.frame!.layout.sections[1].lays.map((l) => l[0].lay as unknown as GridLay);
    expect(lays.length).toBe(36);
    for (const lay of lays) {
      expect(lay.images.length, `${lay.system} × ${lay.ySystem}`).toBeGreaterThan(0);
      expect(lay.pend.length).toBe(0);
    }
    expect(new Set(lays.map((l) => `${l.system}|${l.ySystem}`)).size).toBe(36);
    expect(drawn(v, 'image')).toBeGreaterThanOrEqual(36);
    expect(v.store.issues.value).toEqual([]);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
    }
    expect((await exported(v, true)).svg).not.toBe((await exported(v, false)).svg);
  });
});
