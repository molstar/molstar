/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Plasmid ring and the Complex's blocks stories, headless on a fake clock, composed as the stories compose them (a
// viewer of the view and its SVG export): the ring drawn through its warp with its window over its detail's range, the
// detail following the ring's window across the origin; a hit on the ring names the feature there; the export's lists
// are the canvas's but for the window; the blocks drawn on their own systems.

import {
  BlockPlacement,
  bundledTextMetrics,
  DarkTokens,
  hitView,
  type RingPlacement,
  recordDrawLists,
} from '@molstar/seqx-core/render';
import { ComplexViews, PlasmidViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

describe('Plasmid ring, headless', () => {
  it('the ring and its detail: the window the detail’s range, over the origin in two pieces; a hit names the feature', async () => {
    viewer = await openExample(PlasmidViews.find((v) => v.id === 'ring')!, { gene: '' });
    const v = viewer;
    v.view.setRange('/views/main/ring', 2600, 2760);
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const [ring, detail] = v.frame!.layout.sections;
    const p = ring.placement as RingPlacement;
    expect([p.windowShown, p.windowActive, p.windowStart, p.windowEnd]).toEqual([true, true, 2600, 2760]);
    expect(detail.axis!.pieces).toEqual([2600, 2686, 0, 74]);
    // bla runs over the origin on the minus ring: a hit there at position 0 names it
    const keys = new Set<string>();
    for (let y = p.strips.bandY[1]; y < p.strips.bandY[1] + p.strips.bandH[1]; y += 2) {
      p.warp.map(p.frameX(0.5), y);
      const k = hitView(v.frame!.layout, p.warp.px, p.warp.py)?.element?.key;
      if (k !== undefined) keys.add(String(k));
    }
    expect([...keys]).toContain('bla');
    // the export draws the ring and its labels, not its window (an exported map is a figure)
    const svg = v.renderSvg({ idPrefix: 'p-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(r.issues).toEqual([]);
    // the ring's window is paths through the warp; the detail's rail, idle, is rects
    expect(r.lists.passes.window.paths.n).toBe(0);
    expect(r.lists.passes.labels.texts.n).toBeGreaterThan(10);
  });
});

describe('Plasmid ring, zoomed out past its fit', () => {
  it('a zoom out over the ring at its fit stays the fit; then `=` three times on the detail: both drawn', async () => {
    viewer = await openExample(PlasmidViews.find((v) => v.id === 'ring')!, { gene: '' });
    const v = viewer;
    v.input.mac = false;
    const [ring, detail] = v.frame!.layout.sections;
    const p = ring.placement as RingPlacement;
    p.warp.map(p.frameX(500), p.strips.bandY[0] + 3);
    v.input.wheel(p.warp.px, p.warp.py, { deltaX: 0, deltaY: 300, ctrlKey: true, deltaMode: 0 });
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    // a window as wide as the ring is the fit, never a whole turn from elsewhere
    expect(v.view.get('/views/main/ring')).toMatchObject({ start: 0, end: 2686 });
    v.input.click(v.dataX + 100, detail.y + detail.height / 2);
    for (let k = 0; k < 3; k++) {
      expect(v.input.key({ key: '=' })).toBe(true);
      await vi.advanceTimersByTimeAsync(400);
      v.draw();
      for (const id of ['/views/main/ring', '/views/main/detail']) {
        const w = v.view.get(id)!;
        expect(Number.isFinite(w.start) && w.end > w.start && w.end - w.start <= 2686).toBe(true);
      }
      const f = v.frame!;
      expect(f.layout.sections.every((s) => s.placement.strips.n > 0)).toBe(true);
      expect(f.lists.passes.data.opCount).toBeGreaterThan(40);
    }
    const svg = v.renderSvg({ idPrefix: 'z-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(r.issues).toEqual([]);
    expect(r.lists.passes.data.opCount).toBeGreaterThan(40);
  });
});

describe('Complex blocks, headless', () => {
  it('the row: a strip a block on its own system, the copies merged ×2', async () => {
    viewer = await openExample(ComplexViews.find((v) => v.id === 'blocks')!, { gene: '' });
    const sec = viewer.frame!.layout.sections[0];
    const p = sec.placement as BlockPlacement;
    expect(p.strips.axes.map((a) => a.system)).toEqual(['/data/cdk2', '/data/cyclin', '/data/atp']);
    expect(p.labels).toEqual(['CDK2 (A) ×2', 'cyclin A (B) ×2', 'ATP ×2']);
  });

  // a one-position block is a chip: its label's width at the fit view and zoomed, its label whole, canvas and SVG
  const chips = [
    { id: 'blocks', labels: ['ATP ×2'], zoom: [500, 559] },
    { id: 'blocks-split', labels: ['ATP', 'ATP'], zoom: [1050, 1118] },
    { id: 'blocks-stacked', labels: ['ATP ×2'], zoom: [100, 160] },
    { id: 'blocks-stacked-split', labels: ['ATP', 'ATP'], zoom: [100, 160] },
    { id: 'blocks-ring', labels: ['ATP ×2'] },
  ] as const;
  for (const c of chips)
    for (const dark of [false, true])
      it(`${c.id} (${dark ? 'dark' : 'light'}): the ATP chips at their label’s width, never clipped`, async () => {
        viewer = await openExample(ComplexViews.find((v) => v.id === c.id)!, { gene: '' });
        const v = viewer;
        if (dark) v.setTokens(DarkTokens);
        const text = bundledTextMetrics();
        const windows = 'zoom' in c ? [undefined, c.zoom] : [undefined];
        for (const w of windows) {
          if (w) v.view.setRange('/views/main/complex', w[0], w[1]);
          await vi.advanceTimersByTimeAsync(50);
          v.draw();
          const sec = v.frame!.layout.sections[0];
          const p = sec.placement as BlockPlacement | RingPlacement;
          const shown = p.strips.axes.map((_, i) => i).filter((i) => p.strips.chip[i]);
          expect(shown.map((i) => p.strips.axes[i].system)).toEqual(c.labels.map(() => '/data/atp'));
          if (p instanceof BlockPlacement) expect(shown.map((i) => p.labels[i])).toEqual([...c.labels]);
          for (const i of shown)
            expect(p.strips.axes[i].width).toBeCloseTo(text.measure('12px sans', c.labels[0]) + 10, 6);
          if (!(p instanceof BlockPlacement)) continue;
          // the rail skips the chips as it skips the gaps: the section's units across it, its window the view's
          const rail = sec.rail;
          if (rail && w) expect([rail.total, rail.start, rail.end]).toEqual([sec.axis!.total, w[0], w[1]]);
          if (c.id === 'blocks') expect(rail).toBeDefined();
          // each chip's label drawn in its pill, its width within the chip's
          const texts = recordDrawLists(v.frame!.lists, ['data'])
            .split('\n')
            .filter((l) => l.startsWith(`text "${c.labels[0]}"`));
          expect(texts.length).toBe(shown.length);
          for (const l of texts)
            expect(Number(/ width ([\d.]+)/.exec(l)![1])).toBeLessThanOrEqual(Number(/ max ([\d.]+)/.exec(l)![1]));
        }
        const svg = v.renderSvg({ idPrefix: 'c-' });
        await vi.advanceTimersByTimeAsync(400);
        const r = await svg;
        expect(r.issues).toEqual([]);
        expect(r.svg).toContain(`>${c.labels[0]}<`);
      });
});
