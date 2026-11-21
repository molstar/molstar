/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { type BlockPlacement, DarkTokens, LightTokens, type RingPlacement } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { GenomeViews } from '../views.js';

// The Genome coordinates stories, headless on a fake clock, composed as the stories compose them (a viewer of the view
// and its SVG export, light and dark): the blocks laid out as the story says, the export settled with no issue, its
// blocks' attributes written, the dark SVG another than the light.

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

const view = (id: string) => GenomeViews.find((v) => v.id === id)!;

async function open(id: string) {
  viewer = await openExample(view(id), { gene: '', delayMs: 0, wait: 600 });
  viewer.draw();
  return viewer;
}

/** The view's SVG export in a theme, as the story's twin draws it (the viewer's own export). */
async function exported(v: SeqViewer, dark: boolean) {
  const r = v.renderSvg({ idPrefix: 'g-', tokens: dark ? DarkTokens : LightTokens });
  await vi.advanceTimersByTimeAsync(400);
  return r;
}

describe('Genome coordinates stories, headless', () => {
  it('One assembly in a row: its blocks, each its label, the gaps at most a quarter; exported light and dark', async () => {
    const v = await open('row');
    const p = v.frame!.layout.sections[0].placement as BlockPlacement;
    // chrM (16.6 kb of 3.1 Gb) is under a device pixel at the fit: no strip, its gap with it
    expect(p.strips.n).toBe(24);
    expect(p.labels.slice(0, 3)).toEqual(['chr1', 'chr2', 'chr3']);
    let used = 0;
    for (let i = 0; i < p.strips.n; i++) used += p.strips.axes[i].width;
    expect(used / p.row.width).toBeGreaterThanOrEqual(0.75 - 1e-6);
    const light = await exported(v, false);
    const dark = await exported(v, true);
    expect([light.settled, light.issues]).toEqual([true, []]);
    expect(light.svg).toContain('data-system="/data/toy/(sequence :name chrM)"');
    expect(dark.svg).not.toBe(light.svg);
  });

  it('One assembly on a ring: its title, and the detail following its window with no issue', async () => {
    const v = await open('ring');
    const [ring, detail] = v.frame!.layout.sections;
    const p = ring.placement as RingPlacement;
    expect(p.title).toBe('Toy38');
    expect(p.subtitle).toMatch(/Gb$/);
    v.view.setRange('/views/main/ring', 1e8, 4e8);
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const d = v.frame!.layout.sections[1];
    expect([d.axis!.u0, d.axis!.u1]).toEqual([1e8, 4e8]);
    expect(v.layoutIssues.value.map((i) => i.code)).not.toContain('seq.render.follow.blocks');
    void detail;
    const r = await exported(v, true);
    expect(r.issues).toEqual([]);
  });

  it('Two assemblies on one ring: two groups, their names placed', async () => {
    const v = await open('two-on-a-ring');
    const p = v.frame!.layout.sections[0].placement as RingPlacement;
    expect(p.groupN).toBe(2);
    expect(p.title).toBe('peach · grape');
    const r = await exported(v, false);
    expect(r.svg).toContain('data-assembly="/data/peach"');
    expect(r.svg).toContain('data-assembly="/data/grape"');
  });

  it('B in A’s order: the inverted block reversed, chr4 at four times its length', async () => {
    const v = await open('b-in-a-order');
    const b = v.frame!.layout.sections[1].placement as BlockPlacement;
    expect(b.labels).toEqual(['chr1', 'chr2', 'chr2 inv ←', 'chr2', 'chr3', 'chr4 ×4']);
    const inv = b.strips.axes[2];
    // its left edge shows the inversion's end, its right edge its start: it counts down
    const [left, right] = [inv.positionAt(inv.x0 + 0.01), inv.positionAt(inv.x0 + inv.width - 0.01)];
    expect(left > 69.9e6 && left < 70e6 && right >= 40e6 && right < 40.1e6).toBe(true);
    const r = await exported(v, false);
    expect(r.svg).toContain('data-direction="reversed" data-scale="1"');
    expect(r.svg).toContain('data-scale="4"');
  });

  it('3,000 contigs: runs on the ring and in the row', async () => {
    const v = await open('contigs');
    const [ring, row] = v.frame!.layout.sections;
    expect((ring.placement as RingPlacement).runN).toBeGreaterThan(0);
    const p = row.placement as BlockPlacement;
    expect(p.runN).toBeGreaterThan(0);
    // a run's blocks have no labels of their own; most of the contigs are under a device pixel: not read
    expect(p.strips.n).toBeLessThan(400);
    const r = await exported(v, true);
    expect(r.settled).toBe(true);
  });

  it('every story’s view renders headless with no issue (the export alone, a context of its own)', async () => {
    const ctx = createSeqContext({});
    for (const g of GenomeViews) {
      const r = await ctx.renderSvg(g.build(), { size: { width: 900, height: 'auto' } });
      expect([g.id, r.issues]).toEqual([g.id, []]);
    }
    ctx.retire();
  });
});
