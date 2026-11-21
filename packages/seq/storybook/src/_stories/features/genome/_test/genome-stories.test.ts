/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { type BlockPlacement, DarkTokens, LightTokens, type RingPlacement } from '@molstar/seqx-core/render';
import { actionsShown } from '@molstar/seqx-core/render/headers-draw';
import { Describer } from '@molstar/seqx-viewer/describe';
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

  it('One assembly in a row: a hover over the ruler names the block and its position, at the fit and zoomed', async () => {
    const v = await open('row');
    const s = v.frame!.layout.sections[0];
    const r = s.rail!;
    const p = s.placement as BlockPlacement;
    // chr3's middle, over the ruler where the rail's grab zone is
    const a = p.strips.axes[2];
    const x = a.x0 + a.width / 2;
    // the ruler band's foot: inside the rail's grab zone too
    const y = Math.min(p.strips.bandY[2 * p.strips.bands] + p.strips.bandH[2 * p.strips.bands] - 2, r.bottom - 2);
    expect(y).toBeGreaterThanOrEqual(r.top);
    // at the fit the rail stays a target: its body a place (as outside a window), never the whole window lit
    const fit = v.hit(x, y)!;
    expect([fit.target, fit.part, fit.system, fit.block]).toEqual([
      'window',
      'outside',
      '/data/toy/(sequence :name chr3)',
      2,
    ]);
    expect(fit.position! > 0 && fit.position! < 230e6).toBe(true);
    v.input.move(x, y);
    v.draw();
    expect(v.hover.current).toMatchObject({ target: 'window', system: '/data/toy/(sequence :name chr3)', block: 2 });
    expect(v.emphasis.hover?.part).toBe('outside');
    // drawn: the position mark (1.5 px across the line at the pointer), the rail never grown
    const rects = v.frame!.lists.passes.window.rects;
    const grownPx = v.context.renderConfig.navigator.rail.grownPx;
    let marks = 0;
    for (let i = 0; i < rects.n; i++) {
      expect(rects.h[i]).not.toBe(grownPx);
      if (rects.w[i] === 1.5 && Math.abs(rects.x[i] + 0.75 - x) < 1e-6) marks++;
    }
    expect(marks).toBe(1);
    // a sweep across it at the fit zooms to the swept region
    const x2 = p.strips.axes[4].x0 + p.strips.axes[4].width / 2;
    const u0 = v.view.get('/views/main/genome')!;
    expect(v.input.down(x, y)).toBe(true);
    v.input.move((x + x2) / 2, y);
    v.input.move(x2, y);
    v.input.up();
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const swept = v.view.get('/views/main/genome')!;
    const b2 = s.section.blocks![2];
    const b4 = s.section.blocks![4];
    expect(swept.end - swept.start).toBeLessThan((u0.end - u0.start) / 2);
    expect(swept.start).toBeGreaterThan(b2.offset!);
    expect(swept.end).toBeLessThan(b4.offset! + b4.axis.ranges![1]);
    v.view.fit('/views/main/genome');
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    // zoomed in: the rail is a target, and it names the block's system and position too, never the section's units
    v.view.setRange('/views/main/genome', 4e8, 9e8);
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const r2 = v.frame!.layout.sections[0].rail!;
    const h = v.hit(r2.x + r2.width / 2, r2.bottom - 2)!;
    expect(h.navigator).toBe('rail');
    expect(h.system).toMatch(/^\/data\/toy\/\(sequence :name chr\d+\)$/);
    expect(h.position!).toBeLessThan(260e6);
  });

  it('a pointer sweep along a ruler’s row keeps its header’s actions shown on every step (gaps, the rail, its header)', async () => {
    for (const id of ['row', 'b-in-a-order', 'square']) {
      const v = await open(id);
      const s = v.frame!.layout.sections.find((x) => x.tracks.some((t) => t.track.address.endsWith('/ruler')))!;
      const t = s.tracks.find((x) => x.track.address.endsWith('/ruler'))!;
      const config = v.context.renderConfig;
      for (const y of [
        t.box.y + 1,
        t.box.y + t.box.height / 2,
        t.box.y + t.box.height - 1,
        ...(s.rail ? [s.rail.y] : []),
      ]) {
        let steps = 0;
        for (let x = 2; x < v.width - 2; x += 3) {
          v.input.move(x, y);
          v.draw();
          expect(actionsShown(t.track.address, s.id, v.emphasis, config, true), `${id} ${x},${y}`).toBe(true);
          steps++;
        }
        expect(steps).toBeGreaterThan(100);
      }
      // off the row: they go
      v.input.move(5, s.y + s.height + 40);
      v.draw();
      expect(actionsShown(t.track.address, s.id, v.emphasis, config, true)).toBe(false);
      v.input.leave();
      v.retire();
      viewer = undefined;
    }
  });

  it('One assembly in a row: a double-click on a block’s label fits it, again the whole', async () => {
    const v = await open('row');
    const s = v.frame!.layout.sections[0];
    const p = s.placement as BlockPlacement;
    const a = p.strips.axes[2];
    const x = a.x0 + a.width / 2;
    const y = p.strips.top[2] + 4;
    const total = s.axis!.u1;
    const chr3 = s.section.blocks![2];
    v.input.doubleClick(x, y);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const fitted = v.frame!.layout.sections[0].axis!;
    expect([fitted.u0, fitted.u1]).toEqual([chr3.offset, chr3.offset! + chr3.axis.ranges![1]]);
    v.input.doubleClick(v.frame!.layout.sections[0].axis!.x0 + 20, y);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    expect(v.frame!.layout.sections[0].axis!.u1 - v.frame!.layout.sections[0].axis!.u0).toBe(total);
  });

  it('One assembly in a row: a double-click on a block’s ruler fits the block, again the whole (never one base)', async () => {
    const v = await open('row');
    const s = v.frame!.layout.sections[0];
    const p = s.placement as BlockPlacement;
    const a = p.strips.axes[2];
    const x = a.x0 + a.width / 2;
    // on the ruler's band: a position under the point, no item
    const y = p.strips.bandY[0] + 4;
    expect(v.hit(x, y)?.element).toBeUndefined();
    const total = s.axis!.u1 - s.axis!.u0;
    const chr3 = s.section.blocks![2];
    v.input.doubleClick(x, y);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const fitted = v.frame!.layout.sections[0].axis!;
    expect([fitted.u0, fitted.u1]).toEqual([chr3.offset, chr3.offset! + chr3.axis.ranges![1]]);
    const p2 = v.frame!.layout.sections[0].placement as BlockPlacement;
    v.input.doubleClick(fitted.x0 + fitted.width / 2, p2.strips.bandY[0] + 4);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const back = v.frame!.layout.sections[0].axis!;
    expect(back.u1 - back.u0).toBe(total);
  });

  it('B in A’s order: the hover card names a block’s row by its system (a sequence, never a chain)', async () => {
    const v = await open('b-in-a-order');
    const s = v.frame!.layout.sections[1];
    const owner = {};
    v.hover.set(
      { section: s.section.address, system: s.section.blocks![2].axis.system, position: 5e7, block: 2 } as never,
      owner,
    );
    const card = new Describer({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    }).describe(v.hover.current);
    v.hover.clear(owner);
    expect(card?.rows).toContainEqual({ label: 'Sequence', value: 'chr2 inv' });
    expect(card?.rows.some((r) => r.label === 'Chain')).toBe(false);
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
    // each group's name past the arcs and the block labels: its box meets no placed label and stays off the arcs
    // ("peach" and "grape" straddled their arcs among the block labels)
    for (let k = 0; k < p.groupN; k++) {
      const [cx, cy] = [p.groupTX[k], p.groupTY[k]];
      const w = 6.7 * p.groupText[k].length;
      const h = p.titleSize;
      const nx = Math.max(cx - w / 2, Math.min(p.cx, cx + w / 2)) - p.cx;
      const ny = Math.max(cy - h / 2, Math.min(p.cy, cy + h / 2)) - p.cy;
      expect(Math.hypot(nx, ny), p.groupText[k]).toBeGreaterThan(p.r + 3);
      for (let q = 0; q < p.labels.m; q++) {
        const lx = p.labels.right[q] ? p.labels.x[q] : p.labels.x[q] - p.labels.width[q];
        const ly = p.labels.y[q] - 0.35 * p.labelSize;
        const meets =
          lx < cx + w / 2 && lx + p.labels.width[q] > cx - w / 2 && ly - 6 < cy + h / 2 && ly + 6 > cy - h / 2;
        expect(meets, `${p.groupText[k]} over label ${q}`).toBe(false);
      }
    }
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
