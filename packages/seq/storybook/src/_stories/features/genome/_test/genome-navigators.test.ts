/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { GenomeViews } from '../views.js';

// The genome navigators' stories, headless on a fake clock as the stories compose them (a viewer with its navigator,
// and the view's SVG export light and dark): a navigator over blocks (its window, zoom strip, readout and slider text,
// a click a jump, a double-click a block then the whole), `navigate.acrossBlocks: stop`, `block.next`, and a ring's
// navigator with `switch: include` (EC-184).

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const view = (id: string) => GenomeViews.find((v) => v.id === id)!;
const Genome = '/views/main/genome';

async function open(id: string, config?: Parameters<typeof openExample>[1]['config']) {
  viewer = await openExample(view(id), { gene: '', delayMs: 0, wait: 600, ...(config ? { config } : {}) });
  viewer.draw();
  return viewer;
}

async function settle(v: SeqViewer, ms = 1000) {
  await vi.advanceTimersByTimeAsync(ms);
  v.draw();
}

describe('A navigator over blocks (EC-184)', () => {
  it('draws the genome’s blocks, its window the whole picture at the fit, with no issue', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    const f = n.draw()!;
    expect(n.issues).toEqual([]);
    const w = f.windows[0];
    expect(w.shown).toBe(true);
    expect([w.x0, w.x1]).toEqual([0, n.width]);
    expect(n.values.value[0].text).toBe('chr1:1 – chrM:16,569');
  });

  it('a window inside one chromosome: the zoom strip of that chromosome, the slider naming it and the assembly', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    const s = v.view.resolved.value!.sections[0];
    const chr2 = s.blocks![1];
    const at = chr2.offset! + 147_690_000;
    v.view.setRange(Genome, at, at + 10_000);
    await settle(v, 100);
    const h0 = n.height;
    const f = n.draw()!;
    // the zoom strip below the box: the window a few kb of a 3 Gb genome
    expect(n.height).toBeGreaterThan(h0 - 1);
    expect(n.height).toBe(28 + 16);
    const w = f.windows[0];
    expect([w.s0, w.s1]).toEqual([chr2.axis.system, chr2.axis.system]);
    expect([w.p0, w.p1]).toEqual([147_690_000, 147_700_000]);
    n.setFocus(true);
    n.setFocus(false);
    expect(n.values.value[0].text).toBe('chr2: 147,690,001–147,700,000 of Toy38');
  });

  it('the readout names the block and the position there; a click on another block is a jump there, back by Alt+Left', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    n.draw();
    const s = v.view.resolved.value!.sections[0];
    v.view.setRange(Genome, 0, 50e6);
    await settle(v, 100);
    n.draw();
    // over chr8 on the picture: its name and a position on it
    const x = 300;
    n.move(x, 5);
    expect(n.readout.value).toMatch(/^chr\d+: [\d,]+$/);
    const name = n.readout.value!.split(':')[0];
    expect(n.down(x, 5)).toBe(true);
    n.up();
    expect(n.click(x, 5)).toBe(true);
    await settle(v);
    const now = v.view.get(Genome)!;
    const c = (now.start + now.end) / 2;
    const block = s.blocks!.find((b) => c >= b.offset! && c < b.offset! + b.axis.ranges![1])!;
    expect(v.spec.model.value.systems.get(block.axis.system)?.name).toBe(name);
    expect(v.view.back()).toBe(true);
    await settle(v);
    expect(v.view.get(Genome)!.start).toBe(0);
  });

  it('a double-click fits the block under the pointer, again the whole (toggle-fit)', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    n.draw();
    const s = v.view.resolved.value!.sections[0];
    const total = s.axis.ranges ? s.axis.ranges[1] : v.view.get(Genome)!.end;
    expect(n.doubleClick(20, 5)).toBe(true);
    await settle(v);
    const chr1 = s.blocks![0];
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual([0, chr1.axis.ranges![1]]);
    n.draw();
    expect(n.doubleClick(20, 5)).toBe(true);
    await settle(v);
    expect(v.view.get(Genome)!.end).toBe(total);
  });

  it('PageDown on its slider and Alt+PageDown on the view: the next block, fitted (block.next)', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    n.draw();
    const s = v.view.resolved.value!.sections[0];
    const [chr1, chr2, chr3] = s.blocks!;
    v.view.setRange(Genome, 1e6, 2e6);
    await settle(v, 100);
    expect(n.key('PageDown')).toBe(true);
    await settle(v);
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual([
      chr2.offset,
      chr2.offset! + chr2.axis.ranges![1],
    ]);
    expect(v.input.key({ key: 'PageDown', altKey: true })).toBe(true);
    await settle(v);
    expect(v.view.get(Genome)!.start).toBe(chr3.offset);
    expect(v.runCommand('block.previous', { section: Genome })).toBe('Previous block');
    await settle(v);
    expect(v.view.get(Genome)!.start).toBe(chr2.offset);
    void chr1;
  });

  it('exports with its navigator, light and dark, with no issue', async () => {
    const v = await open('chromosome-at-a-time');
    for (const tokens of [LightTokens, DarkTokens]) {
      const r = v.renderSvg({ idPrefix: 'n-', tokens, navigators: 'include' } as never);
      await vi.advanceTimersByTimeAsync(400);
      const done = await r;
      expect([done.settled, done.issues]).toEqual([true, []]);
      expect(done.svg).toContain('data-navigator="genome-nav"');
    }
  });
});

describe('navigate.acrossBlocks: stop', () => {
  it('a pan stops at the block’s end, a zoom out at the block, a fit is the block; a jump moves to another', async () => {
    const v = await open('chromosome-at-a-time', [[viewerItem('navigate.acrossBlocks'), 'stop']]);
    const s = v.view.resolved.value!.sections[0];
    const [chr1, chr2] = s.blocks!;
    const end1 = chr1.axis.ranges![1];
    v.view.setRange(Genome, end1 - 20e6, end1 - 10e6);
    await settle(v, 100);
    v.view.pan(Genome, 50e6);
    await settle(v, 100);
    expect(v.view.get(Genome)!.end).toBeCloseTo(end1, 0);
    v.view.zoom(Genome, 1 / 1000);
    await settle(v, 100);
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual([0, end1]);
    // a jump into chr2 moves there, and a fit is chr2 then
    v.view.setRange(Genome, chr2.offset! + 1e6, chr2.offset! + 2e6);
    await settle(v, 100);
    v.view.fit(Genome);
    await settle(v, 100);
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual([
      chr2.offset,
      chr2.offset! + chr2.axis.ranges![1],
    ]);
  });
});

describe('A ring with navigator.switch: include', () => {
  it('a click on a chromosome shows it alone on the ring (one undo step); Home shows every one again', async () => {
    const v = await open('ring-switch');
    const n = v.navigatorModel('chromosomes')!;
    n.draw();
    expect(n.issues).toEqual([]);
    const ringNode = () =>
      v.spec.root.value.children!.find((c) => c.kind === 'views')!.children![0].children!.find((c) => c.id === 'ring')!;
    // chr2's place on the navigator
    const x = 60;
    const name = v.spec.model.value.systems.get(n.hit(x, 5)!.system!)!.name;
    expect(name).toBe('chr2');
    expect(n.down(x, 5)).toBe(true);
    n.up();
    expect(n.click(x, 5)).toBe(true);
    await settle(v);
    expect(ringNode().params?.include).toEqual([name]);
    expect(v.view.resolved.value!.sections[0].blocks!.length).toBe(1);
    // a click on what it already shows alone: no change
    n.draw();
    expect(n.down(x, 5)).toBe(true);
    n.up();
    expect(n.click(x, 5)).toBe(false);
    expect(n.key('Home')).toBe(true);
    await settle(v);
    expect(ringNode().params?.include).toBeUndefined();
    expect(v.view.resolved.value!.sections[0].blocks!.length).toBe(25);
  });
});

describe('The square over blocks × blocks', () => {
  it('draws the assembly’s pairs as a lattice, the window on the diagonal; a double-click fits a cell’s pair, again the whole', async () => {
    const v = await open('square');
    const n = v.navigatorModel('square')!;
    const f = n.draw()!;
    expect(n.issues).toEqual([]);
    expect([n.width, n.height]).toEqual([220, 220]);
    // 13 lines each way for 12 chromosomes, the window one rectangle
    expect(f.lists.passes.data.size).toBeGreaterThanOrEqual(26);
    expect(f.lists.passes.window.size).toBe(1);
    const s = v.view.resolved.value!.sections[0];
    const blocks = s.blocks!;
    const total = blocks[11].offset! + blocks[11].axis.ranges![1];
    // a cell (chr1, chr2): left, a few px below chr1's row
    const y = 1 + (218 * blocks[1].offset!) / total + 3;
    const h = n.hit(5, y)!;
    expect([h.block, h.yBlock]).toEqual([0, 1]);
    expect(n.doubleClick(5, y)).toBe(true);
    await settle(v);
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual([
      0,
      blocks[1].offset! + blocks[1].axis.ranges![1],
    ]);
    n.draw();
    expect(n.doubleClick(5, y)).toBe(true);
    await settle(v);
    expect(v.view.get(Genome)!.end).toBe(total);
  });
});
