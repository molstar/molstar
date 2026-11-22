/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { blockSpan, placeOfUnit } from '@molstar/seqx-core/view';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { GenomeViews, OneAtATime, openSquare } from '../views.js';

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
    expect(n.height).toBe(28 + 6 + 22);
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

  it('PageDown on its slider and Alt+PageDown on the view: the next block by the switching rule (block.next)', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    n.draw();
    const s = v.view.resolved.value!.sections[0];
    const [chr1, chr2, chr3] = s.blocks!;
    const len = (b: typeof chr1) => b.axis.ranges![1];
    // zoomed in: 1 Mb at 30% of chr1 → 1 Mb at 30% of chr2, then of chr3
    const at30 = (b: typeof chr1) => b.offset! + 0.3 * len(b);
    v.view.setRange(Genome, at30(chr1) - 5e5, at30(chr1) + 5e5);
    await settle(v, 100);
    expect(n.key('PageDown')).toBe(true);
    await settle(v);
    const now = () => v.view.get(Genome)!;
    expect(now().end - now().start).toBeCloseTo(1e6, 0);
    expect((now().start + now().end) / 2).toBeCloseTo(at30(chr2), 0);
    expect(v.input.key({ key: 'PageDown', altKey: true })).toBe(true);
    await settle(v);
    expect((now().start + now().end) / 2).toBeCloseTo(at30(chr3), 0);
    // from a block whole: the previous whole
    v.view.setRange(Genome, chr3.offset!, chr3.offset! + len(chr3));
    await settle(v, 100);
    expect(v.runCommand('block.previous', { section: Genome })).toBe('Previous block');
    await settle(v);
    expect([now().start, now().end]).toEqual([chr2.offset, chr2.offset! + len(chr2)]);
  });

  it('a click on another chromosome: whole from whole, the zoom kept at the same relative place, whole past its length', async () => {
    const v = await open('chromosome-at-a-time');
    const n = v.navigatorModel('genome-nav')!;
    const s = v.view.resolved.value!.sections[0];
    const named = (name: string) =>
      s.blocks!.find((b) => v.spec.model.value.systems.get(b.axis.system)?.name === name)!;
    const chr1 = named('chr1');
    const chr21 = named('chr21');
    const len = (b: typeof chr1) => b.axis.ranges![1];
    const now = () => v.view.get(Genome)!;
    const clickOn = async (b: typeof chr1) => {
      n.draw();
      let x = -1;
      for (let px = 1; px < n.width && x < 0; px++) if (n.hit(px, 5)?.system === b.axis.system) x = px + 1;
      expect(n.down(x, 5)).toBe(true);
      n.up();
      expect(n.click(x, 5)).toBe(true);
      await settle(v);
    };
    // a whole chr21 → chr1 whole
    v.view.setRange(Genome, chr21.offset!, chr21.offset! + len(chr21));
    await settle(v, 100);
    await clickOn(chr1);
    expect([now().start, now().end]).toEqual([chr1.offset, chr1.offset! + len(chr1)]);
    // 1 Mb at 30% of chr1 → 1 Mb at 30% of chr21
    const c1 = chr1.offset! + 0.3 * len(chr1);
    v.view.setRange(Genome, c1 - 5e5, c1 + 5e5);
    await settle(v, 100);
    await clickOn(chr21);
    expect(now().end - now().start).toBeCloseTo(1e6, 0);
    expect((now().start + now().end) / 2).toBeCloseTo(chr21.offset! + 0.3 * len(chr21), 0);
    // a window wider than chr21 → chr21 whole
    expect(len(chr21)).toBeLessThan(80e6);
    v.view.setRange(Genome, chr1.offset! + 10e6, chr1.offset! + 90e6);
    await settle(v, 100);
    await clickOn(chr21);
    expect([now().start, now().end]).toEqual([chr21.offset, chr21.offset! + len(chr21)]);
  });

  it('each navigator story exports with its navigator, light and dark, with no issue (the stories’ SVG twins)', async () => {
    for (const [id, nav] of [
      ['chromosome-at-a-time', 'genome-nav'],
      ['ring-switch', 'chromosomes'],
      ['square', 'square'],
    ] as const) {
      const v = await open(id);
      const svgs: string[] = [];
      for (const tokens of [LightTokens, DarkTokens]) {
        const r = v.renderSvg({ idPrefix: 'n-', tokens, navigators: 'include' } as never);
        await vi.advanceTimersByTimeAsync(400);
        const done = await r;
        expect([id, done.settled, done.issues]).toEqual([id, true, []]);
        expect(done.svg).toContain(`data-navigator="${nav}"`);
        svgs.push(done.svg);
      }
      expect(svgs[1]).not.toBe(svgs[0]);
      v.retire();
      viewer = undefined;
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
    // a navigator part bound to a ring section reads as any horizontal one: the window's region at rest, the mark's region
    // while it hovers the mark, the window's again when it leaves
    const rest = n.values.value[0].text;
    expect(rest).toMatch(/^chr1:1 – chr\w+:[\d,]+$/);
    expect(n.readout.value).toBe(rest);
    n.move(60, 5);
    expect(n.readout.value).toBe('chr1 – chrM');
    n.leave();
    expect(n.readout.value).toBe(rest);
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
    // a cell (chr1, chr2): chr1's column, chr2's row (the first y over it)
    let x = 0;
    while (n.hit(x, 150)?.block !== 0) x++;
    x += 2;
    let y = 0;
    while (n.hit(x, y)?.yBlock !== 1) y++;
    y += 2;
    const h = n.hit(x, y)!;
    expect([h.block, h.yBlock]).toEqual([0, 1]);
    expect(n.doubleClick(x, y)).toBe(true);
    await settle(v);
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual([
      0,
      blocks[1].offset! + blocks[1].axis.ranges![1],
    ]);
    n.draw();
    expect(n.doubleClick(x, y)).toBe(true);
    await settle(v);
    expect(v.view.get(Genome)!.end).toBe(total);
  });
});

describe('The story as composed: one chromosome at a time (navigate.acrossBlocks: stop, the zoom strip at 200 px)', () => {
  it('a click on another chromosome on the navigator (its name or its strip) is a jump there: the view and the window on it', async () => {
    const v = await open('chromosome-at-a-time', OneAtATime);
    const n = v.navigatorModel('genome-nav')!;
    n.draw();
    const s = v.view.resolved.value!.sections[0];
    const sys = (k: number) => s.blocks![k].axis.system;
    // the view opens on chr1 (the stop's fit), the zoom strip under the genome, its track alone (no label row)
    expect(v.view.get(Genome)!.end).toBe(s.blocks![0].axis.ranges![1]);
    expect(n.height).toBe(28 + 6 + 22);
    for (const [k, y] of [
      [6, 4],
      [4, 14],
      [9, 6],
    ] as const) {
      // an x on the navigator over block k (its name at the top, its strip lower)
      let x = -1;
      for (let px = 2; px < n.width && x < 0; px++) if (n.hit(px, y)?.system === sys(k)) x = px + 2;
      expect(x, `block ${k}`).toBeGreaterThan(0);
      expect(n.down(x, y)).toBe(true);
      n.up();
      expect(n.click(x, y)).toBe(true);
      await settle(v);
      // centred there at the window's width, kept inside that chromosome (the stop)
      const now = v.view.get(Genome)!;
      const b = s.blocks![k];
      expect(now.start, `block ${k}`).toBeGreaterThanOrEqual(b.offset!);
      expect(now.end, `block ${k}`).toBeLessThanOrEqual(b.offset! + blockSpan(b));
      const f = n.draw()!;
      expect([f.windows[0].s0, f.windows[0].s1]).toEqual([sys(k), sys(k)]);
    }
  });
});

/** Every text a navigator's frame draws whose box leaves the navigator (its baseline less the font's size, its ends). */
function textsOutside(f: { lists: import('@molstar/seqx-core/render').DrawLists }, width: number, height: number) {
  const out: string[] = [];
  for (const pass of Object.values(f.lists.passes)) {
    const t = pass.texts;
    for (let i = 0; i < t.n; i++) {
      const font = f.lists.string(t.font[i]);
      const size = Number(/(\d+(?:\.\d+)?)px/.exec(font)?.[1] ?? 12);
      const w = Math.min(t.width[i], t.maxWidth[i]);
      const x0 = t.anchor[i] === 0 ? t.x[i] : t.anchor[i] === 1 ? t.x[i] - w / 2 : t.x[i] - w;
      const top = t.y[i] - size * 0.8;
      const bottom = t.y[i] + size * 0.2;
      if (top < -0.5 || bottom > height + 0.5 || x0 < -0.5 || x0 + w > width + 0.5)
        out.push(`${pass.name}: ${f.lists.string(t.str[i])} at ${t.x[i].toFixed(1)},${t.y[i].toFixed(1)}`);
    }
  }
  return out;
}

describe('A zoom strip over blocks fits its navigator (the strip’s ticks and names inside its box)', () => {
  it('on a chromosome at a time, at several zooms and widths', async () => {
    const v = await open('chromosome-at-a-time', OneAtATime);
    const n = v.navigatorModel('genome-nav')!;
    const s = v.view.resolved.value!.sections[0];
    for (const width of [1200, 700, 360])
      for (const [k, w] of [
        [6, 1],
        [1, 0.1],
        [9, 0.002],
      ] as const) {
        v.resize(width);
        const b = s.blocks![k];
        const span = blockSpan(b) * w;
        v.view.setRange(s.address, b.offset! + span * 0.1, b.offset! + span * 1.1);
        await settle(v, 100);
        const f = n.draw()!;
        expect(n.height, `${k} ${w}`).toBeGreaterThan(28);
        expect(textsOutside(f, n.width, n.height), `${width} ${k} ${w}`).toEqual([]);
      }
  });

  it('on the ring that switches: a chromosome switched in, its zoom strip inside the navigator', async () => {
    const v = await open('ring-switch');
    const n = v.navigatorModel('chromosomes')!;
    for (const name of ['chr5', 'chr7', 'chr21']) {
      await v.execute('section.include', { address: '/views/main/ring', include: [name] });
      await settle(v);
      const f = n.draw()!;
      expect(n.height, name).toBeGreaterThan(28);
      expect(textsOutside(f, n.width, n.height), name).toEqual([]);
    }
  });
});

describe('The square as composed (it opens on chr3)', () => {
  it('names its blocks on both axes, its cells sized by length; a click on another cell centres there; its title names the window', async () => {
    const v = await open('square');
    openSquare(v);
    await settle(v, 100);
    const n = v.navigatorModel('square')!;
    const f = n.draw()!;
    const texts: string[] = [];
    const t = f.lists.passes.labels.texts;
    for (let i = 0; i < t.n; i++) texts.push(f.lists.string(t.str[i]));
    // along y the names whole, along x the shared prefix dropped where a column is narrow
    for (let k = 1; k <= 12; k++) expect(texts).toContain(`chr${k}`);
    expect(texts).toContain('1');
    const s = v.view.resolved.value!.sections[0];
    const chr3 = s.blocks![2];
    n.setFocus(true);
    n.setFocus(false);
    expect(n.values.value[0].text).toBe(`chr3: 1–${chr3.axis.ranges![1].toLocaleString('en-US')} of Toy38`);
    // its readout is short for its narrow part: the block and its range, no assembly; fixed at rest
    expect(n.readout.value).toBe(`chr3: 1–${chr3.axis.ranges![1].toLocaleString('en-US')}`);
    expect(n.readout.value!.length).toBeLessThan(28);
    // cells by length: chr1's column wider than chr12's
    const col = (k: number) => {
      const xs: number[] = [];
      for (let px = 0; px < n.width; px += 0.25) if (n.hit(px, 150)?.block === k) xs.push(px);
      return xs[xs.length - 1] - xs[0];
    };
    expect(col(0)).toBeGreaterThan(col(11) * 1.5);
    // a click on chr8's column, chr1's row (outside the window): the view centred on it, one step back
    let x = -1;
    for (let px = 0; px < n.width; px++) if (n.hit(px, 18)?.block === 7) x = px;
    expect(n.hit(x, 18)!.part).toBe('outside');
    expect(n.down(x, 18)).toBe(true);
    n.up();
    expect(n.click(x, 18)).toBe(true);
    await settle(v);
    const now = v.view.get(Genome)!;
    expect(placeOfUnit(s, (now.start + now.end) / 2)?.system).toBe(s.blocks![7].axis.system);
    expect(v.view.back()).toBe(true);
  });
});

describe('The square’s sweep (a drag from a cell outside the window)', () => {
  /** The rects of the window pass: the window, then a sweep's. */
  const windowRects = (n: NonNullable<ReturnType<SeqViewer['navigatorModel']>>) => {
    const r = n.draw()!.lists.passes.window.rects;
    return Array.from({ length: r.n }, (_, i) => ({ x: r.x[i], y: r.y[i], w: r.w[i], h: r.h[i] }));
  };

  async function opened() {
    const v = await open('square');
    openSquare(v);
    await settle(v, 100);
    const n = v.navigatorModel('square')!;
    n.draw();
    return { v, n };
  }

  it('draws the swept rectangle live on the diagonal, snapped to whole cells; release fits exactly those', async () => {
    const { v, n } = await opened();
    expect(windowRects(n)).toHaveLength(1);
    const px = (k: number) => {
      const xs: number[] = [];
      for (let i = 0; i < n.width; i += 0.25) if (n.hit(i, 18)?.block === k) xs.push(i);
      return xs;
    };
    const c7 = px(7);
    const c8 = px(8);
    expect(c7.length).toBeGreaterThan(0);
    // pressed inside chr8's column (off the window), ending mid-way over chr9's
    const x = c7[Math.floor(c7.length / 2)];
    const end = c8[Math.floor(c8.length / 2)];
    expect(n.hit(x, 18)!.part).toBe('outside');
    expect(n.down(x, 18)).toBe(true);
    n.move(end, 60);
    const [win, sweep] = windowRects(n);
    expect(win).toBeDefined();
    // both cells whole: from chr8's left edge to chr9's right edge
    expect(sweep.x).toBeCloseTo(c7[0], 0);
    expect(sweep.x + sweep.w).toBeCloseTo(c8[c8.length - 1], 0);
    expect(sweep.w).toBeGreaterThan(end - x);
    // one section on both axes: the same range on y, so its square sits on the diagonal
    expect(sweep.y).toBeGreaterThan(0);
    expect(sweep.h).toBeGreaterThan(0);
    n.up();
    await settle(v);
    const s = v.view.resolved.value!.sections[0].blocks!;
    const now = v.view.get(Genome)!;
    expect(now.start).toBeCloseTo(s[7].offset!, 3);
    expect(now.end).toBeCloseTo(s[8].offset! + blockSpan(s[8]), 3);
  });

  it('Escape cancels it: no sweep drawn, the view unchanged', async () => {
    const { v, n } = await opened();
    const before = [v.view.get(Genome)!.start, v.view.get(Genome)!.end];
    let x = -1;
    for (let px = 0; px < n.width; px++) if (n.hit(px, 18)?.block === 7) x = px;
    expect(n.down(x, 18)).toBe(true);
    n.move(x + 30, 60);
    expect(windowRects(n)).toHaveLength(2);
    expect(n.key('Escape')).toBe(true);
    expect(windowRects(n)).toHaveLength(1);
    n.up();
    await settle(v);
    expect([v.view.get(Genome)!.start, v.view.get(Genome)!.end]).toEqual(before);
  });
});
