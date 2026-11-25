/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type RangesLay, recordDrawLists } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { DefaultChromosome, Tp53 } from '@molstar/seqx-examples/p53/chromosome-data';
import { HoverCardProvider } from '@molstar/seqx-viewer/behaviors';
import { viewerItem } from '@molstar/seqx-viewer/config';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { clickAt } from '../../features/selection/compose.js';
import { type FollowedSvg, followSvg } from '../../follow-svg.js';

// The Chromosome story, headless on a fake clock, composed as `example-root.ts` composes it: the viewer with the
// chromosome's delayed runtime, and its SVG twin exporting through the viewer's store each time the view rests (two
// consumers of one store). The twin is the export of what the canvas shows at the fit and at TP53, and the export's
// draw lists the canvas's; TP53's exons, a node its region brought, stay drawn on the canvas once the exports that
// read them have ended; a double-click on TP53 focuses its span; the hover card names the gene and its span; a click
// selects it, 0 fits again; a drag on the ideogram's window pans the view.

const DelayMs = 150;
const view = P53Views.find((v) => v.id === 'chromosome')!;
const raf = DomHost.requestAnimationFrame;
let viewer: SeqViewer | undefined;
let followed: FollowedSvg | undefined;

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  followed?.dispose();
  followed = undefined;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const strip = (x: string) => x.replace(/t\d+-/g, '');

async function open() {
  const base = DefaultSeqViewerSpec();
  const v = (viewer = createSeqViewer({
    ...base,
    ui: undefined,
    // the export has no rail row, so the twin compares with none drawn beside the ideogram (the old replacement rule)
    config: [[viewerItem('navigator.positionsRail'), 'unlessNavigator']] as never,
    registry: [...(base.registry ?? []), ...(view.registry?.({ delayMs: DelayMs }) ?? [])],
  }));
  v.resize(720);
  const host = { innerHTML: '' } as unknown as HTMLElement;
  const svg = (followed = followSvg(v, host, { prefix: 't', svg: { navigators: 'include' } }));
  void v
    .load(view.build())
    .then(() => v.settle())
    .then(() => svg.draw());
  await vi.advanceTimersByTimeAsync(3000);
  const fresh = async () =>
    strip(
      (
        await Promise.all([
          v.renderSvg({ idPrefix: 't1-', navigators: 'include' }),
          vi.advanceTimersByTimeAsync(10 * DelayMs),
        ])
      )[0].svg,
    );
  // the export's draw lists against the canvas's, the data pass (the navigators are the export's alone)
  const lists = async () => {
    const r = (await Promise.all([v.renderSvg(), vi.advanceTimersByTimeAsync(10 * DelayMs)]))[0];
    v.draw();
    return { canvas: recordDrawLists(v.frame!.lists, ['data']), svg: recordDrawLists(r.lists, ['data']) };
  };
  return { v, host, fresh, lists };
}

/** A track's first layer as laid out on the canvas, by its label. */
const layerOf = (v: SeqViewer, label: string) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === label)?.layers[0].lay as RangesLay | undefined;

/** The centre of TP53's mark on the genes' track. */
function tp53Point(v: SeqViewer): [number, number] {
  const lay = layerOf(v, 'Genes')!;
  const k = Array.from({ length: lay.hn }, (_, j) => j).find((j) => lay.keys[lay.hitem[j]] === 'TP53')!;
  expect(k, 'TP53 hit box').toBeGreaterThanOrEqual(0);
  return [(lay.hx0[k] + lay.hx1[k]) / 2, lay.laneY + (lay.hlane[k] + 0.5) * lay.laneHeight];
}

const Section = '/views/main/chr17';
/** About TP53: 35 kb, its gene a bar. */
const Around = { start: Tp53.start - 8000, end: Tp53.end + 8000 };

describe('p53’s chromosome, the story headless', () => {
  it('the canvas and its SVG twin show the same view at the fit and at TP53; its exons stay once the exports end', async () => {
    const { v, host, fresh, lists } = await open();
    const r = v.view.get(Section)!;
    expect([r.start, r.end]).toEqual([0, DefaultChromosome.length]);
    expect(layerOf(v, 'Genes')!.mode).toBe('density');
    expect(strip(host.innerHTML)).toBe(await fresh());
    const fit = await lists();
    expect(fit.canvas.length).toBeGreaterThan(0);
    expect(fit.svg).toBe(fit.canvas);

    v.view.setRange(Section, Around.start, Around.end);
    // the region answers, TP53's exons arrive as a node; the view rests and the twin exports through the same store
    await vi.advanceTimersByTimeAsync(3000);
    const exons = layerOf(v, 'Exons TP53')!;
    expect(exons.pending).toEqual([]);
    expect(exons.m).toBeGreaterThan(0);
    expect(strip(host.innerHTML)).toBe(await fresh());
    expect(strip(host.innerHTML)).toContain('TP53');
    const at = await lists();
    expect(at.svg).toBe(at.canvas);
    // every export has ended: the canvas still draws TP53's exons, nothing pending
    await vi.advanceTimersByTimeAsync(5000);
    v.draw();
    const after = layerOf(v, 'Exons TP53')!;
    expect(after.pending).toEqual([]);
    expect(after.m).toBeGreaterThan(0);
  });

  it('a double-click on TP53 focuses its span; the hover card names it and its span', async () => {
    const { v } = await open();
    v.view.setRange(Section, Around.start, Around.end);
    await vi.advanceTimersByTimeAsync(3000);
    v.draw();

    // the card says an item's label, extent and strand (TP53 is on the minus strand)
    v.input.move(...tp53Point(v));
    await vi.advanceTimersByTimeAsync(600);
    const card = v.behavior(HoverCardProvider)!.shown.value!.card;
    expect(card.title).toBe('TP53');
    expect(card.rows).toContainEqual({ label: 'Extent', value: 'chr17 7,668,402–7,687,550' });
    expect(card.rows).toContainEqual({ label: 'Strand', value: '−' });

    v.input.doubleClick(...tp53Point(v));
    await vi.advanceTimersByTimeAsync(1000);
    const r = v.view.get(Section)!;
    expect(r.start).toBeCloseTo(Tp53.start, 0);
    expect(r.end).toBeCloseTo(Tp53.end, 0);
  });

  it('a click selects a gene; 0 fits the chromosome again', async () => {
    const { v } = await open();
    v.view.setRange(Section, Around.start, Around.end);
    await vi.advanceTimersByTimeAsync(3000);
    v.draw();
    clickAt(v, tp53Point(v));
    expect(v.selection.value.value.elements.keys.map((k) => k.key)).toEqual(['TP53']);
    expect(v.input.key({ key: '0' })).toBe(true);
    await vi.advanceTimersByTimeAsync(1000);
    const r = v.view.get(Section)!;
    expect([r.start, r.end]).toEqual([0, DefaultChromosome.length]);
  });

  it('a drag on the ideogram’s window pans the view by as much of the chromosome', async () => {
    const { v } = await open();
    const span = 4_000_000;
    v.view.setRange(Section, 20_000_000, 20_000_000 + span);
    await vi.advanceTimersByTimeAsync(3000);
    v.draw();
    const nav = v.navigatorModel('chr17')!;
    const w = nav.draw()!.windows[0];
    expect(w.p0).toBeCloseTo(20_000_000, -3);
    const mid = (w.dx0 + w.dx1) / 2;
    const dx = 0.1 * (w.hi - w.lo);
    expect(nav.down(mid, 12)).toBe(true);
    nav.move(mid + dx, 12);
    nav.up();
    const r = v.view.get(Section)!;
    expect(r.start).toBeCloseTo(20_000_000 + 0.1 * DefaultChromosome.length, -5);
    expect(r.end - r.start).toBeCloseTo(span, -2);
  });
});
