/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesLay } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { Tp53 } from '@molstar/seqx-examples/p53/chromosome-data';
import { HoverCardProvider } from '@molstar/seqx-viewer/behaviors';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type FollowedSvg, followSvg } from '../../follow-svg.js';

// The Chromosome story, headless on a fake clock, composed as `example-root.ts` composes it: the viewer with the
// chromosome's delayed runtime, and its SVG twin exporting through the viewer's store each time the view rests (two
// consumers of one store). The twin is the export of what the canvas shows at the fit and at TP53; TP53's exons, a
// node its region brought, stay drawn on the canvas once the exports that read them have ended; a double-click on
// TP53 focuses its span; the hover card names the gene and its span.

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
  return { v, host, fresh };
}

/** A track's first layer as laid out on the canvas, by its label. */
const layerOf = (v: SeqViewer, label: string) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === label)?.layers[0].lay as RangesLay | undefined;

/** The centre of TP53's mark on the genes' track. */
function tp53Point(v: SeqViewer): [number, number] {
  const lay = layerOf(v, 'Genes')!;
  const k = Array.from({ length: lay.hn }, (_, j) => j).find((j) => lay.keys[lay.hitem[j]] === 'TP53')!;
  expect(k, 'TP53 hit box').toBeGreaterThanOrEqual(0);
  return [(lay.hx0[k] + lay.hx1[k]) / 2, lay.box.y + (lay.hlane[k] + 0.5) * lay.laneHeight];
}

const Section = '/views/main/chr17';
/** About TP53: 35 kb, its gene a bar. */
const Around = { start: Tp53.start - 8000, end: Tp53.end + 8000 };

describe('p53’s chromosome, the story headless', () => {
  it('the canvas and its SVG twin show the same view at the fit and at TP53; its exons stay once the exports end', async () => {
    const { v, host, fresh } = await open();
    const r = v.view.get(Section)!;
    expect([r.start, r.end]).toEqual([0, 83_257_441]);
    expect(layerOf(v, 'Genes')!.items).toBeGreaterThan(1100);
    expect(strip(host.innerHTML)).toBe(await fresh());

    v.view.setRange(Section, Around.start, Around.end);
    // the region answers, TP53's exons arrive as a node; the view rests and the twin exports through the same store
    await vi.advanceTimersByTimeAsync(3000);
    const exons = layerOf(v, 'Exons TP53')!;
    expect(exons.pending).toEqual([]);
    expect(exons.m).toBeGreaterThan(0);
    expect(strip(host.innerHTML)).toBe(await fresh());
    expect(strip(host.innerHTML)).toContain('TP53');
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

    // the card says an item's label and extent (a ranges row's strand it does not say: the viewer's description)
    v.input.move(...tp53Point(v));
    await vi.advanceTimersByTimeAsync(600);
    const card = v.behavior(HoverCardProvider)!.shown.value!.card;
    expect(card.title).toBe('TP53');
    expect(card.rows).toContainEqual({ label: 'Extent', value: 'chr17 7,668,402–7,687,550' });

    v.input.doubleClick(...tp53Point(v));
    await vi.advanceTimersByTimeAsync(1000);
    const r = v.view.get(Section)!;
    expect(r.start).toBeCloseTo(Tp53.start, 0);
    expect(r.end).toBeCloseTo(Tp53.end, 0);
  });
});
