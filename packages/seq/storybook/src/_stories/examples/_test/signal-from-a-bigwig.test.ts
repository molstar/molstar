/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { ConservationFile } from '@molstar/seqx-examples/p53/conservation';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { type ValuesLay } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { GeneOffset, RegionLength } from '@molstar/seqx-examples/p53/region';
import { GeneLength } from '@molstar/seqx-examples/p53/data';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';

// The Signal from a bigWig story, headless on a fake clock, composed as `example-root.ts` composes it: the viewer with
// the `bigwig` runtime delayed as the story's twin delays it. The whole region is drawn from the coarse levels first,
// and a pan to a stretch of the gene is read by region and level again: a bin a few pixels wide while the stretch is
// wide, a value a base at the closest zoom. Read by this test, the signal is what M5 measures a bigWig pan on.

const DelayMs = 200;
const view = P53Views.find((v) => v.id === 'signal-from-a-bigwig')!;
const Section = '/views/main/region';
const raf = DomHost.requestAnimationFrame;
let viewer: SeqViewer | undefined;

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const file = new Uint8Array(readFileSync(`${PackageFolder}${ConservationFile}`));
const bytes = {
  readText: async () => '',
  readBytes: async () => file,
};

async function open() {
  const base = DefaultSeqViewerSpec();
  const v = (viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), ...(view.registry?.({ delayMs: DelayMs }) ?? [])],
    // the file read once, up front: the fake clock does not drive a disk read, so an answer waits only on the twin's delay
    ports: { assets: () => ({ reader: bytes }) },
  }));
  v.resize(720);
  void v.load(view.build()).then(() => v.settle());
  await vi.advanceTimersByTimeAsync(5000);
  return v;
}

/** The signal track's values as laid out on the canvas. */
const signal = (v: SeqViewer) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Conservation')?.layers[0].lay as
    | ValuesLay
    | undefined;

describe('p53’s conservation from a bigWig, the story headless', () => {
  it('draws the whole region in bins, then a value a base on a stretch of the gene, each read after its delay', async () => {
    const v = await open();
    const r = v.view.get(Section)!;
    expect([r.start, r.end]).toEqual([0, RegionLength]);
    const fit = signal(v)!;
    expect(fit.pending).toEqual([]);
    expect(fit.failed).toEqual([]);
    // binned: about a bar a pixel, far fewer than the region has positions (the layout's object is reused: read now)
    const bins = fit.n;
    expect(bins).toBeGreaterThan(100);
    expect(bins).toBeLessThan(RegionLength / 4);

    // a pan to 60 bases of the gene: the closest zoom, one bar a position
    const start = GeneOffset + Math.floor(GeneLength / 2);
    v.view.setRange(Section, start, start + 60);
    await vi.advanceTimersByTimeAsync(5000);
    v.draw();
    const near = signal(v)!;
    expect(near.pending).toEqual([]);
    expect(near.failed).toEqual([]);
    expect(near.n).toBe(60);
  });

  it('a pan back to the fit is served from what it read: no answer is waited on', async () => {
    const v = await open();
    v.view.setRange(Section, GeneOffset, GeneOffset + 400);
    await vi.advanceTimersByTimeAsync(5000);
    v.view.setRange(Section, 0, RegionLength);
    await vi.advanceTimersByTimeAsync(2 * 16);
    v.draw();
    expect(signal(v)!.pending).toEqual([]);
  });
});
