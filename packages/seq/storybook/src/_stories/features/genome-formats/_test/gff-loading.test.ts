/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { GenomeFormatsViews } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { withViewRegistry } from '../../../examples/view-registry.js';

const raf = DomHost.requestAnimationFrame;
const viewers: SeqViewer[] = [];
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  for (const v of viewers.splice(0)) v.retire();
  vi.useRealTimers();
});

const file = 'fixtures/genomes/Homo_sapiens.GRCh38.116.chr.TP53.gff3.gz.pages.json';
const bytes = new Uint8Array(readFileSync(`${PackageFolder}${file}`));
const decoder = new TextDecoder();
const reader = { readText: async () => decoder.decode(bytes), readBytes: async () => bytes };
const view = GenomeFormatsViews.find((v) => v.id === 'gff-tp53')!;

function watched(entries: readonly SeqRegistryEntry[], opens: Promise<unknown>[]): SeqRegistryEntry[] {
  return entries.map((e) =>
    e.kind === 'format-runtime'
      ? {
          ...e,
          runtime: {
            ...e.runtime,
            open: (source, ctx) => {
              const p = e.runtime.open(source, ctx);
              opens.push(p);
              return p;
            },
          } satisfies FormatRuntime,
        }
      : e,
  );
}

const Section = '/views/main/chr17';
const heightOf = (v: SeqViewer) =>
  v.frame?.layout.sections[0]?.tracks.find((t) => t.track.label === 'Transcripts')?.box.height;

/** The viewer over TP53's GFF3 on a fake clock, its Transcripts track's height recorded each time it changes. */
async function loading(window?: readonly [number, number]) {
  const base = DefaultSeqViewerSpec();
  const opens: Promise<unknown>[] = [];
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: withViewRegistry(base.registry ?? [], watched(view.registry?.({ delayMs: 150 }) ?? [], opens)),
    ports: { assets: () => ({ reader }) },
  });
  viewers.push(v);
  v.resize(900);
  const heights: number[] = [];
  const step = async (ms: number) => {
    for (let k = 0; k < ms; k += 25) {
      await vi.advanceTimersByTimeAsync(25);
      const h = heightOf(v);
      if (h !== undefined && heights[heights.length - 1] !== h) heights.push(h);
    }
  };
  void v.load(view.build()).then(() => v.settle());
  await step(100);
  if (window) v.view.setRange(Section, window[0], window[1]);
  await step(100);
  for (let n = -1; n !== opens.length; ) {
    n = opens.length;
    await Promise.allSettled(opens);
    await step(200);
  }
  await step(6000);
  return { v, heights, step };
}

describe('TP53’s GFF3 while loading: the transcripts track’s height', () => {
  it('grows to the content of the reachable region at the first settle and keeps it across a pan (never shrunk, never grown)', async () => {
    // zoomed to a stretch of a few transcripts: the lanes in view alone would give a short track
    const { v, heights, step } = await loading([0, 1000]);
    const first = heightOf(v)!;
    for (const [a, b] of [
      [8000, 9000],
      [14_000, 15_000],
      [20_000, 21_000],
      [26_000, 27_000],
      [0, 1000],
    ]) {
      v.view.setRange(Section, a, b);
      await step(2000);
      expect(heightOf(v)).toBe(first);
    }
    // the loading lane, then the content: one change of height, and the whole region's (47 lanes, the cap)
    expect(heights.slice(-2)).toEqual([22, first]);
    expect(heights.length).toBe(2);
    v.view.setRange(Section, 0, 30_000);
    await step(2000);
    expect(heightOf(v)).toBe(first);
  });

  it('a track grows once, never shrinks: the tracks below move at most once', async () => {
    const { heights } = await loading();
    expect(heights).toEqual([22, 240]);
  });
});
