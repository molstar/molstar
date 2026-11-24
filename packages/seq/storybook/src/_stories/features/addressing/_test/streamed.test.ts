/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { bookmarkRows, composeViewer, focusBookmark } from '../compose.js';

// Patches on streamed things, in the Bookmarks story's viewer (headless, on a fake clock): the cache holds about one
// stretch of the region, so moving past the gene evicts its transcripts. A bookmark on one then waits (never dangles);
// a jump on it says what it waits for and asks the source where the transcript is, reads that region and runs when the
// transcript arrives; a key the source does not have is said once, and the jump is dropped. Nothing is read for a
// bookmark that waits until something jumps to it.

const DelayMs = 300;
const view = FeatureViews.find((v) => v.id === 'bookmarks')!;
/** Bytes the store keeps: the region's file opens (16 kB), and the whole region's tiles do not fit. */
const CacheBytes = 17_000;

let viewer: SeqViewer | undefined;

const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
    return 0;
  };
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function make(): Promise<SeqViewer> {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const text = regionSeqText();
  const v = composeViewer(view, {
    delayMs: DelayMs,
    config: [
      [viewerItem('store.cacheBytes'), CacheBytes],
      [viewerItem('store.tileBins'), 256],
    ],
    ports: {
      assets: () => ({ reader: { readText: async () => text, readBytes: async () => new TextEncoder().encode(text) } }),
    },
  });
  viewer = v;
  v.resize(720);
  await v.load(view.build());
  await vi.advanceTimersByTimeAsync(10 * DelayMs);
  return v;
}

const statusOf = (v: SeqViewer, address: string) => bookmarkRows(v).find((r) => r.address === address)?.status;

/** Moves the view past the gene: the gene's tiles leave the cache, and T2 with them. */
async function evictGene(v: SeqViewer) {
  const section = v.view.resolved.value!.sections[0].address;
  v.view.setViewport(section, { start: 5700, end: 6000 });
  await vi.advanceTimersByTimeAsync(10 * DelayMs);
  expect(v.store.store.getNode('/data/region', 'T2')).toBeUndefined();
}

describe('a bookmark on an evicted transcript', () => {
  it('waits, and nothing is read for it until a jump asks (EC-176)', async () => {
    const v = await make();
    expect(statusOf(v, '/refs/t2')).toBe('resolved');
    await evictGene(v);
    expect(statusOf(v, '/refs/t2')).toBe('waiting');
    // a waiting ref holds no data: the store reads nothing more while nothing jumps to it
    const bytes = v.store.store.bytes;
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    expect(v.store.store.bytes).toBe(bytes);
    expect(statusOf(v, '/refs/t2')).toBe('waiting');
  });

  it('a jump says what it waits for, asks the source where it is, and runs when it arrives (EC-176)', async () => {
    const v = await make();
    await evictGene(v);
    expect(focusBookmark(v, '/refs/t2')).toBe('Waiting for T2');
    // the source finds T2 (after its delay), its region is read (after the read's), and the waiting focus runs
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(statusOf(v, '/refs/t2')).toBe('resolved');
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual([2100, 2780]);
    expect(v.addresses.issues.value).toEqual([]);
  });

  it('a link on the node itself, not a bookmark, waits and arrives the same way', async () => {
    const v = await make();
    await evictGene(v);
    const run = v.runLink('!seq-focus=/data/region/(node :key T2)');
    expect(run?.outcomes.map((o) => o.message)).toEqual(['Waiting for T2']);
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual([2100, 2780]);
  });
});

describe('a bookmark on a key the source does not have', () => {
  it('is added as one patch; a jump on it is told once that the source has none, and is dropped', async () => {
    const v = await make();
    expect(await v.addresses.addRef({ id: 't9', to: '/data/region/(node :key T9)', label: 'T9' })).toBe(true);
    expect(statusOf(v, '/refs/t9')).toBe('waiting');
    const before = v.viewport && [v.viewport.start, v.viewport.end];
    expect(focusBookmark(v, '/refs/t9')).toBe('Waiting for T9');
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    // said once, the bookmark then names nothing it can wait for, and the view stays where it was
    expect(v.addresses.issues.value.map((i) => [i.code, i.message])).toEqual([
      ['seq.address.node.unknown', 'Its source has no T9'],
    ]);
    expect(statusOf(v, '/refs/t9')).toBe('missing');
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual(before);
    // asked again, the answer is known: no wait, nothing said twice
    expect(focusBookmark(v, '/refs/t9')).toBe('Its source has no T9');
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(v.addresses.issues.value).toHaveLength(1);
  });
});
