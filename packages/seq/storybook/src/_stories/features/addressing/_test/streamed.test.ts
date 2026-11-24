/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import { regionData, regionSeqText } from '@molstar/seqx-examples/p53/region';
import { openSeq } from '@molstar/seqx-core/data';
import { viewerItem } from '@molstar/seqx-viewer/config';
import { HistoryParticipants, hostHistory } from '@molstar/seqx-viewer/ui-history';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { HistoryModel } from '@molstar/ui/history/history-model';
import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
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
let stack: HistoryModel<Readonly<Record<string, unknown>>> | undefined;

// the viewer's node seeking and the source's key index are lazy chunks, loaded on first use: loaded here, off the fake
// clock (a chunk's load is no timer the clock can run)
beforeAll(async () => {
  await import('@molstar/seqx-viewer/models/node-seek');
  await openSeq(regionData('seq')).findNames?.('T1', {}, new AbortController().signal);
});

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
  stack?.dispose();
  stack = undefined;
  vi.useRealTimers();
});

/** The story's viewer on a small cache; with `history`, a host's undo stack it runs its edits through. */
async function make(history = false): Promise<SeqViewer> {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const text = regionSeqText();
  const participants = new HistoryParticipants();
  if (history) stack = new HistoryModel<Readonly<Record<string, unknown>>>(participants);
  const v = composeViewer(view, {
    delayMs: DelayMs,
    ...(stack ? { history: hostHistory(stack, participants) } : {}),
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

describe('undoing a bookmark on a streamed transcript', () => {
  const bookmarked = (v: SeqViewer) => bookmarkRows(v).some((r) => r.address === '/refs/iso2');

  it('undoes and redoes a bookmark whose transcript was evicted after it was added (EC-162, EC-175)', async () => {
    const v = await make(true);
    expect(await v.addresses.addRef({ id: 'iso2', to: '/data/region/(node :key T2)', label: 'T2' })).toBe(true);
    expect(statusOf(v, '/refs/iso2')).toBe('resolved');
    await evictGene(v);
    expect(statusOf(v, '/refs/iso2')).toBe('waiting');
    // undo and redo move the spec's root alone: nothing is read, nothing is dropped from the store
    const bytes = v.store.store.bytes;
    await stack!.undo();
    expect(bookmarked(v)).toBe(false);
    await stack!.redo();
    expect(statusOf(v, '/refs/iso2')).toBe('waiting');
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(v.store.store.bytes).toBe(bytes);
    // and the bookmark redone still brings its transcript when jumped to
    expect(focusBookmark(v, '/refs/iso2')).toBe('Waiting for T2');
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(statusOf(v, '/refs/iso2')).toBe('resolved');
  });

  it('undoes and redoes a bookmark added while its transcript was not held', async () => {
    const v = await make(true);
    await evictGene(v);
    expect(await v.addresses.addRef({ id: 'iso2', to: '/data/region/(node :key T2)', label: 'T2' })).toBe(true);
    expect(stack!.steps.value.map((s) => s.label)).toEqual(['Add a bookmark']);
    expect(statusOf(v, '/refs/iso2')).toBe('waiting');
    await stack!.undo();
    expect(bookmarked(v)).toBe(false);
    expect(stack!.redoable.value.map((s) => s.label)).toEqual(['Add a bookmark']);
    await stack!.redo();
    expect(statusOf(v, '/refs/iso2')).toBe('waiting');
    // its transcript's region read again (the view back on the gene): the bookmark resolves
    const section = v.view.resolved.value!.sections[0].address;
    v.view.setViewport(section, { start: 2000, end: 3200 });
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    expect(statusOf(v, '/refs/iso2')).toBe('resolved');
  });
});

/** The regions a waiting jump reads now, by their starts (the seek's own wants). */
const seeking = (v: SeqViewer) => v.addresses.seeking().map((w) => w.start);

describe('a jump on a node nested in a streamed one', () => {
  const nested = (key: string) => `/data/region/(node :key T2-domains)/(node :key ${key})`;

  it('waits for the outer node, then the inner one, and runs when both have arrived', async () => {
    const v = await make();
    await evictGene(v);
    const run = v.runLink(`!seq-focus=${nested('domains')}`);
    expect(run?.outcomes.map((o) => o.message)).toEqual(['Waiting for T2-domains, domains']);
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    // the outer node's region read, its nested source opened, the domains held: the jump ran (on the protein's
    // residues, which this view does not show) and let go of its reads
    expect(v.addresses.resolve(nested('domains')).status).toBe('resolved');
    expect(seeking(v)).toEqual([]);
  });

  it('lets go of its reads when the inner node cannot be asked for and does not arrive', async () => {
    const v = await make();
    await evictGene(v);
    expect(v.runLink(`!seq-focus=${nested('nothing')}`)?.outcomes.map((o) => o.message)).toEqual([
      'Waiting for T2-domains, nothing',
    ]);
    await vi.advanceTimersByTimeAsync(DelayMs + 100);
    expect(seeking(v)).toHaveLength(1);
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    // the outer node arrived; its nested source finds no names, so nothing more can be asked: its region is let go
    expect(v.store.store.getNode('/data/region', 'T2-domains')).toBeDefined();
    expect(v.addresses.resolve(nested('nothing')).status).toBe('waiting');
    expect(seeking(v)).toEqual([]);
  });
});

describe('a waiting jump replaced by another', () => {
  it('lets go of the first jump’s reads', async () => {
    const v = await make();
    await evictGene(v);
    expect(focusBookmark(v, '/refs/t2')).toBe('Waiting for T2');
    await vi.advanceTimersByTimeAsync(DelayMs + 100);
    expect(seeking(v)).toHaveLength(1);
    // another waiting jump (on a key the source then says it has none of): T2's region is no longer read
    expect(v.runLink('!seq-focus=/data/region/(node :key T9)')?.outcomes.map((o) => o.message)).toEqual([
      'Waiting for T9',
    ]);
    expect(seeking(v)).toEqual([]);
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(seeking(v)).toEqual([]);
  });
});

describe('a bookmark on a transcript still held after the view moves on', () => {
  it('stays resolved: its region is cached (not wanted, not evicted) though its track has gone', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const text = regionSeqText();
    // the story's own cache (64 MB): nothing is evicted
    const v = composeViewer(view, {
      delayMs: DelayMs,
      ports: {
        assets: () => ({
          reader: { readText: async () => text, readBytes: async () => new TextEncoder().encode(text) },
        }),
      },
    });
    viewer = v;
    v.resize(720);
    await v.load(view.build());
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    expect(focusBookmark(v, '/refs/t2')).toBe('Focused');
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    expect(focusBookmark(v, '/refs/wrap53')).toBe('Focused');
    await vi.advanceTimersByTimeAsync(10 * DelayMs);
    expect(v.store.store.getNode('/data/region', 'T2')).toBeDefined();
    expect(statusOf(v, '/refs/t2')).toBe('resolved');
  });
});
