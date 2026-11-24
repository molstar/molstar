/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { locationsEqual } from '@molstar/seqx-core/query';
import { recordDrawLists } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { bookmarkRows, bookmarkView, composeViewer, focusBookmark, followLink } from '../compose.js';

// The Bookmarks story, headless on a fake clock, composed as the story composes it: the bookmarks on transcripts wait
// until their region's answer brings them; a focus asked of one meanwhile runs when it arrives; bookmarking the view is
// one patch.

const DelayMs = 300;
const view = FeatureViews.find((v) => v.id === 'bookmarks')!;

let viewer: SeqViewer | undefined;

// A frame tick asked for by an earlier file on timers since dropped would never run, and the frame lane would wait for
// it: a host function of this file's own asks for a fresh one (16 ms on the faked clock).
const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** The window the canvas draws (the first leading section's viewport), and whether it is placed. */
function drawn(v: SeqViewer) {
  const section = v.view.resolved.value!.sections[0].address;
  return { window: v.viewport && [v.viewport.start, v.viewport.end], pending: v.view.get(section)?.pending };
}

async function make(which = view) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const text = regionSeqText();
  const v = composeViewer(which, {
    delayMs: DelayMs,
    ports: {
      assets: () => ({ reader: { readText: async () => text, readBytes: async () => new TextEncoder().encode(text) } }),
    },
  });
  viewer = v;
  v.resize(720);
  await v.load(which.build());
  await vi.advanceTimersByTimeAsync(32);
  return v;
}

describe('Bookmarks, headless', () => {
  it('places the source-backed region once its length arrives, and focuses the drawn window', async () => {
    const v = await make();
    // the region's length comes with its source's coverage: the view is placed and drawn whole
    expect(drawn(v)).toEqual({ window: [0, 6000], pending: false });
    expect(focusBookmark(v, '/refs/gene')).toBe('Focused');
    await vi.advanceTimersByTimeAsync(64);
    expect(drawn(v)).toEqual({ window: [2000, 3500], pending: false });
    // a bookmark of the view stores the region drawn, never the ref focused: bookmarks never chain
    expect(bookmarkView(v)).toBe('Bookmarked');
    expect(bookmarkRows(v).at(-1)).toMatchObject({
      address: '/refs/view-1',
      to: '/data/tp53-region/(range :from 2001 :to 3500)',
      status: 'resolved',
    });
  });

  it('runs a waiting focus on an address that is no ref, once its node arrives', async () => {
    const v = await make(P53Views.find((x) => x.id === 'gene-by-region')!);
    const run = followLink(v, '!seq-focus=/data/region/(node :key T3)');
    expect(run?.outcomes.map((o) => o.message)).toEqual(['Waiting for T3']);
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(v.addresses.resolve('/data/region/(node :key T3)').status).toBe('resolved');
    await vi.advanceTimersByTimeAsync(64);
    // T3: the gene's third exon, [960, 1422) of the gene at 2000 in the region
    expect(drawn(v)).toEqual({ window: [2960, 3422], pending: false });
  });

  it('waits for a transcript’s bookmark, focuses it when it arrives, and bookmarks the view in one patch', async () => {
    const v = await make();
    expect(bookmarkRows(v).map((r) => [r.address, r.status])).toEqual([
      ['/refs/gene', 'resolved'],
      ['/refs/t2', 'waiting'],
      ['/refs/wrap53', 'waiting'],
    ]);
    // the focus is asked before the transcript is held: it waits
    expect(focusBookmark(v, '/refs/t2')).toBe('Waiting for T2');

    // the region answers: the transcripts are held, the bookmarks resolve and the waiting focus runs
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(bookmarkRows(v).map((r) => [r.address, r.status])).toEqual([
      ['/refs/gene', 'resolved'],
      ['/refs/t2', 'resolved'],
      ['/refs/wrap53', 'resolved'],
    ]);
    await vi.advanceTimersByTimeAsync(64);
    expect(drawn(v)).toEqual({ window: [2100, 2780], pending: false });

    // a bookmark of what is shown: one patch, listed at once
    const revision = v.spec.revision.value;
    expect(bookmarkView(v)).toBe('Bookmarked');
    expect(v.spec.revision.value).not.toBe(revision);
    expect(bookmarkRows(v).at(-1)).toMatchObject({
      address: '/refs/view-1',
      label: 'View 1',
      to: '/data/tp53-region/(range :from 2101 :to 2780)',
      status: 'resolved',
    });
  });

  it('two consumers of one ref to a dynamic node: the viewer and an export of its own resolve it alike', async () => {
    const v = await make();
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(v.addresses.select('/refs/t2').done).toBe(true);
    await vi.advanceTimersByTimeAsync(64);
    // the export: its own context and store, reading the region itself, resolving the selection's address itself
    const text = regionSeqText();
    const ctx = createSeqContext({
      registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view.registry?.({ delayMs: 0 }) ?? [])],
      ports: {
        assets: () => ({
          reader: { readText: async () => text, readBytes: async () => new TextEncoder().encode(text) },
        }),
      },
    });
    const exported = ctx.renderSvg(v.spec.root.value, {
      size: { width: v.width, height: v.height },
      ...(v.viewport ? { viewport: v.viewport } : {}),
      dpr: v.dpr,
      selection: v.addresses.selectionAddress()!,
    });
    await vi.advanceTimersByTimeAsync(1000);
    const theirs = await exported;
    ctx.retire();
    expect(theirs.settled).toBe(true);
    expect(theirs.selection && locationsEqual(theirs.selection, v.selection.value.value.shown)).toBe(true);
    const drawn = recordDrawLists(v.frame!.lists, ['selection']);
    expect(drawn).toContain('rect');
    expect(recordDrawLists(theirs.lists, ['selection'])).toBe(drawn);
  });
});
