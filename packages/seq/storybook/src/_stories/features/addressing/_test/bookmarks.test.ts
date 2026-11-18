/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { bookmarkRows, bookmarkView, composeViewer, focusBookmark } from '../compose.js';

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

describe('Bookmarks, headless', () => {
  it('waits for a transcript’s bookmark, focuses it when it arrives, and bookmarks the view in one patch', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const text = regionSeqText();
    const v = composeViewer(view, {
      delayMs: DelayMs,
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    viewer = v;
    v.resize(720);
    await v.load(view.build());
    await vi.advanceTimersByTimeAsync(32);

    expect(bookmarkRows(v).map((r) => [r.address, r.status])).toEqual([
      ['/refs/gene', 'resolved'],
      ['/refs/t2', 'waiting'],
      ['/refs/wrap53', 'waiting'],
    ]);
    // the focus is asked before the transcript is held: it waits
    expect(focusBookmark(v, '/refs/t2')).toBe('Waiting for its data');
    const section = v.view.resolved.value!.sections[0].address;

    // the region answers: the transcripts are held, the bookmarks resolve and the waiting focus runs
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(bookmarkRows(v).map((r) => [r.address, r.status])).toEqual([
      ['/refs/gene', 'resolved'],
      ['/refs/t2', 'resolved'],
      ['/refs/wrap53', 'resolved'],
    ]);
    expect([v.view.get(section)!.start, v.view.get(section)!.end]).toEqual([2100, 2780]);

    // a bookmark of what is shown: one patch, listed at once
    const revision = v.spec.revision.value;
    expect(bookmarkView(v)).toBe('Bookmarked');
    expect(v.spec.revision.value).not.toBe(revision);
    expect(bookmarkRows(v).at(-1)).toMatchObject({
      address: '/refs/view-1',
      label: 'View 1',
      to: '/refs/t2',
      status: 'resolved',
    });
    expect(focusBookmark(v, '/refs/gene')).toBe('Focused');
    expect([v.view.get(section)!.start, v.view.get(section)!.end]).toEqual([2000, 3500]);
  });
});
