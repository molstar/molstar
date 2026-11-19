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
import {
  bookmarkRows,
  bookmarkView,
  composeViewer,
  describeRun,
  followLink,
  LinkStatus,
  RegionText,
} from '../compose.js';

// The Bookmarks story's text, headless on a fake clock, composed as the story composes it: the delayed region and the
// text's links (a transcript's waits for its data, then moves the view; "Bookmark this view" adds a ref).

const bookmarks = FeatureViews.find((v) => v.id === 'bookmarks')!;
const DelayMs = 300;

let viewer: SeqViewer | undefined;

// a frame tick asked for by an earlier file on timers since dropped would never run: this file's own host function
const raf = DomHost.requestAnimationFrame;
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

describe('Bookmarks: the text’s links, headless', () => {
  async function make() {
    const text = regionSeqText();
    const v = composeViewer(bookmarks, {
      delayMs: DelayMs,
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    viewer = v;
    v.resize(720);
    await v.load(bookmarks.build());
    await vi.advanceTimersByTimeAsync(32);
    return v;
  }
  const hrefOf = (text: string) => {
    const r = RegionText.find((x) => typeof x !== 'string' && x.text === text);
    return typeof r === 'object' ? r.href : '';
  };

  it('a link on a transcript that arrives later waits, then moves the view; "bookmark this view" adds a ref', async () => {
    const v = await make();
    const run = followLink(v, hrefOf('T2'));
    expect(describeRun(run)).toBe('Waiting for its data');
    expect(bookmarkRows(v).find((r) => r.address === '/refs/t2')?.status).toBe('waiting');
    await vi.advanceTimersByTimeAsync(4 * DelayMs);
    expect(bookmarkRows(v).find((r) => r.address === '/refs/t2')?.status).toBe('resolved');
    await vi.advanceTimersByTimeAsync(64);
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual([2100, 2780]);

    // a link to what is held moves the view at once
    expect(describeRun(followLink(v, hrefOf('the gene')))).toContain('Focused');
    await vi.advanceTimersByTimeAsync(64);
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual([2000, 3500]);
    expect(await bookmarkView(v)).toBe('Bookmarked');
    expect(bookmarkRows(v).at(-1)).toMatchObject({ address: '/refs/view-1', status: 'resolved' });
  });

  it('the status follows a wait: it says it waits once, then Focused when the view has moved', async () => {
    const v = await make();
    const said: string[] = [];
    const links = new LinkStatus(v, (m) => said.push(m));
    expect(links.follow(hrefOf('T2'))).toBe('Waiting for its data');
    expect(said).toEqual([]);
    await vi.advanceTimersByTimeAsync(4 * DelayMs + 64);
    expect(said).toEqual(['Focused']);
    // a link that is held says its outcome alone, and a later link ends an earlier wait
    expect(links.follow(hrefOf('the gene'))).toBe('Focused');
    links.dispose();
    expect(describeRun(undefined)).toBe('Not a link');
  });

  it('a link followed before the view is loaded is kept by the viewer, and runs once it is', async () => {
    const text = regionSeqText();
    const v = composeViewer(bookmarks, {
      delayMs: DelayMs,
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    viewer = v;
    v.resize(720);
    const said: string[] = [];
    const links = new LinkStatus(v, (m) => said.push(m));
    // one message: the viewer waits for the spec, nothing here holds the link
    expect(links.follow(hrefOf('the gene'))).toBe('Waiting for its data');
    await v.load(bookmarks.build());
    const settled = v.settle();
    await vi.advanceTimersByTimeAsync(4 * DelayMs + 400);
    await settled;
    await vi.advanceTimersByTimeAsync(64);
    expect(said).toEqual(['Focused']);
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual([2000, 3500]);
    links.dispose();
  });
});
