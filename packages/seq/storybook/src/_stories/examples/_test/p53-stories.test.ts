/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  bookmarkRows,
  bookmarkView,
  composeViewer as composeAddressing,
  describeRun,
  followLink,
} from '../../features/addressing/compose.js';
import {
  cancelDrag,
  composeEditing,
  dragHeight,
  type EditingSession,
  endDrag,
  historyRows,
  proteinSection,
  setDomainsColor,
  toggleTrack,
} from '../../features/editing/compose.js';
import {
  addTrack,
  AddedTrack,
  describeLinkRun,
  LinkStatus,
  RegionText,
  trackHeight,
  trackOrder,
} from '../p53/compose.js';

// The p53 example stories, headless on a fake clock, each composed as its story composes it: the editing story's host
// stack, viewer and SVG export; the bookmarks story's delayed region and its text's links. The selection example is the
// Features story's root with a hint: its drag on the gene is the Features twin's.

const tracks = P53Views.find((v) => v.id === 'tracks')!;
const bookmarks = FeatureViews.find((v) => v.id === 'bookmarks')!;
const DelayMs = 300;

let viewer: SeqViewer | undefined;
let session: EditingSession | undefined;

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
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

/** The viewer's SVG export, its asynchronous work let go on the fake clock. */
async function svgOf(v: SeqViewer, idPrefix: string): Promise<string> {
  const p = v.renderSvg({ idPrefix });
  await vi.advanceTimersByTimeAsync(300);
  return (await p).svg;
}

describe('p53 Edit and undo, headless', () => {
  async function make() {
    const s = composeEditing(tracks);
    session = s;
    s.viewer.resize(720);
    await s.viewer.load(tracks.build());
    await vi.advanceTimersByTimeAsync(400);
    return s;
  }

  it('add a track, hide one, set a theme: three steps, undone and redone, the SVG following', async () => {
    const s = await make();
    const svg0 = await svgOf(s.viewer, 'e-');
    expect(trackOrder(s)).toEqual(['letters', 'domains', 'conservation']);

    expect(await addTrack(s)).toBe('Added');
    expect(await addTrack(s)).toBe('Already added');
    expect(await toggleTrack(s, 'conservation')).toBe('Hidden');
    expect(await setDomainsColor(s, '#d62728')).toBe('Colored');
    expect(trackOrder(s)).toEqual(['letters', 'domains', AddedTrack, 'conservation']);
    expect(historyRows(s).done).toEqual(['Color a layer', 'Hide a track', 'Add a track']);
    await vi.advanceTimersByTimeAsync(400);
    const svg1 = await svgOf(s.viewer, 'e-');
    expect(svg1).not.toBe(svg0);
    expect(svg1).toContain('Domains (added)');
    expect(svg1.toLowerCase()).toContain('#d62728');

    for (let i = 0; i < 3; i++) await s.stack.undo();
    await vi.advanceTimersByTimeAsync(400);
    expect(trackOrder(s)).toEqual(['letters', 'domains', 'conservation']);
    expect(await svgOf(s.viewer, 'e-')).toBe(svg0);
    await s.stack.redo();
    expect(trackOrder(s)).toContain(AddedTrack);
  });

  it('a drag previews a height in the spec shown, Esc ends it in nothing, a release in one step', async () => {
    const s = await make();
    const heightOf = (shown: boolean) =>
      proteinSection(s, shown)?.children?.find((c) => c.id === 'domains')?.params?.height;
    const original = heightOf(true);
    expect(original).toBeUndefined();
    // the slider follows the spec: the laid-out height while nothing is set
    expect(trackHeight(s, 'domains')).toBeGreaterThan(0);
    expect(dragHeight(s, 'domains', 40)).toBe(true);
    expect(dragHeight(s, 'domains', 52)).toBe(true);
    expect(heightOf(true)).toBe(52);
    expect(trackHeight(s, 'domains')).toBe(52);
    expect(heightOf(false)).toBe(original);
    expect(cancelDrag(s)).toBe('Cancelled');
    expect(heightOf(true)).toBe(original);
    expect(historyRows(s).done).toEqual([]);

    dragHeight(s, 'domains', 52);
    await endDrag(s);
    expect(heightOf(false)).toBe(52);
    expect(heightOf(true)).toBe(52);
    expect(historyRows(s).done).toEqual(['Resize a track']);
    await s.stack.undo();
    expect(heightOf(true)).toBe(original);
  });
});

describe('p53 Bookmarks and links, headless', () => {
  async function make() {
    const text = regionSeqText();
    const v = composeAddressing(bookmarks, {
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
    expect(describeRun(run)).toBe('Waiting for its data; The address waits for its data');
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
    links.loaded();
    expect(links.follow(hrefOf('T2'))).toBe('Waiting for its data');
    expect(said).toEqual([]);
    await vi.advanceTimersByTimeAsync(4 * DelayMs + 64);
    expect(said).toEqual(['Focused']);
    // a link that is held says its outcome alone, and a later link ends an earlier wait
    expect(links.follow(hrefOf('the gene'))).toBe('Focused');
    links.dispose();
    expect(describeLinkRun(undefined)).toBe('Not a link');
  });

  it('a link followed before the view is loaded waits, then runs once it is', async () => {
    const text = regionSeqText();
    const v = composeAddressing(bookmarks, {
      delayMs: DelayMs,
      ports: { assets: () => ({ reader: { readText: async () => text } }) },
    });
    viewer = v;
    v.resize(720);
    const said: string[] = [];
    const links = new LinkStatus(v, (m) => said.push(m));
    expect(links.follow(hrefOf('the gene'))).toBe('Waiting for the view to load');
    // the viewer itself, asked before a spec is loaded, names nothing
    expect(describeLinkRun(followLink(v, hrefOf('the gene')))).toContain('Nothing to focus');
    await v.load(bookmarks.build());
    const settled = v.settle();
    await vi.advanceTimersByTimeAsync(4 * DelayMs + 400);
    await settled;
    links.loaded();
    await vi.advanceTimersByTimeAsync(64);
    expect(said).toEqual(['Focused']);
    expect(v.viewport && [v.viewport.start, v.viewport.end]).toEqual([2000, 3500]);
    links.dispose();
  });
});
