/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { recordDrawLists } from '@molstar/seqx-core/render';
import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile } from '@molstar/seqx-examples/p53/data';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
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
  setDomainsColor,
  toggleTrack,
} from '../../features/editing/compose.js';
import {
  cellPoint,
  composeViewer as composeSelection,
  describeSelection,
  dragFrom,
  selectionSvg,
} from '../../features/selection/compose.js';
import { addTrack, AddedTrack, RegionText, trackOrder } from '../p53/compose.js';

// The p53 example stories, headless on a fake clock, each composed as its story composes it: the editing story's host
// stack, viewer and SVG export; the bookmarks story's delayed region and its text's links; the selection story's drag on
// the gene, in the canvas's lists and the SVG's.

const tracks = P53Views.find((v) => v.id === 'tracks')!;
const onGene = P53Views.find((v) => v.id === 'protein-on-gene')!;
const bookmarks = FeatureViews.find((v) => v.id === 'bookmarks')!;
const DelayMs = 300;

let viewer: SeqViewer | undefined;
let session: EditingSession | undefined;
let gene = '';
beforeAll(async () => {
  gene = await nodeFixtureAssets().reader.readText(GeneFile, new AbortController().signal);
});

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

  it('a drag previews a height, Esc ends it in nothing, a release in one step', async () => {
    const s = await make();
    expect(dragHeight(s, 'domains', 40)).toBe(true);
    expect(dragHeight(s, 'domains', 52)).toBe(true);
    expect(cancelDrag(s)).toBe('Cancelled');
    expect(historyRows(s).done).toEqual([]);
    dragHeight(s, 'domains', 52);
    await endDrag(s);
    expect(historyRows(s).done).toEqual(['Resize a track']);
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
});

describe('p53 Selection, headless', () => {
  it('a drag on the gene selects nucleotides, projected onto the protein in the canvas and the SVG alike', async () => {
    const v = composeSelection(onGene, { ports: { assets: () => ({ reader: { readText: async () => gene } }) } });
    viewer = v;
    v.input.mac = false;
    v.resize(720);
    await v.load(onGene.build());
    await vi.advanceTimersByTimeAsync(400);
    v.view.setRange('/views/main/protein', 170, 180);
    await vi.advanceTimersByTimeAsync(400);
    v.draw();
    expect(describeSelection(v)).toBe('selection: — · cursor: —');

    const s = v.frame!.layout.sections[1];
    const y = s.tracks[0].box.y + 6;
    const a = Math.floor(s.axis!.visible([])[0]) + 3;
    dragFrom(v, [cellPoint(v, 1, 0, a)[0], y], [cellPoint(v, 1, 0, a + 8)[0], y]);
    expect(v.selection.value.value.shown.parts.map((p) => p.system)).toEqual(['/data/tp53']);
    expect(describeSelection(v)).not.toContain('selection: —');
    v.draw();

    const lists = recordDrawLists(v.frame!.lists, ['selection']);
    // the protein's section (0) carries the projection, lighter than the gene's own selection
    const rows = lists.split('\n').filter((l) => l.startsWith('rect '));
    const inSection = (n: number) =>
      rows.filter((l) => {
        const sec = v.frame!.layout.sections[n];
        const ry = Number(l.split(' ')[2]);
        return ry >= sec.y && ry < sec.y + sec.height;
      });
    expect(inSection(0).length).toBeGreaterThan(0);
    expect(inSection(1).length).toBeGreaterThan(0);
    const svg = await selectionSvg(v, 'p-');
    expect(recordDrawLists(svg.lists, ['selection'])).toBe(lists);
    expect(await svgOf(v, 'p-')).toContain('<svg');
  });
});
