/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import type { Box } from '@molstar/seqx-core/render';
import { ExampleViews, type ExampleView } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { HeadersProvider } from '@molstar/seqx-viewer/behaviors/headers';
import {
  HeaderEntryNames,
  type HeadersParams,
  type HeaderGeometry,
  type HeaderTarget,
} from '@molstar/seqx-viewer/header-entries';
import { type UiHistory, uiHistory } from '@molstar/seqx-viewer/ui-history';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { fakeClock } from '../../structure-link/_test/open.js';
import { headersOf, headersSpec, type HeadersComposition, menuLine, widgetPlaces } from '../compose.js';

// The track headers' stories, headless on a fake clock, composed as the stories compose them (compose.ts): the p53
// views with a filter, a hidden layer and hidden tracks, the plasmid's ring. The witnesses are facts of the view (what
// is hidden, what is filtered) and of the geometry the viewer publishes: a chip's count against the spec, boxes that
// never overlap, a click and its undo against the resolved view. The look is the stories'.

const view = (id: string): ExampleView => ExampleViews.find((v) => v.id === id)!;

let viewer: SeqViewer | undefined;
let restore: (() => void) | undefined;
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  restore?.();
  restore = undefined;
});

/** A settled headless viewer of a story's composition, on the fake clock. */
async function open(id: string, options: HeadersComposition = {}): Promise<SeqViewer> {
  const v = view(id);
  const spec = headersSpec({ ...options, registry: v.registry?.({ delayMs: 0 }) });
  const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
  const opened = createSeqViewer({
    ...spec,
    ui: undefined,
    history: uiHistory(),
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => read(uri).toString('utf8'),
          readBytes: async (uri) => new Uint8Array(read(uri)),
        },
      }),
    },
  });
  viewer = opened;
  opened.resize(1440);
  await opened.load(v.build());
  await opened.settle();
  restore = fakeClock();
  return opened;
}

const geometry = (v: SeqViewer): readonly HeaderGeometry[] => v.headers!.value;
const trackOf = (v: SeqViewer, id: string) => geometry(v).find((g) => g.target.track?.endsWith(`/${id}`))!;
const sectionOf = (v: SeqViewer, id: string) =>
  geometry(v).find((g) => !g.target.track && g.target.section.endsWith(`/${id}`));
const chips = (g: HeaderGeometry | undefined) => g?.states.map((s) => s.text) ?? [];
const mid = (b: Box) => [b.x + b.width / 2, b.y + b.height / 2] as const;
const overlap = (a: Box, b: Box) =>
  a.x < b.x + b.width && b.x < a.x + a.width && a.y < b.y + b.height && b.y < a.y + a.height;
const visible = (v: SeqViewer) => v.view.resolved.value!.sections.map((s) => s.tracks.map((t) => t.address));

/** The click of a story's pointer on a chip, and the patch it makes landed. */
async function clickChip(v: SeqViewer, g: HeaderGeometry, name: string) {
  const [x, y] = mid(g.states.find((s) => s.name === name)!.box);
  v.input.move(x, y);
  await vi.advanceTimersByTimeAsync(20);
  await v.input.click(x, y);
  await vi.advanceTimersByTimeAsync(20);
  await v.settle();
}
const undo = async (v: SeqViewer) => {
  await (v.spec.history as UiHistory).stack.undo();
  await vi.advanceTimersByTimeAsync(20);
  await v.settle();
};

describe('the built-ins on p53', () => {
  it('say what the view holds: a filter, a layer and two tracks hidden', async () => {
    const v = await open('headers');
    expect(chips(sectionOf(v, 'protein'))).toEqual(['2 hidden']);
    expect(chips(trackOf(v, 'domains'))).toEqual(['1 filtered', '1 layer hidden']);
    expect(chips(trackOf(v, 'letters'))).toEqual([]);
    // the hidden tracks are not in the resolved view, the count says how many are left out
    expect(visible(v)[0].map((a) => a.split('/').pop())).toEqual(['letters', 'domains', 'conservation']);
    // every shown track has its slots, hide and the menu's "⋯", and no box of a header overlaps another of it
    for (const t of ['letters', 'domains', 'conservation']) {
      const g = trackOf(v, t);
      expect(g.slots.map((s) => s.name)).toEqual([HeaderEntryNames.hide, 'seq:menu']);
      const parts = [...g.slots, ...g.states].map((p) => p.box);
      for (const [i, a] of parts.entries()) {
        for (const b of parts.slice(i + 1)) expect(overlap(a, b)).toBe(false);
      }
    }
    expect(menuLine(v, trackOf(v, 'letters').target)).toBe(
      'Hide track, Color, Move down | Reverse, Collapse or expand hidden stretches, Show all',
    );
  });

  it('a chip’s click undoes its state in one step, and one undo brings it back', async () => {
    const v = await open('headers');
    await clickChip(v, trackOf(v, 'domains'), HeaderEntryNames.filtered);
    expect(chips(trackOf(v, 'domains'))).toEqual(['1 layer hidden']);
    await undo(v);
    expect(chips(trackOf(v, 'domains'))).toEqual(['1 filtered', '1 layer hidden']);
    await clickChip(v, sectionOf(v, 'protein')!, HeaderEntryNames.hiddenTracks);
    expect(visible(v)[0].map((a) => a.split('/').pop())).toEqual([
      'letters',
      'domains',
      'conservation',
      'secondary',
      'variants',
    ]);
    expect(chips(sectionOf(v, 'protein'))).toEqual([]);
    await undo(v);
    expect(chips(sectionOf(v, 'protein'))).toEqual(['2 hidden']);
  });

  it('as DOM parts: the same items from the same behavior, the canvas composed with a parts port', async () => {
    const drawn = await open('headers');
    const items = (v: SeqViewer, track: string) =>
      headersOf(v)
        .partsState()
        .headers.find((h) => h.geometry.target.track?.endsWith(`/${track}`))!
        .items.map((i) => `${i.name} ${i.label}`);
    const canvas = items(drawn, 'domains');
    expect(canvas).toContain(`${HeaderEntryNames.hide} Hide track`);
    viewer?.retire();
    restore?.();
    const dom = await open('headers', { parts: 'dom' });
    expect(items(dom, 'domains')).toEqual(canvas);
    const params = (o: HeadersComposition) =>
      headersSpec(o).behaviors.find((b) => b.provider === HeadersProvider)!.params as HeadersParams;
    expect(params({}).parts).toBeUndefined();
    expect(typeof params({ parts: 'dom' }).parts).toBe('function');
  });
});

describe('a host’s entry and a built-in left out', () => {
  it('the host’s slot on every track runs its command; Color is in no slot and no menu', async () => {
    const lines: string[] = [];
    const v = await open('headers', {
      host: true,
      labelWidth: 170,
      without: [HeaderEntryNames.color],
      say: (l) => lines.push(l),
    });
    for (const t of ['letters', 'domains', 'conservation']) {
      const g = trackOf(v, t);
      expect(g.slots.map((s) => s.name)).toEqual([HeaderEntryNames.hide, 'story-host:copy-address', 'seq:menu']);
      const labels = headersOf(v)
        .menuOf(g.target)!
        .groups.flat()
        .map((i) => i.label);
      expect(labels).toContain('Copy the track’s address');
      expect(labels).not.toContain('Color');
    }
    const target: HeaderTarget = trackOf(v, 'domains').target;
    await headersOf(v).run('story-host:copy-address', target);
    expect(lines).toEqual(['Copied /views/main/protein/domains']);
  });

  it('without the host’s composition nothing of it is there', async () => {
    const v = await open('headers');
    expect(trackOf(v, 'letters').slots.map((s) => s.name)).toEqual([HeaderEntryNames.hide, 'seq:menu']);
    expect(menuLine(v, trackOf(v, 'letters').target)).toContain('Color');
  });
});

describe('the hidden-tracks chip', () => {
  it('a section whose tracks are all hidden keeps its header; the chip shows both, one undo hides them', async () => {
    const v = await open('headers-hidden');
    const g = sectionOf(v, 'hidden')!;
    expect(g.shown).toBe(true);
    expect(chips(g)).toEqual(['2 hidden']);
    expect(visible(v)[0]).toEqual([]);
    await clickChip(v, g, HeaderEntryNames.hiddenTracks);
    expect(visible(v)[0].map((a) => a.split('/').pop())).toEqual(['domains', 'conservation']);
    await undo(v);
    expect(visible(v)[0]).toEqual([]);
    expect(chips(sectionOf(v, 'hidden'))).toEqual(['2 hidden']);
  });

  it('its menu lists each hidden track by its label, then "Show all"', async () => {
    const v = await open('headers-hidden');
    const menu = headersOf(v).menuOf(sectionOf(v, 'hidden')!.target)!;
    const labels = menu.groups.flat().flatMap((i) => [i.label, ...(i.items?.map((s) => s.label) ?? [])]);
    expect(labels).toEqual(expect.arrayContaining(['Show all', 'Show Domains', 'Show Conservation']));
  });
});

describe('short tracks', () => {
  it('their states are marks and their slots wait for the pointer; the tall track’s header is whole', async () => {
    const v = await open('headers-short');
    const short = trackOf(v, 'domains');
    expect(short.line.height).toBe(14);
    expect(short.states.map((s) => s.mark)).toEqual([true, true]);
    expect(chips(short)).toEqual(['1 filtered', '1 layer hidden']);
    expect(short.slots).toEqual([]);
    expect(trackOf(v, 'conservation').slots).toEqual([]);
    const tall = trackOf(v, 'tall');
    expect(tall.slots.map((s) => s.name)).toEqual([HeaderEntryNames.hide, 'seq:menu']);
    // the marks stay in reach of the menu: what they say is an item there, each undoing its state
    const labels = headersOf(v)
      .menuOf(short.target)!
      .groups.flat()
      .map((i) => i.label);
    expect(labels).toEqual(expect.arrayContaining(['Hide track', 'Clear the filter', 'Show every layer']));
  });
});

describe('a ring through its menu', () => {
  it('the ring has no header line; its track’s menu is the linear one without Reverse or Collapse', async () => {
    const v = await open('ring');
    expect(geometry(v).some((g) => g.target.section.endsWith('/ring'))).toBe(false);
    const target = { section: '/views/main/ring', track: '/views/main/ring/plus' };
    expect(menuLine(v, target)).toBe('Hide track, Color, Move down');
    // the host's button: the menu opens at a place and is the behavior's state
    expect(headersOf(v).openMenu(target, { x: 120, y: 60 })).toBe(true);
    expect(headersOf(v).menu.value?.menu.target).toEqual(target);
    headersOf(v).closeMenu();
    expect(headersOf(v).menu.value).toBeUndefined();
    // the detail below it has its headers as any linear section's
    expect(geometry(v).filter((g) => g.target.section.endsWith('/detail')).length).toBe(4);
  });
});

describe('a host’s widget from the geometry', () => {
  it('one per shown track, left of its line, never over it; hiding a track takes its widget', async () => {
    const v = await open('headers');
    const places = widgetPlaces(geometry(v));
    expect(places.map((p) => p.target.track?.split('/').pop())).toEqual(['letters', 'domains', 'conservation']);
    for (const p of places) {
      const line = geometry(v).find((g) => g.target.track === p.target.track)!.line;
      expect(p.x + p.size).toBeLessThanOrEqual(line.x);
      expect(p.y).toBeGreaterThanOrEqual(line.y);
      expect(p.y + p.size).toBeLessThanOrEqual(line.y + line.height);
    }
    await v.execute('track.hide', { address: '/views/main/protein/domains' });
    await vi.advanceTimersByTimeAsync(20);
    await v.settle();
    expect(widgetPlaces(geometry(v)).map((p) => p.target.track?.split('/').pop())).toEqual(['letters', 'conservation']);
  });
});
