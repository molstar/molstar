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
import { ColorByEntry, colorByBehavior } from '../../../support/host-adapter.js';
import { fakeClock } from '../../structure-link/_test/open.js';
import { headerLabel, headersOf, headersSpec, type HeadersComposition, menuLine, undoKey, widgetPlaces } from '../compose.js';

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
  const spec = headersSpec({
    ...options,
    registry: [...(v.registry?.({ delayMs: 0 }) ?? []), ...(options.registry ?? [])],
  });
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
    // the column holds one slot, so the two states are both marks, never a chip beside a mark
    expect(trackOf(v, 'domains').states.map((x) => x.mark)).toEqual([true, true]);
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
    const itemCount = () => headersOf(v).menuOf(trackOf(v, 'domains').target)!.groups[0].length;
    expect(itemCount()).toBe(7);
    await clickChip(v, trackOf(v, 'domains'), HeaderEntryNames.filtered);
    expect(chips(trackOf(v, 'domains'))).toEqual(['1 layer hidden']);
    // "Clear the filter" is gone from its menu
    expect(itemCount()).toBe(6);
    await undo(v);
    expect(chips(trackOf(v, 'domains'))).toEqual(['1 filtered', '1 layer hidden']);
    expect(itemCount()).toBe(7);
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

  it('the host’s Cmd/Ctrl+Z keys undo a chip’s click and Shift redoes it (F15)', async () => {
    const v = await open('headers');
    const domains = () => chips(trackOf(v, 'domains'));
    await clickChip(v, trackOf(v, 'domains'), HeaderEntryNames.filtered);
    expect(domains()).toEqual(['1 layer hidden']);
    expect(undoKey(v, { key: 'a', metaKey: true })).toBe(false);
    expect(undoKey(v, { key: 'z' })).toBe(false);
    await undoKey(v, { key: 'z', metaKey: true });
    await vi.advanceTimersByTimeAsync(20);
    await v.settle();
    expect(domains()).toEqual(['1 filtered', '1 layer hidden']);
    await undoKey(v, { key: 'Z', ctrlKey: true, shiftKey: true });
    await vi.advanceTimersByTimeAsync(20);
    await v.settle();
    expect(domains()).toEqual(['1 layer hidden']);
  });

  it('the DOM parts are a composition of the same behavior: the parts port is the story arg, nothing else differs', async () => {
    const params = (o: HeadersComposition) =>
      headersSpec(o).behaviors.find((b) => b.provider === HeadersProvider)!.params as HeadersParams;
    expect(params({}).parts).toBeUndefined();
    expect(typeof params({ parts: 'dom' }).parts).toBe('function');
    // the DOM they draw is checked in jsdom (track-headers-parts.test.ts)
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

  it('ST3’s "Color by this track" entry takes a slot (the box icon) on a track with a color scheme', async () => {
    const v = await open('headers', {
      registry: [ColorByEntry],
      behaviors: [colorByBehavior(() => 'colored')],
      labelWidth: 170,
    });
    expect(trackOf(v, 'letters').slots.map((x) => x.name)).toContain('fake-host:color-by');
    // the entry applies where a layer has schemes: not on a track of ranges without a color channel's choice
    const labels = headersOf(v)
      .menuOf(trackOf(v, 'letters').target)!
      .groups.flat()
      .map((i) => i.label);
    expect(labels).toContain('Color the structure by this track');
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

  it('a secondary click on the chip lists each hidden track in order, then "Show all"; an item shows that one', async () => {
    const v = await open('headers-hidden');
    const [x, y] = mid(sectionOf(v, 'hidden')!.states[0].box);
    expect(headersOf(v).context(x, y)).toBe(true);
    const shown = headersOf(v).menu.value!;
    expect(shown.menu.target).toEqual(sectionOf(v, 'hidden')!.target);
    const items = shown.menu.groups.flat();
    expect(items.map((i) => i.label)).toEqual(['Show Domains', 'Show Conservation', 'Show all']);
    await headersOf(v).runItem(items[1]);
    await vi.advanceTimersByTimeAsync(20);
    await v.settle();
    expect(visible(v)[0].map((a) => a.split('/').pop())).toEqual(['conservation']);
    expect(chips(sectionOf(v, 'hidden'))).toEqual(['1 hidden']);
  });
});

describe('short tracks', () => {
  it('their states are marks and they have no slots, the tall track’s header is whole', async () => {
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
    // the story names the menu's track by its label (F17)
    expect(headerLabel(v, target)).toBe('Forward strand');
    expect(headerLabel(v, { section: '/views/main/ring', track: '/views/main/ring/sites' })).toBe('Sites');
    // the host's button: the menu opens at a place and is the behavior's state
    expect(headersOf(v).openMenu(target, { x: 120, y: 60 })).toBe(true);
    expect(headersOf(v).menu.value?.menu.target).toEqual(target);
    headersOf(v).closeMenu();
    expect(headersOf(v).menu.value).toBeUndefined();
    // a secondary click on the ring's band opens the menu of its section or track, never the browser's
    const ring = v.frame!.layout.sections.find((l) => l.section.address === '/views/main/ring')!;
    const track = ring.tracks[0];
    const at = [track.box.x + track.box.width / 2, track.box.y + track.box.height / 2] as const;
    expect(headersOf(v).context(...at)).toBe(true);
    expect(headersOf(v).menu.value?.menu.target.section).toBe('/views/main/ring');
    headersOf(v).closeMenu();
    // the detail below it has its headers as any linear section's
    expect(geometry(v).filter((g) => g.target.section.endsWith('/detail')).length).toBe(4);
  });
});

describe('a host’s widget from the geometry', () => {
  it('one per shown track, left of its label box, centered on it; hiding a track takes its widget', async () => {
    const v = await open('headers');
    const places = widgetPlaces(geometry(v));
    expect(places.map((p) => p.target.track?.split('/').pop())).toEqual(['letters', 'domains', 'conservation']);
    for (const p of places) {
      const label = geometry(v).find((g) => g.target.track === p.target.track)!.label;
      // the badge is the label box's neighbor on its left, centered on it
      expect(p.x + p.size + 2).toBe(label.x);
      expect(p.y + p.size / 2).toBeCloseTo(label.y + label.height / 2, 9);
    }
    await v.execute('track.hide', { address: '/views/main/protein/domains' });
    await vi.advanceTimersByTimeAsync(20);
    await v.settle();
    expect(widgetPlaces(geometry(v)).map((p) => p.target.track?.split('/').pop())).toEqual(['letters', 'conservation']);
  });
});
