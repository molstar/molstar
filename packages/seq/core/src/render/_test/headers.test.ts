/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headers' line (_spec/README.md "Render", Track headers): the slots in fixed columns, never over each other, a
// state or the "⋯", their x the column's whatever the label; the states beside the label, on a line of their own, or
// marks; a section's line under `section.header` counting each hidden track once; the parts drawn alike in the draw
// lists and their SVG, hit as laid out, the slots only while their track's actions show and never in an export.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  actionColumns,
  bundledTextMetrics,
  drawListsToSvg,
  type Emphasis,
  type HeaderIcon,
  type HeaderParts,
  type HeaderPartsOf,
  hitView,
  recordDrawLists,
  renderConfig,
  type SectionLayout,
} from '../index.js';
import {
  builtInHeaderStates,
  drawsSectionHeader,
  HeaderMenuName,
  sectionHeaderLine,
  trackHeaderLine,
} from '../headers.js';
import { renderSvg } from '../render-svg.js';
import { P53, seeded } from './fixtures.js';
import { settle } from './helpers.js';

const Text = bundledTextMetrics();
/** A cross: two strokes in the 24-unit box. */
const Cross: HeaderIcon = { commands: [0, 1, 0, 1], points: [4, 4, 20, 20, 4, 20, 20, 4] };

/** A host's entries: two slots and the menu on every track, the built-in states, a slot on a section with a state. */
const Host: HeaderPartsOf = (f) => {
  const own = builtInHeaderStates(f);
  if (!f.track) return own && { ...own, menu: true };
  return {
    slots: [
      { name: 'seq:hide', label: 'Hide track', icon: Cross, enabled: true },
      { name: 'seq:color', label: 'Color', enabled: false },
    ],
    ...(own?.states ? { states: own.states } : {}),
    menu: true,
  };
};

const Rows = {
  start: [0, 63, 93, 318, 323],
  end: [61, 94, 292, 360, 355],
  id: ['a', 'b', 'c', 'd', 'e'],
  category: ['activation', 'linker', 'binding', 'linker', 'oligomer'],
};
const Filter = '/(not :of (item :category linker))';

/** p53: a letters track, a filtered ranges track (two lines high), a short filtered one, and the tracks `hidden` hides. */
function spec(hidden: readonly string[] = []): Node {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
  const d = b.data.ranges('d', { target: p, rows: Rows });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    const t = (id: string, params: Record<string, unknown>, layer: Record<string, unknown>) =>
      s.track(id, { ...params, ...(hidden.includes(id) ? { hidden: true } : {}) } as never).layer('l', layer as never);
    t('seq', { label: 'Sequence' }, { type: 'letters', data: p });
    t('dom', { label: 'Domains' }, { type: 'ranges', data: `${d}${Filter}` });
    t('thin', { label: 'Thin', height: 10 }, { type: 'ranges', data: `${d}${Filter}` });
    t('more', { label: 'More domains' }, { type: 'ranges', data: d });
  });
  return b.build();
}

const Size = { width: 720, height: 260 };

async function frame(hidden: readonly string[] = [], headers?: HeaderPartsOf, emphasis?: Emphasis) {
  // a column of 140 px: two slots, the "⋯" and the chips beside the labels (an 80 px one holds one slot, its chips marks)
  const s = await settle(spec(hidden), [], { labelColumn: { min: 140 } });
  const input = { size: Size, settled: true, ...(headers ? { headers } : {}), ...(emphasis ? { emphasis } : {}) };
  return { s, f: await s.frame(input), input };
}

/** Two boxes overlap (their insides). */
const overlap = (a: Box, b: Box) =>
  a.x < b.x + b.width && b.x < a.x + a.width && a.y < b.y + b.height && b.y < a.y + a.height;
interface Box {
  x: number;
  y: number;
  width: number;
  height: number;
}

describe('a header line’s geometry', () => {
  it('slots never overlap each other, a state or the “⋯”; a slot’s x is the column’s, whatever the label and its states', () => {
    const r = seeded(29);
    const words = ['1 filtered', '12 filtered', '1 layer hidden', '3 layers hidden', 'loading', 'x'];
    for (let n = 0; n < 400; n++) {
      const config = renderConfig({
        labelWidth: 80 + Math.floor(r() * 120),
        header: { slots: Math.floor(r() * 5), slotPx: 16 + Math.floor(r() * 8) },
      });
      const slots = Array.from({ length: Math.floor(r() * 5) }, (_, i) => ({
        name: `e:${i}`,
        label: `Entry ${i}`,
        enabled: r() < 0.7,
      }));
      const states = Array.from({ length: Math.floor(r() * 4) }, (_, i) => ({
        name: `s:${i}`,
        text: words[Math.floor(r() * words.length)],
      }));
      const parts: HeaderParts = { slots, states, menu: r() < 0.8 };
      const height = [8, 14, 20, 24, 36, 60][Math.floor(r() * 6)];
      const well = { x: 200, y: 40, width: 400, height };
      const { header } = trackHeaderLine(parts, '/v/s/t', 's', well, config, Text, false);
      const cols = actionColumns(config);
      const boxes: Box[] = [...header.slots, ...(header.menu ? [header.menu] : [])];
      for (let i = 0; i < boxes.length; i++) {
        for (let j = i + 1; j < boxes.length; j++) expect(overlap(boxes[i], boxes[j])).toBe(false);
        for (const s of header.states) expect(overlap(boxes[i], s.box)).toBe(false);
        // inside the label column
        expect(boxes[i].x).toBeGreaterThanOrEqual(config.inset);
        expect(boxes[i].x + boxes[i].width).toBeLessThanOrEqual(config.inset + config.labelWidth + 1e-9);
      }
      // the slots in their columns: the k-th at the columns' x plus k slots, the "⋯" at the column's end
      header.slots.forEach((s, k) => expect(s.x).toBe(cols.x + k * config.header.slotPx));
      if (header.menu) expect(header.menu.x).toBe(config.inset + config.labelWidth - config.header.slotPx);
      // a short track: no slots, its states marks
      if (height < config.header.minHeightPx) {
        expect(boxes).toEqual([]);
        expect(header.states.every((s) => s.mark)).toBe(true);
      }
      // every state laid out once, in its order
      expect(header.states.map((s) => s.name)).toEqual(states.map((s) => s.name));
      // the label's room ends before the states beside it
      if (header.label.width > 0)
        for (const s of header.states)
          if (s.y <= header.label.y + header.label.height)
            expect(s.x).toBeGreaterThanOrEqual(header.label.x + header.label.width);
    }
  });

  it('a slot’s x is one column for every track of a view: labels long and short, states or none', async () => {
    const { f } = await frame([], Host);
    const s = f.layout.sections[0];
    const eyes = s.tracks.filter((t) => t.header!.slots.length).map((t) => t.header!.slots[0].x);
    expect(eyes.length).toBe(3);
    expect(new Set(eyes).size).toBe(1);
  });
});

describe('the states', () => {
  it('beside the label on a track of one line, on a line of their own on a track of two, marks on a short one', async () => {
    const { f } = await frame([], Host);
    const [seq, dom, thin, more] = f.layout.sections[0].tracks;
    expect(seq.header!.states).toEqual([]);
    expect(more.header!.states).toEqual([]);
    // two lines high: the chip under the label, the label whole
    const chip = dom.header!.states[0];
    expect(chip).toMatchObject({ name: 'seq:filtered', text: '2 filtered', mark: false, x: f.layout.config.inset });
    expect(chip.y).toBeGreaterThan(dom.labelY);
    // short: no slots, the state a mark (its words the card's)
    expect(thin.header!.short).toBe(true);
    expect(thin.header!.slots).toEqual([]);
    expect(thin.header!.states[0]).toMatchObject({ name: 'seq:filtered', mark: true, count: 2 });
  });

  it('no headers given: the built-in states alone, no slots and no “⋯”; the label as before', async () => {
    const { f } = await frame();
    for (const t of f.layout.sections[0].tracks) {
      expect(t.header!.slots).toEqual([]);
      expect(t.header!.menu).toBeUndefined();
      expect(t.header!.actions).toBeUndefined();
    }
    expect(f.layout.sections[0].tracks[1].header!.states[0].text).toBe('2 filtered');
  });
});

describe('the section’s header line (EC-237: hidden tracks kept in the spec, counted where they were)', () => {
  it('every hidden track is counted once in its section’s chip; the line only while one is hidden (`auto`)', async () => {
    const ids = ['seq', 'dom', 'thin', 'more'];
    for (let mask = 0; mask < 16; mask++) {
      const hidden = ids.filter((_, i) => mask & (1 << i));
      const { f } = await frame(hidden);
      const s = f.layout.sections[0];
      expect(s.tracks.length).toBe(4 - hidden.length);
      if (!hidden.length) {
        expect(s.header).toBeUndefined();
        continue;
      }
      expect(s.header!.states).toHaveLength(1);
      expect(s.header!.states[0]).toMatchObject({ name: 'seq:hidden-tracks', count: hidden.length });
      expect(s.header!.states[0].text).toBe(`${hidden.length} hidden`);
      // the line above the first track, `header.sectionPx` high; a section with every track hidden keeps it
      expect(s.header!.line.height).toBe(f.layout.config.header.sectionPx);
      if (s.tracks.length) expect(s.tracks[0].well.y).toBeGreaterThanOrEqual(s.header!.line.y + s.header!.line.height);
      expect(s.height).toBeGreaterThanOrEqual(s.header!.line.height);
    }
  });

  it('`section.header`: `always` draws a line for every section, `none` never', () => {
    const config = (h: 'always' | 'none' | 'auto') => renderConfig({ section: { header: h } });
    const empty = { section: {} as never, filtered: [], hiddenLayers: [], hiddenTracks: [] };
    expect(builtInHeaderStates(empty)).toBeUndefined();
    const line = sectionHeaderLine(undefined, { address: '/v/s' } as never, 's', 30, config('always'), Text, false);
    expect(line).toMatchObject({
      line: { y: 30, height: 20 },
      slots: [],
      states: [],
      hit: { target: 'section-header' },
    });
    const hidden: HeaderParts = { states: [{ name: 'seq:hidden-tracks', text: '1 hidden' }] };
    expect(drawsSectionHeader(undefined, config('always'))).toBe(true);
    expect(drawsSectionHeader(undefined, config('auto'))).toBe(false);
    expect(drawsSectionHeader(hidden, config('auto'))).toBe(true);
    expect(drawsSectionHeader(hidden, config('none'))).toBe(false);
  });
});

describe('drawn alike, hit as laid out', () => {
  it('the draw lists and their SVG draw the states, and the actions of the hovered track only (two consumers)', async () => {
    const { s, f, input } = await frame(['more'], Host);
    const sec = f.layout.sections[0];
    const dom = sec.tracks[1];
    const slot = dom.header!.slots[0];
    // at rest: the states drawn (the section's chip, the track's), no slot
    const rest = recordDrawLists(f.lists, ['labels', 'tip']);
    expect(rest).toContain('"1 hidden"');
    expect(rest).toContain('"2 filtered"');
    expect(rest).not.toMatch(/^pass tip[\s\S]*path/m);
    // a hover on the slot: the actions drawn over the label column in the tip pass, the hovered slot in the accent
    const hit = hitView(f.layout, slot.x + slot.width / 2, slot.y + slot.height / 2)!;
    expect(hit).toMatchObject({ target: 'header-slot', state: 'seq:hide', layer: dom.track.address, section: 'p' });
    expect(hit).toBe(slot.hit);
    const hovered = s.renderer.frame({ ...input, emphasis: { hover: hit } });
    const tip = recordDrawLists(hovered.lists, ['tip']);
    expect(tip.match(/^path /gm)).toHaveLength(1);
    expect(tip).toContain('"C"');
    const svg = drawListsToSvg(hovered.lists);
    expect(svg).toContain('2 filtered');
    expect(svg).toContain('1 hidden');
    expect(svg.match(/<path /g)?.length).toBeGreaterThan(0);
    // the other tracks' actions are not drawn: one backdrop, two slots and three dots
    expect(tip.match(/^rect /gm)).toHaveLength(1 + 3);
  });

  it('a hover on a track’s data shows its actions; the cursor on its header and `emphasis.header` too; `always` and `never`', async () => {
    const { s, f, input } = await frame([], Host);
    const t = f.layout.sections[0].tracks[0];
    const tipPaths = (e: Emphasis) => recordDrawLists(s.renderer.frame({ ...input, emphasis: e }).lists, ['tip']);
    const data = hitView(f.layout, t.box.x + 40, t.box.y + t.box.height / 2)!;
    expect(data.layer.startsWith(t.track.address)).toBe(true);
    expect(tipPaths({ hover: data })).toMatch(/^path /m);
    expect(tipPaths({ cursor: { system: '/data/p', position: 0, part: 'header', track: t.track.address } })).toMatch(
      /^path /m,
    );
    expect(tipPaths({ header: t.track.address })).toMatch(/^path /m);
    expect(tipPaths({})).not.toMatch(/^path /m);
    // `always`: every track's actions at rest; `never`: none laid out (the menu only), the states beside the labels
    const at = async (show: 'always' | 'never') => {
      const r = await settle(spec(), [], { labelColumn: { min: 140 }, header: { actions: { show } } });
      return r.frame({ ...input, emphasis: {} });
    };
    const always = await at('always');
    expect(recordDrawLists(always.lists, ['tip']).match(/^path /gm)).toHaveLength(3);
    const never = await at('never');
    expect(recordDrawLists(never.lists, ['tip'])).not.toMatch(/^(path|rect) /m);
    for (const x of never.layout.sections[0].tracks) expect(x.header!.slots).toEqual([]);
  });

  it('a state’s hit names its entry, its track and its count; the section’s line is `section-header`', async () => {
    const { f } = await frame(['more'], Host);
    const sec: SectionLayout = f.layout.sections[0];
    for (const t of sec.tracks)
      for (const st of t.header!.states) {
        const h = hitView(f.layout, st.box.x + st.box.width / 2, st.box.y + st.box.height / 2);
        expect(h).toMatchObject({ target: 'header-state', state: st.name, layer: t.track.address, count: st.count });
      }
    const chip = sec.header!.states[0];
    expect(hitView(f.layout, chip.box.x + 1, chip.box.y + chip.box.height / 2)).toMatchObject({
      target: 'header-state',
      state: 'seq:hidden-tracks',
      layer: '/views/main/p',
      count: 1,
    });
    const line = sec.header!.line;
    expect(hitView(f.layout, line.x + line.width - 30, line.y + line.height / 2)).toBe(sec.header!.hit);
    // the "⋯" names the menu
    const menu = sec.tracks[0].header!.menu!;
    expect(hitView(f.layout, menu.x + 2, menu.y + menu.height / 2)).toMatchObject({ state: HeaderMenuName });
  });

  it('EC-146: a header’s hover lights none of the track’s items: the hover and selection passes draw nothing for it', async () => {
    const { s, f, input } = await frame([], Host);
    const slot = f.layout.sections[0].tracks[1].header!.slots[0];
    const none = recordDrawLists(s.renderer.frame({ ...input, emphasis: {} }).lists, ['hover', 'highlight']);
    const over = recordDrawLists(s.renderer.frame({ ...input, emphasis: { hover: slot.hit } }).lists, [
      'hover',
      'highlight',
    ]);
    expect(over).toBe(none);
  });

  it('an export draws the states, never the slots; `header.states.export: none` leaves the states out', async () => {
    const out = await renderSvg(spec(['more']), {
      size: { width: 720, height: 200 },
      headers: Host,
      emphasis: { header: '/views/main/p/dom' },
    });
    expect(out.svg).toContain('2 filtered');
    expect(out.svg).toContain('1 hidden');
    expect(recordDrawLists(out.lists, ['tip'])).not.toMatch(/^(path|rect) /m);
    const bare = await renderSvg(spec(['more']), {
      size: { width: 720, height: 200 },
      config: renderConfig({ header: { states: { export: 'none' } } }),
    });
    expect(bare.svg).not.toContain('2 filtered');
    expect(bare.svg).not.toContain('1 hidden');
  });
});
