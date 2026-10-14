/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text layout's spacing (seqspec:viewer#text-measuring): the render props `wrappedGroupGap` (between groups),
// `wrappedGap` (between lines), `wrappedTrackGap` (letters to the first annotation lane) and `wrappedLaneHeight` (a
// lane under the text) place what they name; `compact` is the dense preset, shorter than `wrapped` on case 1. A line
// holds whole groups alike: a short last group goes on a line of its own where a full group would not fit.
// `wrappedAlign: center` centers the block of lines by its widest line in the area it ends with, and `wrappedMaxWidth`
// wraps the lines narrower than the area.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { case1 } from '../../core/_test/format/cases.js';
import type { ShownIndex } from '../src/axis.js';
import { DefaultPresets, EntryKind, Hit, KernelRenderers, type LayerRenderer, type RenderProps } from '../src/index.js';
import { opened } from './layout-fixtures.js';

const T = EntryKind.track;
/** Mono 12 px (the letters' cap) in the fixture measurer: 7.2 px a cell. */
const CELL = 7.2;

/** A protein of `n` residues as text, with one region underlined over residues 3 to 30. */
function text(n: number): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const value = 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(Math.ceil(n / 30)).slice(0, n);
  const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value });
  const r = b.data.annotation({ ref: 'r', geometry: 'ranges', target: a, value: { start: [3], stop: [30] } });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
  v.track({ ref: 'v.letters' }).letters();
  v.track({ ref: 'v.r' }).ranges({ data: r, display: 'underline' } as never);
  return b.getDocument();
}

/** The kernel's renderers (their lanes' `lineHeights`), on the built-in presets (`wrapped`, `compact`). */
const BUILT_IN = {
  presets: false as const,
  renderers: [...KernelRenderers] as LayerRenderer[],
  registry: [{ presets: [...DefaultPresets] }],
};

/** Opened on the built-ins, the section `width` px wide. */
function open(doc: SeqViewSpecDoc, width: number, render?: Partial<RenderProps>) {
  return opened(doc, {
    ...BUILT_IN,
    size: { 'v.main': [width, 2000] },
    render: { wrappedGroupLabels: 'none', wrappedGutter: 0, ...render },
  });
}

describe("the text layout's spacing", () => {
  it("each prop places what it names: groups, lines, the letters' gap to the lanes, a lane", async () => {
    const t = await open(text(60), 400, {
      wrappedGroupGap: 11,
      wrappedGap: 7,
      wrappedTrackGap: 3,
      wrappedLaneHeight: 6,
    });
    const view = t.ctx.layout.axis('v.main');
    // Position 10 starts the second group: ten cells and the group gap after the first.
    expect(view && view.x(0, 10 as ShownIndex) - view.x(0, 0 as ShownIndex)).toBeCloseTo(10 * CELL + 11, 4);
    const e = t.entries('v.main');
    // Line 0: the letters, then the region's lane (one lane: 3 + 6); line 1 starts 7 px under it.
    const [letters0, lane0, letters1] = e;
    expect([letters0[0], letters0[1], lane0[1], lane0[4]]).toEqual([T, 'v.letters', 'v.r', 3 + 6]);
    expect(letters1[3]).toBe(lane0[3] + lane0[4] + 7);
  });

  it('a short last group goes on a line of its own where a full group would not fit', async () => {
    // 393 positions; a data area for 13 full groups and 3 letters more (not a 14th full group).
    const groups = 13 * 10 * CELL + 12 * 8;
    const probe = await open(text(393), 1000);
    const lp = probe.layout('v.main');
    const chrome = 1000 - (lp.dataRight - lp.dataLeft);
    const t = await open(text(393), Math.ceil(groups + 8 + 3 * CELL + 2 + chrome));
    const l = t.layout('v.main').lines;
    const lines = Array.from({ length: l.count }, (_, n) => [l.start[n], l.stop[n]]);
    expect(lines).toEqual([
      [0, 130],
      [130, 260],
      [260, 390],
      [390, 393],
    ]);
    // Where a full group fits too, the short one joins the line.
    const wide = await open(text(393), Math.ceil(groups + 8 + 10 * CELL + 2 + chrome));
    expect(wide.layout('v.main').lines.count).toBe(3);
  });

  it("compact is the dense preset: case 1's text is shorter than wrapped's, still in groups of 10", async () => {
    const styled = (style: string): SeqViewSpecDoc => {
      const doc = case1();
      type N = { kind: string; params?: Record<string, unknown>; children?: N[] };
      const walk = (n: N, inText: boolean): N => ({
        ...n,
        params: inText && n.kind === 'section' ? { ...n.params, style } : n.params,
        children: n.children?.map((c) => walk(c, inText || (n.kind === 'view' && n.params?.layout === 'text'))),
      });
      return { ...doc, root: walk(doc.root as N, false) } as SeqViewSpecDoc;
    };
    const height = async (style: string) => {
      const t = await opened(styled(style), { ...BUILT_IN, size: { 'text.main': [800, 4000] } });
      const section = t.ctx.viewports.sections.value[0];
      const l = t.layout(section);
      const last = l.entryCount - 1;
      return { h: l.entryTop[last] + l.entryHeight[last], props: t.ctx.layout.props(section) };
    };
    const wrapped = await height('wrapped');
    const compact = await height('compact');
    expect(compact.props.wrappedGroup).toBe(10);
    expect(compact.props.wrappedGroupGap).toBeLessThan(wrapped.props.wrappedGroupGap);
    expect(compact.h).toBeLessThan(wrapped.h * 0.9);
  });

  /** Mono 12 px groups of 10 with the preset's 8 px gap: `g` groups are `g · 72 + (g − 1) · 8` px. */
  const groupsPx = (g: number) => g * 10 * CELL + (g - 1) * 8;

  it('wrappedAlign center puts the block of lines in the middle of the area; hits and the margins follow', async () => {
    // 60 positions in a 400 px section: a data area of 368 px (8 px inset left, 24 px strip right) holds 4 groups a
    // line (312 px; a fifth would need 392), so two lines, of 40 and 20 positions.
    const start = await open(text(60), 400);
    const centered = await open(text(60), 400, { wrappedAlign: 'center' });
    const a = start.layout('v.main');
    const c = centered.layout('v.main');
    expect([a.dataLeft, a.dataRight]).toEqual([8, 376]);
    expect(Array.from(c.lines.stop.subarray(0, c.lines.count))).toEqual([40, 60]);
    // The block is as wide as its widest line, with the rest of the area split evenly: 28 px either side.
    expect(c.dataRight - c.dataLeft).toBe(groupsPx(4));
    expect([c.dataLeft - a.dataLeft, a.dataRight - c.dataRight]).toEqual([28, 28]);
    // Position 0 is at the centered left edge: a hit there names it, a hit in the margin names nothing.
    const hit = (t: typeof centered, x: number) => {
      const y = c.scrolledTop + c.entryTop[0] + 4;
      const h = Hit.create();
      t.ctx.layout.hitTest('v.main', x, y, h);
      return h.shown;
    };
    const view = centered.ctx.layout.axis('v.main');
    expect(view?.x(0, 0 as ShownIndex)).toBeCloseTo(c.dataLeft, 4);
    expect(hit(centered, c.dataLeft + CELL / 2)).toBe(0);
    expect(hit(centered, a.dataLeft + 14)).toBe(-1);
  });

  it('wrappedMaxWidth wraps the lines there, and centered puts that column in the middle', async () => {
    // 120 positions in a 600 px section: the area (568 px) holds 7 groups a line; a column of 330 px holds 4 (312 px).
    const wide = await open(text(120), 600, { wrappedAlign: 'center' });
    const column = await open(text(120), 600, { wrappedAlign: 'center', wrappedMaxWidth: 330 });
    const w = wide.layout('v.main');
    const c = column.layout('v.main');
    expect(Array.from(w.lines.stop.subarray(0, w.lines.count))).toEqual([70, 120]);
    expect(Array.from(c.lines.stop.subarray(0, c.lines.count))).toEqual([40, 80, 120]);
    expect(c.dataRight - c.dataLeft).toBe(groupsPx(4));
    // Centered in the whole area: (568 − 312) / 2 = 128 px either side of the column.
    const s0 = (await open(text(120), 600)).layout('v.main');
    expect([c.dataLeft - s0.dataLeft, s0.dataRight - c.dataRight]).toEqual([128, 128]);
  });

  it('centered with a vertical overview that the text does not need: centered in the area it ends with', async () => {
    // `verticalNav: overview` reserves the overview while the lines are laid out; text that fits needs no overview,
    // so the area is the whole width less the inset and the strip, and the margins are measured in it.
    const props = { verticalNav: 'overview', verticalOverview: 40 } as const;
    const t = await open(text(60), 400, { ...props, wrappedAlign: 'center' });
    const l = t.layout('v.main');
    expect(l.vscroll).toBe(false);
    expect(l.vnavWidth).toBe(0);
    const right = 400 - 24;
    expect(l.dataRight - l.dataLeft).toBe(groupsPx(4));
    expect([l.dataLeft - 8, right - l.dataRight]).toEqual([28, 28]);
  });

  it('the default is start: nothing moves unless asked', async () => {
    const t = await open(text(60), 400);
    expect(t.ctx.layout.props('v.main').wrappedAlign).toBe('start');
  });
});
