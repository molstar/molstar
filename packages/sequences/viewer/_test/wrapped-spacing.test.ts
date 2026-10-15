/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text layout's spacing (seqspec:viewer#the-text-view): the render props `wrappedGroupGap` (between groups),
// `wrappedGap` (between lines), `wrappedTrackGap` (letters to the first annotation lane) and `wrappedLaneHeight` (a
// lane under the text) place what they name; `compact` is the dense preset, shorter than `wrapped` on case 1. A line
// holds whole groups alike: a short last group goes on a line of its own where a full group would not fit.
// `wrappedAlign: center` centers the block of lines by its widest line in the area it ends with, and `wrappedMaxWidth`
// wraps the lines narrower than the area. The lanes per layer kind: each band's height
// from the render props (letters and emphasis the font's, underlines `wrappedTrackGap` + a `wrappedLaneHeight` a lane, a
// strip `wrappedTrackGap` + `wrappedStripHeight`), `compact`'s own; no band for a kind that draws nothing in text; a
// document of letters alone spaced as before. Chips: a structure's ligands flow after the chains on shared rows as wide
// as their labels, wrapping as words do, at the lines' spacing, with no group labels' band.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { case1 } from '../../core/_test/format/cases.js';
import type { ShownIndex } from '../src/axis.js';
import { chipFontPx, groupLabelBand } from '../src/runs.js';
import { DefaultPresets, EntryKind, Hit, KernelRenderers, type LayerRenderer, type RenderProps } from '../src/index.js';
import { opened } from './layout-fixtures.js';
import { chipsDoc } from './prepare-fixtures.js';

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

describe('the lanes per layer kind', () => {
  /** 60 residues as text: letters, an underlined region, values, color, and the kinds that do not draw in text. */
  function kinds(): SeqViewSpecDoc {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(2) });
    const r = b.data.annotation({ ref: 'r', geometry: 'ranges', target: a, value: { start: [3], stop: [30] } });
    const value = Array.from({ length: 60 }, (_, i) => (i % 10) / 10);
    const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } });
    const p = b.data.annotation({ ref: 'p', geometry: 'pairs', target: a, value: { start: [0], b_start: [20] } });
    const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
    v.track({ ref: 'v.ruler' }).ruler();
    const letters = v.track({ ref: 'v.letters' });
    letters.letters();
    letters.emphasis({ data: r });
    v.track({ ref: 'v.r' }).ranges({ data: r, display: 'underline' } as never);
    v.track({ ref: 'v.values' }).values({ data: d });
    v.track({ ref: 'v.color' }).color({ data: d });
    // Color beside letters on one track: their cells, in the letters' band.
    const behind = v.track({ ref: 'v.behind' });
    behind.letters();
    behind.color({ data: d });
    v.track({ ref: 'v.plot' }).layer({ type: 'dot_plot', data: p } as never);
    return b.getDocument();
  }
  /** The first line's bands as [track, height]. */
  const firstLine = (t: Awaited<ReturnType<typeof open>>) => {
    const l = t.layout('v.main');
    return t
      .entries('v.main')
      .filter((e, k) => e[0] === T && l.entryRow[k] === 0)
      .map((e) => [e[1], e[4]]);
  };

  it("each kind's band from the render props; none for a kind that draws nothing in text", async () => {
    const t = await open(kinds(), 1000);
    // Letters (with their emphasis) the font's 12 × 1.6; the underline 2 + 5; a strip 2 + 3; the ruler and the dot
    // plot nothing (text omits them).
    expect(firstLine(t)).toEqual([
      ['v.letters', 20],
      ['v.r', 7],
      ['v.values', 5],
      ['v.color', 5],
      ['v.behind', 20],
    ]);
    const custom = await open(kinds(), 1000, { wrappedTrackGap: 3, wrappedLaneHeight: 6, wrappedStripHeight: 4 });
    expect(firstLine(custom).map((b) => b[1])).toEqual([20, 9, 7, 7, 20]);
  });

  it('values alone on its track is drawn in text, with no fallback reported', async () => {
    const t = await open(kinds(), 1000);
    const m = t.ctx.graph.assembled.value?.view;
    const types = Array.from(m?.resolved.layers.type ?? [], (x) => (x < 0 ? '' : m?.set.strings.get(x)));
    expect(types).toContain('values');
    expect(t.ctx.graph.issues.value.filter((i) => i.message.includes('"values"'))).toEqual([]);
  });

  it('compact sets the strip and the lanes tighter', async () => {
    const t = await open(kinds(), 1000);
    t.ctx.ui.setStyle('v.main', 'compact');
    await t.ctx.graph.update();
    const props = t.ctx.layout.props('v.main');
    expect([props.wrappedTrackGap, props.wrappedLaneHeight, props.wrappedStripHeight]).toEqual([0, 4, 2]);
    expect(firstLine(t).map((b) => b[1])).toEqual([20, 4, 2, 2, 20]);
  });

  it('a document of letters alone is spaced as before: one band of the letters a line, lines 4 px apart', async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'a', alphabet: 'protein', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(4) });
    b.view({ ref: 'v', label: 'V', layout: 'text' }).track({ ref: 'v.letters' }).letters();
    const t = await open(b.getDocument(), 400);
    const e = t.entries('v.main');
    expect(e.every((x) => x[0] === T && x[1] === 'v.letters' && x[4] === 20)).toBe(true);
    expect(e[1][3] - (e[0][3] + e[0][4])).toBe(4);
  });
});

describe('chips: short units in the text', () => {
  const LIGANDS = ['HEM', 'ERY', 'NI', 'ATP', 'NAG', 'SO4'];
  /** Each line as [block, left, top]. */
  const rows = (t: Awaited<ReturnType<typeof open>>) => {
    const l = t.layout('v.main').lines;
    return Array.from({ length: l.count }, (_, k) => [l.block[k], l.left[k], l.top[k], l.follows[k]]);
  };
  const chip = (t: Awaited<ReturnType<typeof open>>, id: string) =>
    Math.ceil(t.ctx.text.label(id, t.ctx.text.font('ui', chipFontPx(12), 600)).width) + 10;

  it("flow after the chain on one row, each as wide as its label and the group gap apart, at the lines' spacing", async () => {
    const t = await open(chipsDoc(LIGANDS), 1000, { wrappedGroupLabels: 'aboveEnd' });
    const r = rows(t);
    // The chain's line, then one row of chips (the 12 px letters: chips' labels at 10 px).
    expect(r.length).toBe(1 + LIGANDS.length);
    expect(r.slice(1).every((x) => x[2] === r[1][2])).toBe(true);
    expect(r.slice(2).every((x) => x[3] === 1)).toBe(true);
    let x = 0;
    for (let k = 0; k < LIGANDS.length; k++) {
      expect(r[1 + k][1]).toBe(x);
      x += chip(t, LIGANDS[k]) + 8;
    }
    // The chips follow the chain as a new line does: the lines' gap (4 px) under the last line's bottom band, then the
    // group labels' band's room, empty (no labels of their own).
    const l = t.layout('v.main');
    expect(r[1][2]).toBe(l.lines.bottom[0] + 4);
    const labels = t.entries('v.main').filter((e) => e[0] === EntryKind.labels);
    expect(labels.length).toBe(1);
    // Each chip's entry on its row, as tall as the letters', where the line's would start.
    const chips = t.entries('v.main').filter((e, k) => e[0] === T && l.entryRow[k] >= 1);
    const band = groupLabelBand(12);
    expect(chips.map((e) => [e[3], e[4]])).toEqual(LIGANDS.map(() => [r[1][2] + band, 20]));
    // As a new line: from the last line's bottom, the lines' gap and then what a line has above its letters.
    const letters0 = t.entries('v.main').find((e) => e[0] === T)?.[3] ?? 0;
    expect(chips[0][3] - l.lines.bottom[0]).toBe(4 + (letters0 - l.lines.top[0]));
  });

  it("with lines touching (no gap, `compact`'s): the chips' row still a line's room under the last line", async () => {
    const t = await open(chipsDoc(LIGANDS, true), 1000, { wrappedGroupLabels: 'aboveEnd', wrappedGap: 0 });
    const l = t.layout('v.main');
    const first = t.entries('v.main').findIndex((e, k) => e[0] === T && l.entryRow[k] >= 1);
    expect(l.lines.top[1]).toBeGreaterThanOrEqual(l.lines.bottom[0]);
    expect(t.entries('v.main')[first][3] - l.lines.bottom[0]).toBe(groupLabelBand(12));
  });

  it('wrap as words where the next does not fit; rows share the tallest band of each track', async () => {
    const t = await open(chipsDoc([...LIGANDS, ...LIGANDS], true), 260);
    const r = rows(t);
    const chipRows = new Set(r.filter((x) => x[0] >= 1).map((x) => x[2]));
    expect(chipRows.size).toBeGreaterThan(1);
    const l = t.layout('v.main');
    const width = l.dataRight - l.dataLeft;
    for (let k = 0; k < l.lines.count; k++) {
      if (l.lines.block[k] < 1) continue;
      expect(l.lines.left[k] + chip(t, [...LIGANDS, ...LIGANDS][l.lines.block[k] - 1])).toBeLessThanOrEqual(
        width + 0.01,
      );
      // A row starts at the left; a later chip follows on its row.
      expect(l.lines.left[k] === 0).toBe(l.lines.follows[k] === 0);
    }
    // The underline track has no item on a chip: its gap alone, on every chip of a row.
    const under = t.entries('v.main').filter((e, k) => e[1] === 'v.r' && l.lines.block[l.entryRow[k]] >= 1);
    expect(new Set(under.map((e) => e[4]))).toEqual(new Set([2]));
  });

  it('tracks keep one block a line: no chips there', async () => {
    const doc = chipsDoc(LIGANDS) as { root: { children: { kind: string; params: Record<string, unknown> }[] } };
    const view = doc.root.children.find((c) => c.kind === 'view');
    if (view) view.params.layout = 'tracks';
    const t = await open(doc as unknown as SeqViewSpecDoc, 1000);
    expect(t.layout('v.main').lines.count).toBe(0);
  });
});
