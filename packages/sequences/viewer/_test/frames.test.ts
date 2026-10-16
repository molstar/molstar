/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames (seqspec:viewer#frames): named regions outlined across the tracks a `frame` layer spans, in the text (open on
// the side a line break cuts, the title on the line the item starts, on the rows its item targets) and in tracks (over
// the rows' entries, following a row zoom, open where the view cuts); hit as their item on their title's band.

import { describe, expect, it } from 'vitest';
import { Dirty } from '../src/draw/index.js';
import { Hit, HoverKind, type RenderProps } from '../src/index.js';
import { recording } from './composer-fixtures.js';
import { lit } from './overlay-fixtures.js';
import { FIGURE_FRAMES, type FrameRow, figureDoc, sixtyColumns } from './text-rows-fixtures.js';

/** A frame across the first line break (columns 50 to 70), one over rows 2 to 4 alone (a row target each). */
const FRAMES: readonly FrameRow[] = [
  ...FIGURE_FRAMES.slice(0, 1),
  { target: 'globins', start: 50, stop: 70, id: 'f2', label: 'Across', color: '#ff7f0e', line_style: 'dashed' },
  { target: 'globins/2', start: 30, stop: 36, id: 'f3', label: 'Rows', color: '#17becf' },
  { target: 'globins/4', start: 30, stop: 36, id: 'f3', label: 'Rows', color: '#17becf' },
];
/** The segments of `color` in a recording: [x0, y0, x1, y1] (device px), and whether dashed. */
function segments(r: ReturnType<typeof recording>, color: number) {
  const b = r.pathItems;
  const out: { x0: number; y0: number; x1: number; y1: number; dashed: boolean }[] = [];
  for (let p = 0; p < b.n; p++) {
    if (b.color[p] !== color) continue;
    const c = b.cstart[p];
    out.push({ x0: b.coords[c], y0: b.coords[c + 1], x1: b.coords[c + 2], y1: b.coords[c + 3], dashed: b.dash[p] > 0 });
  }
  return out;
}

async function text(render?: Partial<RenderProps>) {
  const probe = await lit(figureDoc(), { width: 1100, height: 1400 });
  const cell = probe.ctx.layout.axis(probe.section)?.axis?.fixedCell ?? 0;
  probe.ctx.dispose();
  const t = await lit(figureDoc({ frames: FRAMES }), {
    width: 1100,
    height: 1400,
    render: render ?? sixtyColumns(cell),
  });
  const l = t.ctx.layout.get(t.section)?.current;
  if (!l) throw new Error('no layout');
  return { ...t, l, cell };
}

describe('frames in the text', () => {
  it('closed where the item ends, open where a line break cuts it; the title on its first line only', async () => {
    const { l, frame } = await text();
    expect([l.lines.start[1], l.lines.count]).toEqual([60, 2]);
    const r = frame(recording(), Dirty.all);
    const lineOf = (y: number) => (y < l.scrolledTop + l.lines.bottom[0] ? 0 : 1);
    const across = segments(r, 0xff7f0e);
    const vertical = across.filter((s) => s.x0 === s.x1);
    const horizontal = across.filter((s) => s.y0 === s.y1);
    // Line 0: its left side (column 50), its top and bottom to the line's end; line 1: its right side (column 70).
    expect(vertical.map((s) => lineOf(s.y0))).toEqual([0, 1]);
    expect(horizontal.map((s) => lineOf(s.y0))).toEqual([0, 0, 1, 1]);
    expect(across.every((s) => s.dashed)).toBe(true);
    const n = segments(r, 0x9467bd);
    expect(n.filter((s) => s.x0 === s.x1).length).toBe(2);
    expect(n.some((s) => s.dashed)).toBe(false);
    const titles = Array.from({ length: r.textItems.n }, (_, i) => r.textItems.text[i]);
    expect(titles.filter((x) => x === 'Across').length).toBe(1);
    expect(titles.filter((x) => x === 'N-segment').length).toBe(1);
  });

  it('an item on rows 2 and 4 frames those rows alone: from row 2 to row 4 of the block', async () => {
    const { l, frame } = await text();
    const tr = l.textRows;
    const top = l.scrolledTop + l.entryTop[l.trackFirst[tr.track]];
    const rows = segments(frame(recording(), Dirty.all), 0x17becf);
    const vertical = rows.filter((s) => s.x0 === s.x1);
    expect(vertical.length).toBe(2);
    for (const s of vertical) {
      expect(s.y0).toBe(Math.round(top + 2 * tr.height));
      expect(s.y1).toBe(Math.round(top + 5 * tr.height));
    }
  });

  it('a frame spans every row of the block (half a track gap outside it); its title band hits its item', async () => {
    const { l, frame, ctx, section } = await text();
    const tr = l.textRows;
    const e = l.trackFirst[tr.track];
    const top = l.scrolledTop + l.entryTop[e];
    const n = segments(frame(recording(), Dirty.all), 0x9467bd).filter((s) => s.x0 === s.x1);
    const gap = ctx.layout.props(section).wrappedTrackGap / 2;
    expect(n[0].y0).toBe(Math.round(top - gap));
    expect(n[0].y1).toBe(Math.round(top + l.entryHeight[e] + gap));
    // The title's band: the frames track's entry on line 0, over the item's columns.
    let band = -1;
    for (let k = 0; k < l.entryCount && band < 0; k++)
      if (l.entryRow[k] === 0 && l.entryHeight[k] > 0 && l.entryTrack[k] === tr.track - 1) band = k;
    expect(band).toBeGreaterThanOrEqual(0);
    const h = Hit.create();
    const x = l.dataLeft + l.lines.left[0] + 5.5 * (ctx.layout.axis(section)?.axis?.fixedCell ?? 0);
    expect(ctx.layout.hitTest(section, x, l.scrolledTop + l.entryTop[band] + l.entryHeight[band] / 2, h)).toBe(true);
    expect(h.kind).toBe(HoverKind.item);
    expect(
      ctx.graph.assembled.value?.set.strings.get(ctx.graph.assembled.value.set.annotations[h.annotation]?.key ?? -1),
    ).toBe('figure.frames');
  });
});

describe('frames in tracks', () => {
  it("over the rows' entries, following a row zoom; open where the view's edge cuts the item", async () => {
    const t = await lit(figureDoc({ style: { preset: 'msa' } }), { width: 1000, height: 600 });
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error('no layout');
    const r = l.rows;
    const at = () => {
      const top = l.scrolledTop + l.entryTop[r.entry] - (t.vp?.live.current.scrollY ?? 0);
      return [Math.round(top), Math.round(top + r.count * r.height)];
    };
    const verticals = () => segments(t.frame(recording(), Dirty.all), 0x9467bd).filter((s) => s.x0 === s.x1);
    const v0 = verticals();
    expect(v0.length).toBe(2);
    expect([v0[0].y0, v0[0].y1]).toEqual(at());
    // A row zoom: the rows' height changes without a build; the frame follows.
    t.vp?.setRowScale(1.5, 0);
    await t.ctx.graph.update();
    const v1 = verticals();
    expect([v1[0].y0, v1[0].y1]).toEqual(at());
    expect(v1[0].y1 - v1[0].y0).toBeGreaterThan(v0[0].y1 - v0[0].y0);
    // Zoomed in on column 10: the frame's start is left of the view, its left side open.
    const cell = 16;
    t.vp?.setView(cell, 10 * cell, 'api');
    const v2 = verticals();
    expect(v2.length).toBe(1);
  });
});
