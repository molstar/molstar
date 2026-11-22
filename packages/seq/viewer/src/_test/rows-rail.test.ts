/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows rail (_spec/README.md "Input", Navigators, The rows rail): a rows band's default navigator, the ruler rail
// upright at the band's right edge beside its last numbers, headless on a fake clock in the platform's event order: a
// drag on its window scrolls the rows (`bandTop`), a grip zooms them (`rowScale`) about the other end, `]` and `[` zoom
// too; a click outside centres the rows; Esc cancels; idle only while every row is shown at the zoom-in limit; drawn
// beside a declared overview unless `navigator.rowsRail` says otherwise.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, overviewNavigator } from '@molstar/seqx-core/lang';
import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { RowsRail, type RowsLay, SlimRailPx, type Rail } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAA';
const RowsTrack = '/views/v/s/rows';

/** An alignment of `n` rows 40 columns wide in a band `height` px tall; an overview over it when asked. */
function spec(n: number, options: { height?: number; overview?: boolean } = {}) {
  const b = createSeqSpec({ title: 'rows rail' });
  const label = Array.from({ length: n }, (_, i) => `r${i}`);
  const letters = Array.from({ length: n }, (_, i) => Letters.slice(i % 40, (i % 40) + 40));
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label, letters } });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '', height: options.height ?? 100 }).layer('r', { type: 'rows', data: aln });
    if (options.overview) v.navigator('nav', overviewNavigator(s));
  });
  return b.build();
}

async function make(root: Node, config: readonly (readonly [unknown, unknown])[] = []) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0], ...config] as never,
  });
  viewers.push(viewer);
  await viewer.load(root);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

// the gesture tests below read exact rows from px: no inset, no floor offset
const Flush = [[SeqCoreConfig.render, { squircleMax: 0 }]] as const;
const band = (v: SeqViewer) => v.view.band(RowsTrack);
const rowsRail = (v: SeqViewer): Rail | undefined => v.frame!.layout.sections[0].rowsRails?.[0];
/** The page y of display row `u` on the rail. */
const yOf = (r: Rail, u: number) => r.x + (u / r.total) * r.width;

describe('the rows rail', () => {
  it('upright at the band’s right edge, beside the numbers, outside the data; its window the rows in view', async () => {
    const v = await make(spec(60));
    const r = rowsRail(v)!;
    expect(r).toBeDefined();
    expect(r.vertical).toBe(true);
    expect(r.active).toBe(true);
    expect(r.bottom - r.top).toBe(SlimRailPx);
    const inset = v.context.renderConfig.inset;
    expect(r.bottom).toBe(800 - inset);
    // beside the data and its last numbers, never over them
    const data = v.frame!.layout.sections[0].tracks[1].box;
    expect(r.top).toBeGreaterThanOrEqual(data.x + data.width);
    // right of the last numbers, as clipped to the band (the widest in view)
    const lay = v.frame!.layout.sections[0].tracks[1].layers.find((l) => l.type.name === 'rows')!.lay as RowsLay;
    expect(lay.lastText.length).toBeGreaterThan(0);
    expect(lay.numberRightX + lay.lastMaxW).toBeLessThanOrEqual(r.top);
    // its ends inset by the well's corner radius
    expect(r.x - r.inset).toBeCloseTo(data.y, 6);
    expect(r.width + 2 * r.inset).toBeCloseTo(data.height, 6);
    expect(r.total).toBe(60);
    expect(r.start).toBe(0);
    expect(r.end).toBeCloseTo(100 / 20, 6);
    // the window floored (5 rows of 60 on 100 px): its drawn middle is body
    const mid = r.x + v.context.renderConfig.navigator.window.minPx / 2;
    expect(v.hit(r.top + 6, mid)?.navigator).toBe(RowsRail);
    expect(v.hit(r.top + 6, mid)?.part).toBe('body');
    expect(v.hit(r.top + 6, yOf(r, 40))?.part).toBe('outside');
  });

  it('a drag on its window scrolls the rows, clamped at the last; a far drag is one step; Esc puts them back', async () => {
    const v = await make(spec(60));
    const r = rowsRail(v)!;
    const x = r.top + 6;
    // the drawn (floored) window's middle; the drag maps by the rows it travels
    const y0 = r.x + v.context.renderConfig.navigator.window.minPx / 2;
    expect(v.hit(x, y0)?.part).toBe('body');
    expect(v.input.down(x, y0)).toBe(true);
    v.input.move(x + 30, y0 + (10 / r.total) * r.width);
    expect(band(v)?.bandTop).toBeCloseTo(10, 6);
    expect(v.input.cursor).toBe('grabbing');
    // past the end: the last rows in view
    v.input.move(x, yOf(r, 200));
    expect(band(v)?.bandTop).toBeCloseTo(55, 6);
    v.input.up();
    expect(v.view.history.value.back).toBe(1);
    v.view.back();
    expect(band(v)?.bandTop ?? 0).toBe(0);
    // Esc during a drag: back where the press found it, no step
    v.draw();
    v.input.down(x, y0);
    v.input.move(x, y0 + 20);
    expect(band(v)?.bandTop).toBeGreaterThan(0);
    v.input.key({ key: 'Escape' });
    expect(band(v)?.bandTop ?? 0).toBe(0);
  });

  it('a grip zooms the rows about the other end; `]` and `[` zoom them too', async () => {
    const v = await make(spec(60), Flush);
    v.view.setBand(RowsTrack, { bandTop: 10 });
    v.draw();
    let r = rowsRail(v)!;
    expect(r.start).toBeCloseTo(10, 6);
    expect(r.end).toBeCloseTo(15, 6);
    const x = r.top + 6;
    // the end grip (a tab just past the window's end) dragged 5 rows down: 10 rows in view, the top kept
    const end = yOf(r, 15) + 2;
    expect(v.hit(x, end)?.part).toBe('end');
    v.input.down(x, end);
    v.input.move(x, end + (yOf(r, 20) - yOf(r, 15)));
    expect(band(v)?.bandTop).toBeCloseTo(10, 6);
    expect(band(v)?.rowScale).toBeCloseTo(100 / 10, 6);
    v.input.up();
    // a zoom is a step
    expect(v.view.history.value.back).toBe(1);
    v.draw();
    r = rowsRail(v)!;
    expect(r.end).toBeCloseTo(20, 6);
    // the start grip dragged 4 rows down: 6 rows in view, the end kept
    const start = yOf(r, 10) - 2;
    expect(v.hit(x, start)?.part).toBe('start');
    v.input.down(x, start);
    v.input.move(x, start + (4 / r.total) * r.width);
    v.input.up();
    expect(band(v)?.rowScale).toBeCloseTo(100 / 6, 6);
    expect(band(v)?.bandTop).toBeCloseTo(14, 6);
    // `]` zooms in, `[` out (the rows behavior's keys), about the band's middle row
    v.draw();
    const before = band(v)!.rowScale!;
    expect(v.input.key({ key: '[' })).toBe(true);
    expect(band(v)?.rowScale).toBeCloseTo(before / 2, 6);
    v.draw();
    expect(v.input.key({ key: ']' })).toBe(true);
    expect(band(v)?.rowScale).toBeCloseTo(before, 6);
  });

  it('with the inset (the well’s corner radius): a press in the inset zone maps to the end; a body drag maps unit-exact', async () => {
    const v = await make(spec(60));
    v.view.setBand(RowsTrack, { bandTop: 10 });
    v.draw();
    let r = rowsRail(v)!;
    expect(r.inset).toBeGreaterThan(0);
    const x = r.top + 6;
    // a body drag of 10 rows' px moves the top exactly 10 rows
    const mid = yOf(r, 12.5);
    v.input.down(x, mid);
    v.input.move(x, mid + (yOf(r, 20) - yOf(r, 10)));
    expect(band(v)?.bandTop).toBeCloseTo(20, 6);
    v.input.up();
    v.draw();
    // a click past the inset axis's end, inside the old extent: the rows centre at the end (clamped to the last rows)
    r = rowsRail(v)!;
    const past = r.x + r.width + r.inset / 2;
    v.input.down(x, past);
    v.input.up();
    v.input.click(x, past);
    v.draw();
    r = rowsRail(v)!;
    expect(r.end).toBeCloseTo(60, 6);
  });

  it('a double-click on its window resets the rows’ zoom (a jump, the top row kept); outside it nothing', async () => {
    const v = await make(spec(60));
    v.view.setBand(RowsTrack, { bandTop: 10, rowScale: 4 });
    v.draw();
    const r = rowsRail(v)!;
    const x = r.top + 6;
    // outside the window: nothing
    expect(v.input.doubleClick(x, yOf(r, r.end + 10))).toBe(true);
    expect(band(v)?.rowScale).toBe(4);
    // on the window: the default zoom, the band's top row kept, one step
    const steps = v.view.history.value.back;
    expect(v.input.doubleClick(x, yOf(r, (r.start + r.end) / 2))).toBe(true);
    expect([band(v)?.rowScale, band(v)?.bandTop]).toEqual([undefined, 10]);
    expect(v.view.history.value.back).toBe(steps + 1);
    v.view.back();
    expect(band(v)?.rowScale).toBe(4);
  });

  it('a double-click on it gives back nothing a click elsewhere cleared', async () => {
    const v = await make(spec(60));
    const r = rowsRail(v)!;
    const data = v.frame!.layout.sections[0].tracks[1].box;
    const cell: [number, number] = [data.x + 20, data.y + 10];
    const click = ([x, y]: [number, number], detail: number) => {
      v.input.down(x, y);
      v.input.up();
      return v.input.click(x, y, { detail });
    };
    click(cell, 1);
    expect(v.selection.value.value.shown.parts.length).toBeGreaterThan(0);
    click(cell, 1);
    expect(v.selection.value.value.shown.parts.length).toBe(0);
    const rail: [number, number] = [r.top + 6, yOf(r, 40)];
    click(rail, 1);
    click(rail, 2);
    v.input.doubleClick(...rail);
    expect(v.selection.value.value.shown.parts.length).toBe(0);
  });

  it('a click outside the window centres the rows there (a jump: back returns); a grip’s click does not', async () => {
    const v = await make(spec(60));
    const r = rowsRail(v)!;
    const x = r.top + 6;
    const y = yOf(r, 40);
    v.input.down(x, y);
    v.input.up();
    expect(v.input.click(x, y, { detail: 1 })).toBe(true);
    expect(band(v)?.bandTop).toBeCloseTo(40 - 2.5, 6);
    expect(v.selection.value.value.shown.parts.length).toBe(0);
    v.view.back();
    expect(band(v)?.bandTop ?? 0).toBe(0);
    v.draw();
    const grip = yOf(rowsRail(v)!, 5) + 2;
    v.input.down(x, grip);
    v.input.up();
    v.input.click(x, grip, { detail: 1 });
    expect(band(v)?.bandTop ?? 0).toBe(0);
  });

  it('idle only when every row is shown and zoom-in is at its limit: a quiet line, no target', async () => {
    // two rows in a 100 px band: shown whole even at the zoom-in limit (40 px a row), the window has nowhere to go
    const v = await make(spec(2));
    const r = rowsRail(v)!;
    expect(r.whole).toBe(true);
    expect(r.active).toBe(false);
    const x = r.top + 6;
    const y = r.x + r.width / 2;
    expect(v.hit(x, y)?.navigator).toBeUndefined();
    v.input.down(x, y);
    expect(v.input.dragHit).toBeUndefined();
    v.input.up();
    // four rows at 20 px: every row shown, but the rows can zoom in (the window can shrink): live, its grips there
    const w = await make(spec(4));
    const four = rowsRail(w)!;
    expect(four.start <= 0 && four.end >= four.total).toBe(true);
    expect(four.whole).toBe(false);
    expect(four.active).toBe(true);
    expect(w.hit(four.top + 6, yOf(four, four.end) - 1)?.part).toBe('end');
  });

  it('a grip dragged to the band’s end shows every row, and dragged back the window shrinks again', async () => {
    const v = await make(spec(60), Flush);
    let r = rowsRail(v)!;
    const x = r.top + 6;
    /**
     * The page y of the end grip within the band (a floored window's is past its true end; a whole one's tab hangs past
     * the band, its inner quarter in it): the last px in the band that hits it.
     */
    const endGrip = (rail: Rail) => {
      let y = Number.NaN;
      for (let p = rail.x; p < rail.x + rail.width; p += 0.5) if (v.hit(x, p)?.part === 'end') y = p;
      return y;
    };
    // the end grip to the rail's end: every row in view
    const end = endGrip(r);
    expect(Number.isNaN(end)).toBe(false);
    v.input.down(x, end);
    v.input.move(x, r.x + r.width + 20);
    v.input.up();
    v.draw();
    r = rowsRail(v)!;
    expect(r.start).toBeCloseTo(0, 6);
    expect(r.end).toBeCloseTo(60, 6);
    // still live, its grips there: the rows can zoom in
    expect(r.active).toBe(true);
    const back = endGrip(r);
    expect(Number.isNaN(back)).toBe(false);
    v.input.down(x, back);
    v.input.move(x, back + (yOf(r, 30) - yOf(r, 60)));
    v.input.up();
    v.draw();
    r = rowsRail(v)!;
    expect(r.end - r.start).toBeLessThan(60);
    expect(r.end).toBeCloseTo(30, 0);
    expect(r.active).toBe(true);
  });

  it('beside a declared overview by default (`navigator.rowsRail: always`); replaced by it under `unlessOverview`; none under `none`', async () => {
    const both = await make(spec(60, { overview: true }));
    const r = rowsRail(both)!;
    expect(r).toBeDefined();
    expect(r.active).toBe(true);
    expect(r.total).toBe(60);
    // the rail and the overview draw one window: a drag on the rail moves the overview's in the same frame
    const nav = both.navigatorModel('nav')!;
    nav.resize(0, 200);
    expect(nav.draw()!.windows[0].r0).toBeCloseTo(0, 6);
    const y0 = r.x + both.context.renderConfig.navigator.window.minPx / 2;
    both.input.down(r.top + 6, y0);
    both.input.move(r.top + 6, y0 + (10 / r.total) * r.width);
    both.input.up();
    both.draw();
    expect(band(both)?.bandTop).toBeCloseTo(10, 6);
    expect(nav.draw()!.windows[0].r0).toBeCloseTo(10, 6);
    const v = await make(spec(60, { overview: true }), [[viewerItem('navigator.rowsRail'), 'unlessOverview']]);
    expect(v.frame!.layout.sections[0].rowsRails).toBeUndefined();
    const none = await make(spec(60), [[viewerItem('navigator.rowsRail'), 'none']]);
    expect(none.frame!.layout.sections[0].rowsRails).toBeUndefined();
  });

  it('the data keeps its width when it is hidden', async () => {
    const w = await make(spec(60));
    expect(rowsRail(w)).toBeDefined();
    const with_ = w.frame!.layout.sections[0].tracks[1].box.width;
    w.setRailsShown(false);
    w.draw();
    expect(w.frame!.layout.sections[0].rowsRails).toBeUndefined();
    expect(w.frame!.layout.sections[0].tracks[1].box.width).toBe(with_ + SlimRailPx);
  });
});
