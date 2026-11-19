/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows rail (_spec/README.md "Input", Navigators, The rows rail): a rows band's default navigator, the ruler rail
// upright at the band's right edge beside its last numbers, headless on a fake clock in the platform's event order: a
// drag on its window scrolls the rows (`bandTop`), a grip zooms them (`rowScale`) about the other end, `]` and `[` zoom
// too; a click outside centres the rows; Esc cancels; idle while every row is shown; a declared overview replaces it.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, overviewNavigator } from '@molstar/seqx-core/lang';
import { RowsRail, SlimRailPx, type Rail } from '@molstar/seqx-core/render';
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

async function make(root: Node) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0]] as never,
  });
  viewers.push(viewer);
  await viewer.load(root);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

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
    expect(r.x).toBeCloseTo(data.y, 6);
    expect(r.width).toBeCloseTo(data.height, 6);
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
    const v = await make(spec(60));
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
    v.input.move(x, end + (5 / r.total) * r.width);
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

  it('idle while every row is shown: a quiet line, no target; live again once the rows are zoomed past the band', async () => {
    const v = await make(spec(4));
    const r = rowsRail(v)!;
    expect(r.whole).toBe(true);
    expect(r.active).toBe(false);
    const x = r.top + 6;
    const y = r.x + r.width / 2;
    expect(v.hit(x, y)?.navigator).toBeUndefined();
    v.input.down(x, y);
    expect(v.input.dragHit).toBeUndefined();
    v.input.up();
    // the rows zoomed past the band: live
    v.view.setBand(RowsTrack, { rowScale: 40 });
    v.draw();
    expect(rowsRail(v)!.active).toBe(true);
  });

  it('a declared overview over the rows replaces it; the data keeps its width when it is hidden', async () => {
    const v = await make(spec(60, { overview: true }));
    expect(v.frame!.layout.sections[0].rowsRails).toBeUndefined();
    const w = await make(spec(60));
    expect(rowsRail(w)).toBeDefined();
    const with_ = w.frame!.layout.sections[0].tracks[1].box.width;
    w.setRailsShown(false);
    w.draw();
    expect(w.frame!.layout.sections[0].rowsRails).toBeUndefined();
    expect(w.frame!.layout.sections[0].tracks[1].box.width).toBe(with_ + SlimRailPx);
  });
});
