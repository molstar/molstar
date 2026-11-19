/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cursor says the gesture (_spec/README.md "Input", Navigators): over a window's grip a resize cursor along the
// navigator's axis (`ew-resize` on a horizontal one, `ns-resize` on a vertical one), held through the grip's drag
// wherever the pointer goes; over its body `grab`, `grabbing` while dragged; outside `pointer`; idle the default. One
// rule for the ruler rail, the rows rail, the strip, the overview and the miniature, headless through the models.

import { createSeqSpec, miniatureNavigator, overviewNavigator, stripNavigator } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { partCursor } from '../navigator.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Letters =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQ';
const S = '/views/v/s';
const RowsTrack = '/views/v/s/rows';

type Kind = 'strip' | 'overview' | 'miniature';

function spec(n: number, nav?: Kind) {
  const b = createSeqSpec({ title: 'cursors' });
  const label = Array.from({ length: n }, (_, i) => `r${i}`);
  const letters = Array.from({ length: n }, (_, i) => Letters.slice(i % 40, (i % 40) + 60));
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label, letters } });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '', height: 100 }).layer('r', { type: 'rows', data: aln });
    if (nav === 'strip') v.navigator('nav', stripNavigator(s));
    if (nav === 'overview') v.navigator('nav', overviewNavigator(s));
    if (nav === 'miniature') v.navigator('nav', miniatureNavigator(s));
  });
  return b.build();
}

async function make(n: number, nav?: Kind) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0]] as never,
  });
  viewers.push(viewer);
  await viewer.load(spec(n, nav));
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

async function drawNav(viewer: SeqViewer) {
  const nav = viewer.navigatorModel('nav')!;
  nav.draw();
  await viewer.settle();
  frames!.flush();
  viewer.draw();
  nav.draw();
  return nav;
}

/** The first coordinate along a line where `part` is hit. */
function find(part: string, to: number, hit: (t: number) => string | undefined): number {
  for (let t = 0; t < to; t += 0.5) if (hit(t) === part) return t;
  throw new Error(`no ${part} found`);
}

describe('the cursor by part', () => {
  it('a grip resizes along its axis, the body grabs, outside points', () => {
    expect(partCursor('start', 'positions', false)).toBe('ew-resize');
    expect(partCursor('end', 'positions', true)).toBe('ew-resize');
    expect(partCursor('start', 'rows', false)).toBe('ns-resize');
    expect(partCursor('end', 'rows', true)).toBe('ns-resize');
    expect(partCursor('body', 'positions', false)).toBe('grab');
    expect(partCursor('body', 'rows', true)).toBe('grabbing');
    expect(partCursor('outside', 'positions', false)).toBe('pointer');
    expect(partCursor('outside', 'positions', true)).toBe('default');
  });
});

describe('the cursor says the gesture', () => {
  for (const kind of ['strip', 'overview', 'miniature'] as const) {
    it(`the ${kind}: grips resize along their axis, the body grabs, outside points, held through a drag`, async () => {
      const v = await make(60, kind);
      v.view.setRange(S, 20, 60);
      v.view.setBand(RowsTrack, { bandTop: 10 });
      v.draw();
      const nav = await drawNav(v);
      expect(nav.cursor).toBe('default');
      const horizontal = kind !== 'overview';
      const probes: { part: string; cursor: 'ew-resize' | 'ns-resize'; vertical: boolean }[] = [];
      for (const part of ['start', 'end']) {
        if (horizontal) probes.push({ part, cursor: 'ew-resize', vertical: false });
        if (kind !== 'strip') probes.push({ part, cursor: 'ns-resize', vertical: true });
      }
      expect(probes.length).toBeGreaterThan(0);
      for (const { part, cursor, vertical } of probes) {
        // a drag resizes the window: each probe starts from the same one
        v.view.setRange(S, 20, 60);
        v.view.setBand(RowsTrack, { bandTop: 10 });
        v.draw();
        nav.draw();
        nav.leave();
        const w = nav.draw()!.windows[0];
        const cx = (w.dx0 + w.dx1) / 2;
        const cy = (w.dy0 + w.dy1) / 2;
        const at: readonly [number, number] = vertical
          ? [cx, find(part, nav.height, (t) => nav.hit(cx, t)?.part)]
          : [find(part, nav.width, (t) => nav.hit(t, cy)?.part), cy];
        nav.move(at[0], at[1]);
        expect(nav.cursor).toBe(cursor);
        expect(nav.down(at[0], at[1])).toBe(true);
        nav.move(at[0] + 300, at[1] + 300);
        expect(nav.cursor).toBe(cursor);
        nav.up();
        nav.leave();
        expect(nav.cursor).toBe('default');
      }
      v.view.setRange(S, 20, 60);
      v.view.setBand(RowsTrack, { bandTop: 10 });
      v.draw();
      nav.draw();
      const g = nav.draw()!.windows[0];
      const bx = (g.dx0 + g.dx1) / 2;
      const by = (g.dy0 + g.dy1) / 2;
      expect(nav.hit(bx, by)?.part).toBe('body');
      nav.move(bx, by);
      expect(nav.cursor).toBe('grab');
      nav.down(bx, by);
      nav.move(bx + 5, by + 5);
      expect(nav.cursor).toBe('grabbing');
      nav.up();
      nav.leave();
      const o = horizontal ? ([3, (g.dy0 + g.dy1) / 2] as const) : ([3, nav.height - 3] as const);
      expect(nav.hit(o[0], o[1])?.part).toBe('outside');
      nav.move(o[0], o[1]);
      expect(nav.cursor).toBe('pointer');
    });
  }

  it('the ruler rail resizes along x, the rows rail along y; each grabs by its body and points outside', async () => {
    const v = await make(60);
    v.view.setRange(S, 20, 60);
    v.view.setBand(RowsTrack, { bandTop: 10 });
    v.draw();
    const sec = v.frame!.layout.sections[0];
    const rail = sec.rail!;
    const y = rail.bottom - 2;
    const partX = (x: number) => (v.input.move(x, y), v.input.windowHit?.part);
    const start = rail.x + find('start', rail.width, (t) => partX(rail.x + t));
    const end = rail.x + find('end', rail.width, (t) => partX(rail.x + t));
    const body = rail.x + find('body', rail.width, (t) => partX(rail.x + t));
    for (const [x, c] of [
      [start, 'ew-resize'],
      [end, 'ew-resize'],
      [body, 'grab'],
    ] as const) {
      v.input.move(x, y);
      expect(v.input.cursor).toBe(c);
    }
    v.input.down(start, y);
    v.input.move(start + 40, y + 200);
    expect(v.input.cursor).toBe('ew-resize');
    v.input.up();
    v.view.setRange(S, 20, 60);
    v.draw();
    v.input.move(rail.x + 1, y);
    expect(v.input.cursor).toBe('pointer');
    // the rows rail, upright at the band's right edge
    const rr = sec.rowsRails![0];
    const rx = rr.top + 6;
    const partY = (t: number) => (v.input.move(rx, t), v.input.windowHit?.part);
    const rs = rr.x + find('start', rr.width, (t) => partY(rr.x + t));
    const rb = rr.x + find('body', rr.width, (t) => partY(rr.x + t));
    v.input.move(rx, rb);
    expect(v.input.cursor).toBe('grab');
    v.input.move(rx, rs);
    expect(v.input.cursor).toBe('ns-resize');
    v.input.down(rx, rs);
    v.input.move(rx + 200, rs + 5);
    expect(v.input.cursor).toBe('ns-resize');
    v.input.up();
  });
});
