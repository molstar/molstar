/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows' zoom limits (_spec/README.md "Configuration", Rows; core's _spec Rows, The rows rail), headless on a fake
// clock: a lazy alignment serving no row levels stops its zoom-out at a device pixel a row (the rail's grip and `[`), so
// the window never reaches the blank state; `rows.zoom` makes the bounds the host's (`minPx`, `maxPx`), and `enabled:
// false` keeps the rail a pan-only bar (no grips, the keys inert, the scroll kept).

import type { Node } from '@molstar/core/repr/tree';
import { cellsSource, type FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { type Rail, type RowsLay } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Track = '/views/v/s/rows';
const Height = 200;

/** A lazy alignment of `n` rows: its rows read by index (it serves no row levels). */
const Lazy: FormatRuntime = {
  name: 'lazy-rows',
  formats: ['lazy-rows'],
  async open(source) {
    const n = (source.value as { rows: number }).rows;
    return cellsSource(n, 60, (start, end) => ({
      codes: new Uint8Array((end - start) * 60).fill(65),
      label: Array.from({ length: end - start }, (_, k) => `L${start + k}`),
      first: new Int32Array(end - start).fill(1),
    }));
  },
};

function lazySpec(n: number): Node {
  const b = createSeqSpec({ title: 'lazy' });
  const src = b.data.source('src', { format: 'lazy-rows', value: { rows: n } });
  const aln = b.data.alignment('aln', { alphabet: 'protein', source: src });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }).track('rows', { label: '', height: Height }).layer('r', { type: 'rows', data: aln }),
  );
  return b.build();
}

/** A held alignment of `n` rows. */
function heldSpec(n: number): Node {
  const b = createSeqSpec({ title: 'held' });
  const label = Array.from({ length: n }, (_, i) => `r${i}`);
  const letters = Array.from({ length: n }, () => 'MEEPQSDPSV'.repeat(4));
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label, letters } });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }).track('rows', { label: '', height: Height }).layer('r', { type: 'rows', data: aln }),
  );
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make(root: Node, config: readonly (readonly [unknown, unknown])[] = []) {
  frames = installFrames();
  const base = DefaultSeqViewerSpec();
  const viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), { kind: 'format-runtime', runtime: Lazy }],
    config: [[viewerItem('view.jumpMs'), 0], ...config],
  } as never);
  viewers.push(viewer);
  await viewer.load(root);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.input.mac = false;
  viewer.draw();
  return viewer;
}

const band = (v: SeqViewer) => v.view.band(Track);
const rail = (v: SeqViewer): Rail => v.frame!.layout.sections[0].rowsRails![0];
const rowsLay = (v: SeqViewer) =>
  v.frame!.layout.sections[0].tracks.flatMap((t) => t.layers).find((l) => l.type.name === 'rows')!.lay as RowsLay;
const yOf = (r: Rail, u: number) => r.x + (u / r.total) * r.width;

/** The y of the rail's end grip (found by hit: a window of few rows of many is drawn floored). */
function endGrip(v: SeqViewer): [number, number] {
  const r = rail(v);
  const x = r.top + 6;
  for (let y = r.x; y < r.x + r.width; y += 0.5) if (v.hit(x, y)?.part === 'end') return [x, y];
  throw new Error('no end grip');
}

/** The end grip dragged to the rail's end: the window as wide as the zoom lets it. */
function dragEndToTheFoot(v: SeqViewer) {
  const r = rail(v);
  const [x, y] = endGrip(v);
  v.input.down(x, y);
  v.input.move(x, r.x + r.width);
  v.input.up();
  v.draw();
}

describe('the rows’ zoom floor on a lazy alignment', () => {
  it('the rail’s grip stops at a device pixel a row: the window never reaches the blank state', async () => {
    const v = await make(lazySpec(1_000_000));
    expect(rowsLay(v).rowFloor).toBe(1 / v.dpr);
    expect(rail(v).rowMin).toBe(1 / v.dpr);
    dragEndToTheFoot(v);
    // the band's px over the rows in view: a device pixel a row at the least
    expect(band(v)?.rowScale).toBeGreaterThanOrEqual(1 / v.dpr - 1e-9);
    expect(band(v)?.rowScale).toBeLessThan(2 / v.dpr);
    const lay = rowsLay(v);
    expect(lay.neutral).toBe(false);
    expect(v.frame!.layout.issues.some((i) => i.code === 'seq.render.rows.reduced')).toBe(false);
  });

  it('`[` stops there too, and says it is so', async () => {
    const v = await make(lazySpec(1_000_000));
    for (let k = 0; k < 12; k++) {
      v.input.key({ key: '[' });
      frames!.flush();
      v.draw();
    }
    expect(band(v)?.rowScale).toBeGreaterThanOrEqual(1 / v.dpr - 1e-9);
    expect(v.input.key({ key: '[' })).toBe(true);
    expect(band(v)?.rowScale).toBeCloseTo(1 / v.dpr, 9);
    expect(rowsLay(v).neutral).toBe(false);
  });

  it('a held alignment has no source floor: its zoom goes down to the configured least', async () => {
    const v = await make(heldSpec(2000));
    expect(rowsLay(v).rowFloor).toBe(0);
    dragEndToTheFoot(v);
    expect(band(v)?.rowScale).toBeLessThan(1 / v.dpr);
  });
});

describe('rows.zoom: the limits the host’s', () => {
  it('`minPx` and `maxPx` bound the keys and the grips; the floor the larger of `minPx` and the source’s', async () => {
    const v = await make(heldSpec(2000), [[viewerItem('rows.zoom'), { minPx: 4, maxPx: 8 }]]);
    expect(rail(v).rowMin).toBe(4);
    dragEndToTheFoot(v);
    expect(band(v)?.rowScale).toBeCloseTo(4, 6);
    for (let k = 0; k < 6; k++) {
      v.input.key({ key: ']' });
      frames!.flush();
      v.draw();
    }
    expect(band(v)?.rowScale).toBeCloseTo(8, 6);
    // a lazy source's floor over a smaller `minPx`
    const w = await make(lazySpec(1_000_000), [[viewerItem('rows.zoom'), { minPx: 0.25 }]]);
    expect(rail(w).rowMin).toBe(1 / w.dpr);
  });

  it('`enabled: false`: the rail is a pan-only bar, the zoom keys are inert, the scroll is kept', async () => {
    const v = await make(heldSpec(200), [[viewerItem('rows.zoom'), { enabled: false }]]);
    const r = rail(v);
    expect(r.pan).toBe(true);
    expect(r.active).toBe(true);
    // a press on a grip's place is the body's: the window scrolls, the rows keep their size
    const x = r.top + 6;
    for (let y = r.x; y < r.x + r.width; y += 0.5) {
      const part = v.hit(x, y)?.part;
      expect(part === 'end' || part === 'start').toBe(false);
    }
    expect(v.input.key({ key: '[' })).toBe(false);
    expect(v.input.key({ key: ']' })).toBe(false);
    expect(band(v)?.rowScale).toBeUndefined();
    // the window's body scrolls the rows
    const mid = yOf(r, (r.start + r.end) / 2);
    v.input.down(x, mid);
    v.input.move(x, mid + 20);
    v.input.up();
    expect(band(v)?.bandTop ?? 0).toBeGreaterThan(0);
    expect(band(v)?.rowScale).toBeUndefined();
  });

  it('on by default: the same rail has its grips and the keys zoom', async () => {
    const v = await make(heldSpec(200));
    expect(rail(v).pan).toBe(false);
    expect(v.input.key({ key: '[' })).toBe(true);
  });
});

describe('a stored rowScale and a bad rows.zoom', () => {
  it('a saved scale under the floor is drawn at it (dpr 2), the rail’s bounds agreeing with the layout', async () => {
    const v = await make(lazySpec(1_000_000), [[viewerItem('canvas.dpr'), 2]]);
    v.view.setBand(Track, { rowScale: 0.1, bandTop: 0 });
    v.draw();
    const lay = rowsLay(v);
    expect(v.dpr).toBe(2);
    expect(lay.rowScale).toBe(0.5);
    expect(lay.neutral).toBe(false);
    expect(rail(v).rowMin).toBe(lay.rowScale);
    // and above the host's most
    v.view.setBand(Track, { rowScale: 500 });
    v.draw();
    expect(rowsLay(v).rowScale).toBe(40);
  });

  it('with the zoom off a stored scale is ignored: the fit’s', async () => {
    const v = await make(heldSpec(200), [[viewerItem('rows.zoom'), { enabled: false }]]);
    const fit = rowsLay(v).rowScale;
    v.view.setBand(Track, { rowScale: 3 });
    v.draw();
    expect(rowsLay(v).rowScale).toBe(fit);
  });

  it('a bad rows.zoom is the default; a double-click resets nothing with the zoom off', async () => {
    for (const bad of [{ minPx: 0 }, { minPx: 9, maxPx: 3 }, { maxPx: Number.NaN }]) {
      const v = await make(heldSpec(200), [[viewerItem('rows.zoom'), bad]]);
      expect(v.config.rows.zoom).toEqual({ minPx: 1 / 64, maxPx: 40, enabled: true });
    }
    const v = await make(heldSpec(200), [[viewerItem('rows.zoom'), { enabled: false }]]);
    v.view.setBand(Track, { rowScale: 3 });
    v.draw();
    const r = rail(v);
    v.input.doubleClick(r.top + 6, yOf(r, (r.start + r.end) / 2));
    expect(band(v)?.rowScale).toBe(3);
  });
});
