/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Track heights and labels beside at their edges (_spec/README.md "Render", Track heights, Labels inside if they fit,
// else beside): the label lane decided from the section's own extent (a pan over a cyclic origin, a follower zoomed),
// one lane block and one label lane a track, lanes never below the band, a tiled source's lanes never moving back on a
// pan, the sizing kept while its key is, a web font arriving, a spliced item on a ring, hidden layers, the
// navigator's pictures, a leader from its item's lane.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  type Frame,
  type FrameInput,
  type RangesLay,
  Renderer,
  type TrackLayout,
  createTextMetrics,
  recordDrawLists,
} from '../index.js';
import { describeItems } from '../items.js';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { settle } from './helpers.js';

type Rows = { start: number[]; end: number[]; id: string[]; label?: string[] };
const Size = { width: 800, height: 600 };

const lay = (f: Frame, track = 0, layer = 0, section = 0) =>
  f.layout.sections[section].tracks[track].layers[layer].lay as RangesLay;
const track = (f: Frame, t = 0, section = 0): TrackLayout => f.layout.sections[section].tracks[t];

describe('the label lane decided from the section’s own extent', () => {
  // 100 positions over the data box (~700 px): a 20-position item is ~140 px, its label fits; at half the scale not
  const cyclic = () => {
    const b = createSeqSpec({ title: 'cyclic' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true });
    const r = b.data.ranges('r', {
      target: s,
      rows: { start: [40], end: [60], id: ['a'], label: ['A fairly long label'] },
    });
    b.views.view('main', (v) => v.section('s', { of: s }).track('t').layer('l', { type: 'ranges', data: r }));
    return b.build();
  };

  it('a cyclic linear track panned across its origin keeps its height and its lanes', async () => {
    const s = await settle(cyclic());
    const at = async (viewport: { start: number; end: number }) => {
      const f = await s.frame({ size: Size, viewport });
      return [track(f).box.height, lay(f).laneY, lay(f).labelY];
    };
    const fit = await at({ start: 0, end: 100 });
    expect(await at({ start: 70, end: 130 })).toEqual(fit);
    expect(await at({ start: 0, end: 100 })).toEqual(fit);
  });

  it('a follower zoomed by its leader keeps its height and its lanes', async () => {
    const b = createSeqSpec({ title: 'follow' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
    // a 5-position item: ~35 px at the fit, its label beside; zoomed in, it would fit
    const r = b.data.ranges('r', { target: s, rows: { start: [40], end: [45], id: ['a'], label: ['Label beside'] } });
    b.views.view('main', (v) => {
      const lead = v.section('lead', { of: s });
      lead.track('t').layer('l', { type: 'letters' });
      v.section('follow', { of: s, follow: lead }).track('t').layer('l', { type: 'ranges', data: r });
    });
    const st = await settle(b.build());
    const at = async (viewport: { start: number; end: number }) => {
      const f = await st.frame({ size: Size, viewports: { '/views/main/lead': viewport } });
      return [track(f, 0, 1).box.height, lay(f, 0, 0, 1).laneY, lay(f, 0, 0, 1).labelY];
    };
    const fit = await at({ start: 0, end: 100 });
    expect(fit[2]).toBeGreaterThanOrEqual(0);
    expect(await at({ start: 35, end: 50 })).toEqual(fit);
  });
});

/** One track of two ranges layers, each its rows. */
function twoLayers(a: Rows, b: Rows): Node {
  const s = createSeqSpec({ title: 'two layers' });
  const seq = s.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
  const ra = s.data.ranges('a', { target: seq, rows: a });
  const rb = s.data.ranges('b', { target: seq, rows: b });
  s.views.view('main', (v) => {
    const t = v.section('s', { of: seq }).track('t');
    t.layer('a', { type: 'ranges', data: ra });
    t.layer('b', { type: 'ranges', data: rb });
  });
  return s.build();
}

describe('one lane block and one label lane a track', () => {
  it('two labeled layers in one track: their labels never overlap, their lanes one block', async () => {
    const s = await settle(
      twoLayers(
        { start: [50], end: [52], id: ['a'], label: ['First label'] },
        { start: [50], end: [52], id: ['b'], label: ['Second label'] },
      ),
    );
    const f = await s.frame({ size: Size, settled: true });
    const [a, b] = [lay(f, 0, 0), lay(f, 0, 1)];
    expect(a.labelY).toBeGreaterThanOrEqual(0);
    expect(b.labelY).toBe(a.labelY);
    expect(b.laneY).toBe(a.laneY);
    const boxes: [number, number][] = [];
    for (const l of [a, b]) for (let k = 0; k < l.ln; k++) boxes.push([l.lx[k], l.lx[k] + l.lw[k]]);
    // the first layer's label first (the layer's order), at its item's start
    expect(a.ln).toBe(1);
    boxes.sort((p, q) => p[0] - q[0]);
    for (let k = 1; k < boxes.length; k++) expect(boxes[k][0]).toBeGreaterThanOrEqual(boxes[k - 1][1]);
  });
});

/** One section of ranges tracks under `heights`. */
function tracks(rows: Rows[], heights?: string, height?: number): Node {
  const b = createSeqSpec({ title: 'tracks' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(50) });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s, ...(heights ? { heights } : {}) });
    rows.forEach((r, k) => {
      const d = b.data.ranges(`r${k}`, { target: s, rows: r });
      sec.track(`t${k}`, height !== undefined ? { height } : {}).layer('l', { type: 'ranges', data: d });
    });
  });
  return b.build();
}

describe('lanes never below the band', () => {
  it('a band taller than its content (uniform): a lane the sticky pass kept past the content is not drawn there', async () => {
    const one: Rows = { start: [10], end: [20], id: ['x'] };
    const three: Rows = { start: [0, 1, 2], end: [100, 100, 100], id: ['p', 'q', 'r'] };
    const s = await settle(tracks([one, three], 'uniform'));
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 200 } });
    const l = lay(f, 0);
    // the item's lane kept at 2, as a layout of another view could have left it
    l.sticky.assign(['k0', 'k1'], [0, 0], [50, 50], [0, 1], 5, false, new Int32Array(2));
    l.sticky.assign(['k0', 'k1', 'sx#0'], [0, 0, 0], [50, 50, 50], [0, 1, 2], 5, false, new Int32Array(3));
    expect(l.sticky.laneOf('sx#0')).toBe(2);
    const g = await s.frame({ size: Size, viewport: { start: 1, end: 199 } });
    const m = lay(g, 0);
    const box = track(g, 0).box;
    for (let j = 0; j < m.m; j++)
      expect(m.laneY + (m.mlane[j] + 1) * m.laneHeight).toBeLessThanOrEqual(box.y + box.height + 1e-9);
  });

  it('an exact fit: a sticky lane that would overflow is packed fresh, never sent to the density row', async () => {
    // A [0,50), B [40,100), C [60,90): two lanes over the whole data (A0, B1, C0)
    const rows: Rows = { start: [0, 40, 60], end: [50, 100, 90], id: ['a', 'b', 'c'] };
    const s = await settle(tracks([rows]));
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 200 } });
    expect(track(f).box.height).toBe(36);
    // C kept at lane 1, as a view with C and something under it would have left it
    lay(f).sticky.assign(['k0'], [60], [90], [0], 5, false, new Int32Array(1));
    lay(f).sticky.assign(['k0', 'sc#0'], [60, 60], [90, 90], [0, 1], 5, false, new Int32Array(2));
    expect(lay(f).sticky.laneOf('sc#0')).toBe(1);
    const g = await s.frame({ size: Size, viewport: { start: 1, end: 199 } });
    expect(lay(g).overflow).toBe(0);
    expect(lay(g).mode).toBe('items');
    expect(lay(g).density).toBeUndefined();
  });
});

describe('a source’s lanes', () => {
  it('never move back on a pan: centred on the most lanes the layer has shown', async () => {
    // a file's rows (not known whole before they are read): one item
    const gff = ['##gff-version 3', 's\tt\tCDS\t11\t20\t.\t+\t.\tID=one'].join('\n');
    const b = createSeqSpec({ title: 'source' });
    const seq = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
    const file = b.data.source('file', { format: 'gff3', value: gff });
    const f = b.data.ranges('f', { target: seq, source: file });
    b.views.view('main', (v) =>
      v.section('s', { of: seq }).track('t', { height: 64 }).layer('l', { type: 'ranges', data: f }),
    );
    const s = await settle(b.build(), [RecordsRuntime]);
    const at = async (start: number, settled = false) => {
      const fr = await s.frame({ size: Size, settled, viewport: { start, end: start + 30 } });
      return { y: lay(fr).laneY, n: lay(fr).laneCount, lay: lay(fr) };
    };
    const first = await at(0);
    expect(first.n).toBe(1);
    // the alone item kept in a third lane, as a view before could have left it: three lanes shown
    first.lay.sticky.assign(['k0', 'k1'], [0, 0], [50, 50], [0, 1], 5, false, new Int32Array(2));
    first.lay.sticky.assign(['k0', 'k1', 'sone#0'], [0, 0, 0], [50, 50, 50], [0, 1, 2], 5, false, new Int32Array(3));
    const busy = await at(1);
    expect(busy.n).toBe(3);
    expect(busy.y).toBeLessThan(first.y);
    // repacked when the view settles (one lane again): the lanes' top stays where three stood
    const back = await at(0, true);
    expect(back.n).toBe(1);
    expect(back.y).toBe(busy.y);
  });
});

describe('the sizing kept while its key is', () => {
  it('the host’s port runs once for a section until its tracks, its width, the config or the metrics change', async () => {
    let calls = 0;
    const s = await settle(
      tracks([{ start: [10], end: [20], id: ['x'] }]),
      [],
      {},
      {
        trackSizing: (_s, _t, named) => {
          calls++;
          return named;
        },
      },
    );
    for (let k = 0; k < 6; k++) await s.frame({ size: Size, viewport: { start: k * 10, end: k * 10 + 100 } });
    expect(calls).toBe(1);
    await s.frame({ size: { width: 700, height: 600 } });
    expect(calls).toBe(2);
  });
});

describe('a web font arriving', () => {
  it('the track labels’ ellipsis and the label lane are found again', async () => {
    let scale = 1;
    const text = createTextMetrics(() => (t) => t.length * 7 * scale);
    const b = createSeqSpec({ title: 'font' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
    const r = b.data.ranges('r', { target: s, rows: { start: [10], end: [30], id: ['a'], label: ['Label inside'] } });
    b.views.view('main', (v) =>
      v.section('s', { of: s }).track('t', { label: 'A track label longer than the column' }).layer('l', {
        type: 'ranges',
        data: r,
      }),
    );
    const st = await settle(b.build());
    const renderer = new Renderer({ ...st.renderer.options, text });
    const frame = (input: FrameInput) => renderer.frame(input);
    const before = frame({ size: Size, settled: true });
    // the column at its most (160 px) before and after: only the metrics change
    const cut = track(before).labelText;
    expect(cut.endsWith('…')).toBe(true);
    expect(lay(before).labelY).toBe(-1);
    // the font arrives wider: the label no longer fits the column, nor the item's mark
    scale = 2;
    text.clear();
    const after = frame({ size: Size, settled: true });
    expect(track(after).labelText.endsWith('…')).toBe(true);
    expect(track(after).labelText.length).toBeLessThan(cut.length);
    expect(lay(after).labelY).toBeGreaterThanOrEqual(0);
  });
});

describe('a spliced item on a ring', () => {
  it('two rows whose gap is inside the circle is no item across the origin: no runs', async () => {
    const b = createSeqSpec({ title: 'spliced' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(250), cyclic: true });
    const r = b.data.ranges('r', { target: s, rows: { start: [100, 700], end: [200, 800], id: ['x', 'x'] } });
    b.views.view('main', (v) => v.section('s', { of: s }).track('t').layer('l', { type: 'ranges', data: r }));
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, settled: true });
    const [item] = describeItems(lay(f));
    expect(item.runs).toBeUndefined();
    expect([item.start, item.end]).toEqual([100, 800]);
  });
});

describe('hidden layers', () => {
  it('a hidden layer neither names the track’s kind nor adds its fixed height', async () => {
    const b = createSeqSpec({ title: 'hidden' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
    const r = b.data.ranges('r', { target: s, rows: { start: [10], end: [20], id: ['x'] } });
    const v = b.data.values('v', { target: s, values: new Float32Array(100) });
    b.views.view('main', (view) => {
      const sec = view.section('s', { of: s });
      const t = sec.track('a');
      t.layer('r', { type: 'ranges', data: r });
      t.layer('v', { type: 'values', data: v, hidden: true });
      sec.track('b').layer('r', {
        type: 'ranges',
        data: b.data.ranges('r2', { target: s, rows: { start: [0, 1, 2], end: [90, 90, 90], id: ['p', 'q', 'r'] } }),
      });
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, settled: true });
    // a ranges track as its neighbour: one kind, the tallest (three lanes), not the values' 32 nor its own 22
    expect(track(f, 0).box.height).toBe(50);
    expect(track(f, 1).box.height).toBe(50);
  });
});

describe('a leader', () => {
  it('starts at its item’s mark, moved clear of the label before it, crossing nothing', async () => {
    // C's label ends just right of B's start: B's label moves a few px right of it, its leader from B's mark
    const rows: Rows = {
      start: [10, 14],
      end: [12, 15],
      id: ['c', 'b'],
      label: ['Cc', 'Bb'],
    };
    const s = await settle(tracks([rows], undefined));
    const f = await s.frame({ size: Size, settled: true });
    const l = lay(f);
    const b = Array.from({ length: l.items }, (_, i) => i).find((i) => l.keys[i] === 'b')!;
    const lane = l.lanes[b];
    const paths = recordDrawLists(f.lists, ['labels'])
      .split('\n')
      .filter((x) => x.startsWith('path '));
    expect(paths.length).toBe(1);
    const y0 = Number(paths[0].replace(/[ML]/g, ' ').split(/\s+/).filter(Boolean)[2]);
    expect(y0).toBeCloseTo(l.laneY + lane * l.laneHeight + l.gap / 2, 6);
  });

  it('a label under another lane’s mark has no leader up through it: dropped', async () => {
    // A in lane 0 over B's start, C before B: B (lane 1) would need a leader through A's mark
    const rows: Rows = {
      start: [0, 10, 14],
      end: [100, 12, 15],
      id: ['a', 'c', 'b'],
      label: ['Long domain', 'Cc', 'Bb'],
    };
    const s = await settle(tracks([rows], undefined));
    const f = await s.frame({ size: Size, settled: true });
    const l = lay(f);
    expect(l.lanes[Array.from({ length: l.items }, (_, i) => i).find((i) => l.keys[i] === 'b')!]).toBeGreaterThan(0);
    const placed = Array.from({ length: l.ln }, (_, k) => l.keys[l.litem[k]]);
    expect(placed).toEqual(['c']);
  });
});
