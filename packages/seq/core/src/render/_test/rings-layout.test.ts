/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rings layout: placement as the inverse of hit-testing over random rotations (every position placed once, `at`
// the inverse of `place`, `offsetAt` monotone in angle), a linear detail following the ring's window (over the origin
// too) and a ring following a detail, the window out of the ring's key (a pan of the detail draws the ring from its
// kept lays, its labels not placed again), drawing against the data, labels outside hit as their items, the window's
// parts, reused buffers against fresh ones, Canvas2D-free lists equal in SVG.

import { describe, expect, it } from 'vitest';
import {
  type Frame,
  type Hit,
  PlacePoint,
  type RingPlacement,
  RingsNav,
  type ValuesLay,
  type NavContext,
  NoColor,
  type RangesLay,
  recordDrawLists,
  drawListsToSvg,
  hitView,
  openGroupOf,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { FakeClock } from '../../data/_test/synthetic.js';
import { rangesLocation } from '@molstar/seqx-core/query';
import { settle } from './helpers.js';
import { plasmidRingSpec, RingFeatures, RingLength, ringLetters, ringSpec } from './rings-fixtures.js';
import { BigRuntime, bigValue } from './fixtures.js';
import { layerLay, type ReadPath, readPaths, rng } from './small-helpers.js';

const Size = { width: 640, height: 900 };

const ringOf = (f: Frame) => f.layout.sections.find((s) => s.layout.name === 'rings')!;

describe('the rings layout', () => {
  it('places the section round its whole axis, the first track outermost, the hole kept, the title in it', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 200 } } });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    expect(p.strips.n).toBe(1);
    expect(ring.axis!.u0).toBe(0);
    expect(ring.axis!.u1).toBe(RingLength);
    // the first band from the outer radius in; the rings' depth keeps a hole of at least 0.35 of the radius
    expect(p.strips.bandY[0]).toBe(0);
    expect(p.depth).toBeLessThanOrEqual(p.r * 0.65 + 1e-9);
    expect(p.title).toBe('construct');
    expect(p.subtitle).toBe('600 bp');
    expect(f.layout.issues.filter((i) => i.severity === 'error')).toEqual([]);
  });

  it('a section’s `center: none` leaves the hole empty, over the render config’s `rings.center`', async () => {
    const none = ringOf(await (await settle(ringSpec({ center: 'none' }))).frame({ size: Size }));
    expect([(none.placement as RingPlacement).title, (none.placement as RingPlacement).subtitle]).toEqual(['', '']);
    // the param wins over the preset either way
    const title = await settle(ringSpec({ center: 'title' }), [], { rings: { center: 'none' } });
    expect((ringOf(await title.frame({ size: Size })).placement as RingPlacement).title).toBe('construct');
    const preset = await settle(ringSpec(), [], { rings: { center: 'none' } });
    expect((ringOf(await preset.frame({ size: Size })).placement as RingPlacement).title).toBe('');
  });

  it('placement is the inverse of hit-testing over rotations: every position placed once, `at` gives it back, `offsetAt` monotone in angle', async () => {
    const r = rng(3);
    for (let n = 0; n < 4; n++) {
      const rotation = Math.floor(r() * 360);
      const s = await settle(ringSpec({ rotation, detail: false }));
      const f = await s.frame({ size: Size });
      const p = ringOf(f).placement as RingPlacement;
      const pt = new PlacePoint();
      const seen = new Set<number>();
      let last = -Infinity;
      for (let u = 0; u < RingLength; u++) {
        expect(p.place(u + 0.5, pt)).toBe(true);
        const { x, y } = pt;
        expect(p.at([x, y], pt)).toBe(true);
        expect(pt.kind).toBe('data');
        const back = Math.floor(pt.offset);
        expect(back).toBe(u);
        seen.add(back);
        const o = p.offsetAt(x, y);
        expect(o).toBeGreaterThan(last);
        last = o;
      }
      expect(seen.size).toBe(RingLength);
    }
  });

  it('a linear detail follows the ring’s window, over the origin as two pieces; the ring follows a detail the same', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 580, end: 620 } } });
    const detail = f.layout.sections.find((x) => x.section.address === '/views/v/detail')!;
    expect(detail.axis!.pieces).toEqual([580, 600, 0, 20]);
    const p = ringOf(f).placement as RingPlacement;
    expect([p.windowShown, p.windowActive, p.windowStart, p.windowEnd]).toEqual([true, true, 580, 620]);
    // the ring follows: the detail's range is its window
    const t = await settle(ringSpec({ ringFollows: true }));
    const g = await t.frame({ size: Size, viewports: { '/views/v/detail': { start: 100, end: 160 } } });
    const q = ringOf(g).placement as RingPlacement;
    expect(ringOf(g).axis!.total).toBe(RingLength);
    expect([q.windowShown, q.windowStart, q.windowEnd]).toEqual([true, 100, 160]);
  });

  it('a pan of the detail does not lay the ring out again: its lays and labels kept, the window moved', async () => {
    const s = await settle(ringSpec());
    const f1 = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 200 } } });
    const ring1 = ringOf(f1);
    const lay = ring1.lays[0][0].lay;
    const labels = (ring1.placement as RingPlacement).labels.m;
    const f2 = s.renderer.frame({ size: Size, viewports: { '/views/v/ring': { start: 140, end: 240 } } });
    expect(f2.rebuilt).toBe(true);
    const ring2 = ringOf(f2);
    expect(ring2).toBe(ring1);
    expect(ring2.lays[0][0].lay).toBe(lay);
    expect((ring2.placement as RingPlacement).labels.m).toBe(labels);
    expect((ring2.placement as RingPlacement).windowStart).toBe(140);
    // drawn as a fresh layout draws it
    const fresh = await settle(ringSpec());
    const g = await fresh.frame({ size: Size, viewports: { '/views/v/ring': { start: 140, end: 240 } } });
    expect(recordDrawLists(f2.lists)).toBe(recordDrawLists(g.lists));
  });

  it('labels outside, each hit as its item; the feature over the origin one label, at its middle angle', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const f = await s.frame({ size: Size });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    const texts = Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]).sort();
    expect(texts).toEqual(['AmpR', 'MCS', 'ori']);
    for (let k = 0; k < p.labels.m; k++) {
      const i = p.labels.placed[k];
      const x = p.labels.right[k] ? p.labels.x[k] + 2 : p.labels.x[k] - 2;
      const h = hitAtPage(f, x, p.labels.y[k] - 3);
      expect(h?.element?.key).toBe(p.labels.hits[i]?.element?.key);
    }
    // ori runs 560 to 40 over the origin: its label's angle is its middle, position 0's
    const ori = p.labels.text.indexOf('ori');
    expect(p.labels.angle[ori]).toBeCloseTo(0, 6);
    // no label is drawn inside a mark on a ring
    expect(f.lists.passes.labels.texts.n).toBeGreaterThanOrEqual(3);
  });

  it('draws the features against the data: each feature’s arc where its positions are, over the origin one arc', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const f = await s.frame({ size: Size });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    // across the features ring at each position: the features holding it, and nothing where none is
    const keysAt = (pos: number) => {
      const keys = new Set<string | number>();
      const x = p.frameX(pos + 0.5);
      for (let y = 1; y < p.strips.bandH[0]; y += 2) {
        p.warp.map(x, y);
        const k = hitAtPage(f, p.warp.px, p.warp.py)?.element?.key;
        if (k !== undefined) keys.add(k);
      }
      return [...keys].sort();
    };
    const rows = RingFeatures;
    for (let pos = 0; pos < RingLength; pos += 7) {
      const want = rows.id.filter((_, i) =>
        rows.start[i] < rows.end[i]
          ? pos >= rows.start[i] && pos < rows.end[i]
          : pos >= rows.start[i] || pos < rows.end[i],
      );
      expect(keysAt(pos)).toEqual(want.sort());
    }
    // the letters the detail would show are the data's
    expect(ringLetters().length).toBe(RingLength);
  });

  it('the window’s parts on the rim: body and ends hit as the window target; idle with no follower: no window', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 200 } } });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    const at = (frameX: number) => {
      p.warp.map(frameX, -5);
      return hitAtPage(f, p.warp.px, p.warp.py);
    };
    expect(at(p.frameX(150))).toMatchObject({ target: 'window', part: 'body', section: 'ring' });
    expect(at(p.frameX(100))).toMatchObject({ target: 'window', part: 'start' });
    expect(at(p.frameX(200))).toMatchObject({ target: 'window', part: 'end' });
    // a window drawn in the window pass, nothing of it in an export by default
    expect(f.lists.passes.window.opCount).toBeGreaterThan(0);
    const alone = await settle(ringSpec({ detail: false }));
    const g = await alone.frame({ size: Size });
    expect((ringOf(g).placement as RingPlacement).windowShown).toBe(false);
    expect(g.lists.passes.window.opCount).toBe(0);
  });

  it('a window with nowhere to go (the whole ring) is quiet: no wash, no grips, no hover', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 0, end: RingLength } } });
    const p = ringOf(f).placement as RingPlacement;
    expect(p.windowShown).toBe(true);
    expect(p.windowActive).toBe(false);
    const paths = readPaths(f.lists, 'window');
    // its outline alone: nothing filled
    expect(paths.length).toBe(1);
    expect(paths[0].fill).toBe(0);
    p.warp.map(p.frameX(300), -5);
    expect(hitAtPage(f, p.warp.px, p.warp.py)).toBeUndefined();
  });

  it('the window’s grips hit across the origin: a start just past it, an end just before it', async () => {
    const s = await settle(ringSpec());
    for (const [start, end, edge, part] of [
      [1, 60, 1, 'start'],
      [520, 599, 599, 'end'],
    ] as const) {
      const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start, end } } });
      const p = ringOf(f).placement as RingPlacement;
      const x = p.frameX(edge);
      for (const dx of [-4, 4]) {
        p.warp.map(x + dx, -5);
        expect(hitAtPage(f, p.warp.px, p.warp.py)).toMatchObject({ target: 'window', part });
      }
    }
  });

  it('the window’s rim is the design’s 12 px hit; a narrow window’s sector grows about its middle, as its rim does', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 101 } } });
    const p = ringOf(f).placement as RingPlacement;
    const mid = p.frameX(100.5);
    p.warp.map(mid, -11.5);
    expect(hitAtPage(f, p.warp.px, p.warp.py)).toMatchObject({ target: 'window' });
    // the sector (the window pass's first path, inside the ring): its angular middle is the window's
    const angle = (x: number, y: number) => Math.atan2(x - p.cx, -(y - p.cy));
    const sector = readPaths(f.lists, 'window')[0];
    let lo = Infinity;
    let hi = -Infinity;
    for (let i = 0; i + 1 < sector.xy.length; i += 2) {
      const t = angle(sector.xy[i], sector.xy[i + 1]);
      lo = Math.min(lo, t);
      hi = Math.max(hi, t);
    }
    p.warp.map(mid, 1);
    expect((lo + hi) / 2).toBeCloseTo(angle(p.warp.px, p.warp.py), 4);
  });

  it('the ring’s lists in SVG: every primitive bent onto the page (no rect survives the warp)', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const f = await s.frame({ size: Size });
    const svg = drawListsToSvg(f.lists);
    expect(svg).toContain('<path');
    expect(f.lists.passes.data.rects.n).toBe(0);
  });
});

describe('a ring’s rotation and what it refuses', () => {
  it('a rotation turns the ring (laid out again) and leaves the window’s offsets, so the detail does not move', async () => {
    const viewports = { '/views/v/ring': { start: 100, end: 160 } };
    const a = await (await settle(ringSpec())).frame({ size: Size, viewports });
    const b = await (await settle(ringSpec({ rotation: 90 }))).frame({ size: Size, viewports });
    const pa = ringOf(a).placement as RingPlacement;
    const pb = ringOf(b).placement as RingPlacement;
    expect([pb.windowStart, pb.windowEnd]).toEqual([pa.windowStart, pa.windowEnd]);
    expect(b.layout.sections[1].axis!.pieces).toEqual(a.layout.sections[1].axis!.pieces);
    // position 0 at the top unturned, at the right turned a quarter clockwise
    const pt = new PlacePoint();
    pa.place(0, pt);
    expect(pt.x).toBeCloseTo(pa.cx, 6);
    expect(pt.y).toBeLessThan(pa.cy);
    pb.place(0, pt);
    expect(pt.x).toBeGreaterThan(pb.cx);
    expect(pt.y).toBeCloseTo(pb.cy, 6);
  });

  it('a ranges layer past its lanes draws its density as a ring band: its runs as sectors, never refused', async () => {
    const starts = Array.from({ length: 60 }, (_, i) => (i * 7) % 590);
    const b = createSeqSpec({ title: 'dense' });
    const c = b.data.sequence('c', { alphabet: 'dna', letters: ringLetters(), cyclic: true });
    const f = b.data.ranges('f', {
      target: c,
      rows: { start: starts, end: starts.map((s) => s + 200), id: starts.map((_, i) => `i${i}`) },
    });
    b.views.view('v', (v) =>
      v.section('ring', { of: c, layout: 'rings' }).track('f', { height: 12 }).layer('f', {
        type: 'ranges',
        data: f,
      }),
    );
    const s = await settle(b.build(), [], { ranges: { overflow: 'layer' } });
    const frame = await s.frame({ size: Size });
    expect(frame.layout.issues).toEqual([]);
    const lay = frame.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    expect(lay.mode).toBe('density');
    // the strip's image is not on the page: its runs of one color are, bent, a sector each, every one inside the ring
    const data = frame.lists.passes.data;
    expect(data.images.n).toBe(0);
    const p = ringOf(frame).placement as RingPlacement;
    const pts = data.paths.points;
    let sectors = 0;
    for (let k = 0; k < data.paths.n; k++) {
      if (data.paths.fill[k] === NoColor) continue;
      sectors++;
      const at = 2 * data.paths.pointStart[k];
      expect(Math.hypot(pts[at] - p.cx, pts[at + 1] - p.cy)).toBeLessThanOrEqual(p.r + 1);
    }
    expect(sectors).toBeGreaterThan(1);
  });
});

describe('a ring in two consumers and in an export', () => {
  it('renderSvg draws the lists a renderer draws (the window only with `window`)', async () => {
    const spec = ringSpec();
    const viewports = { '/views/v/ring': { start: 500, end: 640 } };
    const out = await renderSvg(spec, { size: Size, viewports });
    const s = await settle(spec);
    const f = await s.frame({ size: Size, viewports, figure: true, settled: true });
    expect(recordDrawLists(out.lists)).toBe(recordDrawLists(f.lists));
    expect(out.lists.passes.window.opCount).toBe(0);
    const windowed = await renderSvg(spec, { size: Size, viewports, window: true });
    expect(windowed.lists.passes.window.opCount).toBeGreaterThan(0);
    expect(windowed.svg).toContain('<path');
  });

  it('a selection over the origin is one arc on the ring: two runs whose square ends meet at the seam', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const selection = rangesLocation('/data/construct', [0, 30, 570, 600]);
    const f = await s.frame({ size: Size, emphasis: { selection } });
    const paths = readPaths(f.lists, 'selection');
    // a tint and an outline a run, two runs: each a bent path on the page
    expect(paths.length).toBe(4);
    const p = ringOf(f).placement as RingPlacement;
    // the runs meet at the origin's angle: a point of each at the seam
    p.warp.map(p.frameX(0), 1);
    const seam = { x: p.warp.px, y: p.warp.py };
    const near = (q: ReadPath) => {
      for (let i = 0; i + 1 < q.xy.length; i += 2)
        if (Math.hypot(q.xy[i] - seam.x, q.xy[i + 1] - seam.y) < 2) return true;
      return false;
    };
    expect(paths.filter(near).length).toBeGreaterThanOrEqual(2);
  });
});

describe('the rings layout’s navigation', () => {
  const c: NavContext = { total: 600, width: 400, height: 400, cyclic: true, maxScale: 20 };
  it('pans round the origin, zooms about a unit, reveals by the shortest way, fits', () => {
    const v = { start: 580, end: 620 };
    expect(RingsNav.pan(c, v, -200, 0)).toMatchObject({ start: 0, end: 40 });
    expect(RingsNav.reveal(c, v, 590, 610)).toBe(v);
    expect(RingsNav.reveal(c, v, 5, 15)).toBe(v);
    expect(RingsNav.reveal(c, v, 300, 310)).toMatchObject({ start: 285, end: 325 });
    expect(RingsNav.fit(c, v)).toMatchObject({ start: 0, end: 600 });
    expect(RingsNav.step(c, 599, 'right')).toBe(0);
  });
});

describe('a plasmid’s labels', () => {
  /** The item id a plasmid label names (`feature 3` → `f3`, `Site3` → `s3`; `ori`, `bla` their own). */
  const idOf = (text: string) =>
    text.startsWith('feature ') ? `f${text.slice(8)}` : text.startsWith('Site') ? `s${text.slice(4)}` : text;

  it('each placed label hits its own item (about 30 labels)', async () => {
    const s = await settle(plasmidRingSpec());
    const f = await s.frame({ size: Size });
    const p = ringOf(f).placement as RingPlacement;
    const l = p.labels;
    expect(l.m).toBeGreaterThanOrEqual(25);
    for (let k = 0; k < l.m; k++) {
      const x0 = l.right[k] ? l.x[k] : l.x[k] - l.width[k];
      const h = hitAtPage(f, x0 + l.width[k] / 2, l.y[k] - 0.35 * p.labelSize);
      expect(String(h?.element?.key)).toBe(idOf(l.text[l.placed[k]]));
    }
  });

  it('no two labels’ hit boxes overlap', async () => {
    const s = await settle(plasmidRingSpec());
    const f = await s.frame({ size: Size });
    const p = ringOf(f).placement as RingPlacement;
    const l = p.labels;
    let overlaps = 0;
    for (let y = Math.floor(p.box.y); y < p.box.y + p.box.height; y += 0.5) {
      for (let x = 0; x < Size.width; x += 2) {
        let n = 0;
        for (let k = 0; k < l.m; k++) if (p.inLabel(k, x, y)) n++;
        if (n > 1) overlaps++;
      }
    }
    expect(overlaps).toBe(0);
  });
});

describe('a dense fan of labels', () => {
  /** A 2,686 bp ring: `dense` sites 4 bp apart (a 2° fan at most), and three features across the ring. */
  function denseSpec(dense: number) {
    const n = 2686;
    const sites = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
    for (let k = 0; k < dense; k++) {
      sites.start.push(1081 + 4 * k);
      sites.end.push(1085 + 4 * k);
      sites.id.push(`s${k}`);
      sites.label.push(`Site${k}`);
    }
    const far = { start: [100, 700, 2000], end: [180, 780, 2080], id: ['a', 'b', 'c'], label: ['a', 'b', 'c'] };
    const b = createSeqSpec({ title: 'dense ring' });
    const c = b.data.sequence('plasmid', { alphabet: 'dna', letters: 'A'.repeat(n), cyclic: true });
    const fs = b.data.ranges('sites', { target: c, rows: sites });
    const ff = b.data.ranges('far', { target: c, rows: far });
    b.views.view('v', (v) => {
      const ring = v.section('ring', { of: c, layout: 'rings' });
      ring.track('sites').layer('sites', { type: 'ranges', data: fs });
      ring.track('far').layer('far', { type: 'ranges', data: ff });
    });
    return b.build();
  }
  const placed = async (dense: number, config = {}) => {
    const s = await settle(denseSpec(dense), [], config);
    const p = ringOf(await s.frame({ size: Size })).placement as RingPlacement;
    return Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]).sort();
  };

  it('more than six labels within 12 degrees are one grouped label with the count', async () => {
    const texts = await placed(13);
    expect(texts).toEqual(['13 items', 'a', 'b', 'c']);
  });

  it('six or fewer keep their own labels; a fan wider than the span splits into groups of its own', async () => {
    expect((await placed(6)).filter((t) => t.startsWith('Site')).length).toBe(6);
    // 4 bp apart is 0.54° a site: 13 sites span 6.5°; a 5° span groups the first 10 (over 6) and leaves the last 3 apart
    const texts = await placed(13, { rings: { cluster: { over: 6, degrees: 5 } } });
    expect(texts.filter((t) => t.endsWith('items')).length).toBe(1);
  });

  it('degrees 0 never groups', async () => {
    expect(
      (await placed(13, { rings: { cluster: { over: 6, degrees: 0 } } })).filter((t) => t.startsWith('Site')).length,
    ).toBe(13);
  });

  /** The label at the placed slot `k` of a set: a page point inside its box. */
  const labelPoint = (p: RingPlacement, set: RingPlacement['labels'], k: number) => ({
    x: set.right[k] ? set.x[k] + 3 : set.x[k] - 3,
    y: set.y[k] - 0.35 * p.labelSize,
  });
  const groupSlot = (p: RingPlacement) => {
    for (let k = 0; k < p.labels.m; k++) if (p.labels.of[p.labels.placed[k]] >= 0) return k;
    return -1;
  };

  it('both forms are placed with the layout: the grouped label and the fan of its members, nothing laid out on a hover', async () => {
    const s = await settle(denseSpec(13));
    const f = await s.frame({ size: Size });
    const p = ringOf(f).placement as RingPlacement;
    expect(p.fan.ng).toBe(1);
    expect(p.fan.groupSize[0]).toBe(13);
    const members = Array.from(p.fan.members.subarray(0, 13), (i) => p.fan.text[i]).sort();
    expect(members).toEqual(Array.from({ length: 13 }, (_, k) => `Site${k}`).sort());
    for (let q = 0; q < 13; q++) expect(p.fan.slot[p.fan.members[q]]).toBeGreaterThanOrEqual(0);
    const k = groupSlot(p);
    const hit = p.labels.hits[p.labels.placed[k]]!;
    expect(hit).toMatchObject({ target: 'group', count: 13 });
    expect(hit.members!.map((m) => m.element!.key).sort()).toEqual(
      Array.from({ length: 13 }, (_, q) => `s${q}`).sort(),
    );
    // a hover redraws the emphasis passes alone: the same layout, no rebuild
    const over = hitAtPage(f, labelPoint(p, p.labels, k).x, labelPoint(p, p.labels, k).y)!;
    expect(over.target).toBe('group');
    const hovered = s.renderer.frame({ size: Size, emphasis: { hover: { ...over } } });
    expect(hovered.rebuilt).toBe(false);
    expect(ringOf(hovered).placement).toBe(p);
  });

  it('a hover on a group of up to `fanMax` draws the fan (the grouped form at rest), and its labels are hit as their items', async () => {
    const s = await settle(denseSpec(10));
    const f = await s.frame({ size: Size });
    const p = ringOf(f).placement as RingPlacement;
    const rest = recordDrawLists(f.lists);
    expect(rest).toContain('10 items');
    expect(rest).not.toContain('Site3');
    const k = groupSlot(p);
    const at = labelPoint(p, p.labels, k);
    const over = hitAtPage(f, at.x, at.y)!;
    const hovered = s.renderer.frame({ size: Size, emphasis: { hover: { ...over } } });
    const lists = recordDrawLists(hovered.lists);
    for (let q = 0; q < 10; q++) expect(lists).toContain(`Site${q}`);
    expect(openGroupOf(p, over)).toBe(0);
    // a fan label is hit as its item, while the group is open
    const i = p.fan.members[3];
    const slot = p.fan.slot[i];
    const pt = labelPoint(p, p.fan, slot);
    // the hit reads the hover it is given, as the draw does: no draw between (`f` was never drawn with the hover)
    const member = hitView(f.layout, pt.x, pt.y, false, { ...over })!;
    expect(member.element?.key).toBe(p.fan.hits[i]!.element!.key);
    expect(hitView(f.layout, pt.x, pt.y)?.element?.key).not.toBe(p.fan.hits[i]!.element!.key);
    // the member's own hover keeps its group open
    expect(openGroupOf(p, { ...member })).toBe(0);
    // a hover on that item keeps the fan open; off the group the grouped form is drawn again
    const again = s.renderer.frame({ size: Size, emphasis: { hover: { ...member } } });
    expect(recordDrawLists(again.lists)).toContain('Site9');
    const off = s.renderer.frame({ size: Size, emphasis: {} });
    expect(recordDrawLists(off.lists)).not.toContain('Site3');
    expect(openGroupOf(p, undefined)).toBe(-1);
    // closed, the fan's place is the grouped labels' again
    expect(hitAtPage(off, pt.x, pt.y)?.target).not.toBe('group');
  });

  it('an export draws the grouped form alone', async () => {
    const out = await renderSvg(denseSpec(13), { size: Size });
    const texts = [...out.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    expect(texts).toContain('13 items');
    expect(texts.some((t) => t.startsWith('Site'))).toBe(false);
  });

  it('a fan across the origin (355° to 5°) is one group, its label at the origin', async () => {
    const n = 2686;
    const sites = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
    // 13 sites 3 bp apart from 2650 (355°) through the origin to 33 bp (4.4°)
    for (let k = 0; k < 13; k++) {
      const st = (2650 + 6 * k) % n;
      sites.start.push(st);
      sites.end.push(st + 3 > n ? n : st + 3);
      sites.id.push(`w${k}`);
      sites.label.push(`Wrap${k}`);
    }
    const b = createSeqSpec({ title: 'seam ring' });
    const c = b.data.sequence('plasmid', { alphabet: 'dna', letters: 'A'.repeat(n), cyclic: true });
    const fs = b.data.ranges('sites', { target: c, rows: sites });
    b.views.view('v', (v) => {
      v.section('ring', { of: c, layout: 'rings' }).track('sites').layer('sites', { type: 'ranges', data: fs });
    });
    const s = await settle(b.build());
    const p = ringOf(await s.frame({ size: Size })).placement as RingPlacement;
    expect(p.fan.ng).toBe(1);
    expect(p.fan.groupSize[0]).toBe(13);
    expect(p.labels.m).toBe(1);
    const a = p.fan.groupAngle[0];
    // the middle of 355° and about 4°: within a few degrees of the top, not the bottom (180°)
    const off = Math.min(a, 2 * Math.PI - a);
    expect(off).toBeLessThan((6 * Math.PI) / 180);
  });
});

describe('a ring over tiles', () => {
  function bigRing() {
    const b = createSeqSpec({ title: 'big ring' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const r = b.data.ranges('items', { target: chr, source: src });
    const v = b.data.values('signal', { target: chr, source: src });
    b.views.view('v', (view) => {
      const s = view.section('ring', { of: chr, layout: 'rings', range: [0, 2500], labels: 'none' });
      s.track('signal').layer('signal', { type: 'values', data: v });
      s.track('items').layer('items', { type: 'ranges', data: r });
    });
    return b.build();
  }

  it('a column over a tile edge reduces both tiles’ data; an item over the edge is one item', async () => {
    const s = await settle(bigRing(), [BigRuntime]);
    const f = await s.frame({ size: Size, settled: true });
    const p = ringOf(f).placement as RingPlacement;
    const lay = layerLay<ValuesLay>(f, 0);
    const a = lay.axis;
    // a column a few positions: reduced per angular column
    expect(lay.n).toBeGreaterThan(100);
    expect(lay.n).toBeLessThan(2500);
    let edges = 0;
    for (let j = 0; j < lay.n; j++) {
      // every column read (no tile missing), its high the maximum of the positions it holds
      expect(Number.isNaN(lay.v[j])).toBe(false);
      const u0 = a.unitAt(lay.x0[j] + 1e-6);
      const u1 = a.unitAt(lay.x1[j] - 1e-6);
      if (!(Math.floor(u0) < 1024 && Math.ceil(u1) > 1024) && !(Math.floor(u0) < 2048 && Math.ceil(u1) > 2048))
        continue;
      edges++;
      let inner = -Infinity;
      let outer = -Infinity;
      for (let q = Math.floor(u0); q <= Math.floor(u1); q++) outer = Math.max(outer, bigValue(q));
      for (let q = Math.ceil(u0); q < Math.floor(u1); q++) inner = Math.max(inner, bigValue(q));
      expect(lay.v[j]).toBeGreaterThanOrEqual(inner - 1e-6);
      expect(lay.v[j]).toBeLessThanOrEqual(outer + 1e-6);
    }
    expect(edges).toBeGreaterThanOrEqual(2);
    // `r1000` runs 1000 to 1070 over the edge at 1024: the same item either side of it
    const items = layerLay<RangesLay>(f, 1);
    // the first lane's middle (a tiled source's lanes centred by those in view)
    const band = items.laneY + items.laneHeight / 2;
    const keyAt = (u: number) => {
      p.warp.map(p.frameX(u + 0.5), band);
      return hitAtPage(f, p.warp.px, p.warp.py)?.element?.key;
    };
    expect(String(keyAt(1010))).toBe('r1000');
    expect(String(keyAt(1040))).toBe('r1000');
  });
});

describe('a linear view of a cyclic system across its origin', () => {
  it('draws the window past the total as its two runs, the units on over the origin; a follower the same', async () => {
    const b = createSeqSpec({ title: 'linear' });
    const c = b.data.sequence('construct', { alphabet: 'dna', letters: ringLetters(), cyclic: true });
    b.views.view('v', (v) => {
      const lin = v.section('lin', { of: c });
      lin.track('letters').layer('letters', { type: 'letters' });
      v.section('ring', { of: c, layout: 'rings', follow: lin, labels: 'none' })
        .track('r')
        .layer('r', { type: 'ruler' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: Size, viewports: { '/views/v/lin': { start: 590, end: 630 } } });
    const lin = f.layout.sections[0];
    const a = lin.axis!;
    expect([a.u0, a.u1]).toEqual([590, 630]);
    expect(a.visible([])).toEqual([590, 600, 0, 30]);
    // a letter a position, 40 of them, the unit past the total the origin's
    const y = lin.tracks[0].box.y + 2;
    expect(hitAtPage(f, a.xOfUnit(595.5), y)).toMatchObject({ position: 595 });
    expect(hitAtPage(f, a.xOfUnit(600.5), y)).toMatchObject({ position: 0 });
    expect(hitAtPage(f, a.xOfUnit(629.5), y)).toMatchObject({ position: 29 });
    // the ring's window is the linear view's range, over the origin
    const p = ringOf(f).placement as RingPlacement;
    expect([p.windowStart, p.windowEnd]).toEqual([590, 630]);
  });
});

describe('a ring with no strips', () => {
  const emptyRing = (range?: [number, number]) => {
    const b = createSeqSpec({ title: 'empty ring' });
    const src = b.data.source('never', { format: 'synthetic', uri: 'synthetic:never' });
    const c = range
      ? b.data.sequence('construct', { alphabet: 'dna', letters: ringLetters(), cyclic: true })
      : b.data.sequence('construct', { alphabet: 'dna', source: src, cyclic: true });
    b.views.view('v', (v) => {
      const ring = v.section('ring', { of: c, layout: 'rings', ...(range ? { range } : {}) });
      ring.track('ruler').layer('ruler', { type: 'ruler' });
      v.section('detail', { of: c, follow: ring }).track('letters').layer('letters', { type: 'letters' });
    });
    return b.build();
  };

  it('a pending ring draws its frame without its rings (no strip read past the last)', async () => {
    const clock = new FakeClock();
    const never: FormatRuntime = { name: 'never', formats: ['synthetic'], open: () => new Promise(() => {}) };
    const out = renderSvg(emptyRing(), { size: Size, runtimes: [never], clock, timeoutMs: 1000, window: true });
    await clock.advance(1000);
    const r = await out;
    expect(r.settled).toBe(false);
    expect(r.issues.map((i) => i.code)).toContain('seq.render.axis.pending');
  });

  it('an empty range draws its frame, and hits nothing', async () => {
    const s = await settle(emptyRing([400, 400]));
    const f = await s.frame({ size: Size, windows: true });
    const p = ringOf(f).placement as RingPlacement;
    expect(p.strips.n).toBe(0);
    // its frame drawn: the title in the hole, the lists written as SVG
    expect(f.lists.passes.labels.texts.n).toBeGreaterThan(0);
    expect(drawListsToSvg(f.lists)).toContain('<svg');
    expect(hitView(f.layout, p.cx, p.cy - p.r + 2)).toBeUndefined();
  });
});

/** The hit at a page point of a frame. */
function hitAtPage(f: Frame, x: number, y: number): Hit | undefined {
  return hitView(f.layout, x, y);
}
