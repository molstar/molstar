/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track's well (_spec/README.md "Render", Tracks layout and Shapes): its padding inside it, vertical only (the data's
// band inset by it: what layers lay out, hit and emphasize in), in rows, text and rings (radially), a grid's plot; the
// well a squircle by the squircle rule over its height (square below the minimum, or with `trackShape: square`), its
// data clipped to it with the same path in the lists and the SVG; and every emphasis a squircle: round where it really
// ends, square where cut or below the minimum, an outline concentric.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import {
  DrawLists,
  drawListsToSvg,
  type Emphasis,
  type Frame,
  hitView,
  type PassName,
  recordDrawLists,
  renderSvg,
} from '../index.js';
import { P53, p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';
import { gridSpec } from './pairs-fixtures.js';
import { blocksSpec } from './rings-fixtures.js';

const Size = { width: 800, height: 400 };
const Protein = '/data/p53';

type Padding = number | { top?: number; bottom?: number };

/** p53's letters and domains as tracks, each track's padding (and height) as given, on one section or as rings. */
function padded(letters?: Padding, domains?: Padding, options: { layout?: 'rings' | 'text'; height?: number } = {}) {
  const b = createSeqSpec({ title: 'wells' });
  const p53 = b.data.sequence('p53', {
    alphabet: 'protein',
    letters: P53,
    ...(options.layout === 'rings' ? { cyclic: true } : {}),
  });
  const d = b.data.ranges('domains', {
    target: p53,
    rows: { start: [10, 100], end: [60, 290], label: ['TAD', 'DBD'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53, ...(options.layout ? { layout: options.layout } : {}) });
    s.track('letters', {
      label: 'p53',
      ...(letters !== undefined ? { padding: letters } : {}),
      ...(options.height !== undefined ? { height: options.height } : {}),
    }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains', ...(domains !== undefined ? { padding: domains } : {}) }).layer('ranges', {
      type: 'ranges',
      data: d,
    });
  });
  return b.build();
}

async function frameOf(spec: Node, emphasis?: Emphasis, config = {}, viewport = { start: 0, end: 30 }): Promise<Frame> {
  const s = await settle(spec, [], config);
  const input = { size: Size, viewport };
  await s.frame(input);
  return s.renderer.frame({ ...input, ...(emphasis ? { emphasis } : {}) });
}

/** The spec with a track's params added (found by its id). */
function withParams(spec: Node, track: string, params: Record<string, unknown>): Node {
  const walk = (n: Node): Node =>
    n.kind === 'track' && n.id === track
      ? { ...n, params: { ...n.params, ...params } as Node['params'] }
      : n.children
        ? { ...n, children: n.children.map(walk) }
        : n;
  return walk(spec);
}

const lines = (f: { lists: DrawLists }, pass: PassName) => recordDrawLists(f.lists, [pass]).split('\n').slice(2);

describe('a track’s padding', () => {
  it('moves the data’s band inside the well and not the well’s top; the next track below the well', async () => {
    const plain = await frameOf(padded());
    const pad = await frameOf(padded(8));
    const [a0, b0] = plain.layout.sections[0].tracks;
    const [a, b] = pad.layout.sections[0].tracks;
    // the well starts where the track did; the data 8 px down, as tall as before
    expect(a.well.y).toBe(a0.box.y);
    expect(a.box.y).toBe(a0.box.y + 8);
    expect(a.box.height).toBe(a0.box.height);
    expect(a.well.height).toBe(a0.box.height + 16);
    // the domains below the grown well
    expect(b.box.y).toBe(b0.box.y + 16);
    // the band the layers lay out in: the letters' cells on the data band
    const cell = lines(pad, 'data').find((l) => l.startsWith('rect ') && !l.startsWith('rect 0 '))!;
    expect(Number(cell.split(' ')[2])).toBe(a.box.y);
  });

  it('a top-only padding grows the well above the data alone', async () => {
    const f = await frameOf(padded({ top: 12 }));
    const t = f.layout.sections[0].tracks[0];
    expect(t.box.y - t.well.y).toBe(12);
    expect(t.well.y + t.well.height).toBe(t.box.y + t.box.height);
    const g = await frameOf(padded({ bottom: 6 }));
    const u = g.layout.sections[0].tracks[0];
    expect(u.box.y).toBe(u.well.y);
    expect(u.well.height - u.box.height).toBe(6);
  });

  it('a hit at the padding is no item; the data band hits as before', async () => {
    const f = await frameOf(padded(8));
    const t = f.layout.sections[0].tracks[0];
    const x = t.box.x + 30;
    expect(hitView(f.layout, x, t.well.y + 3)).toBeUndefined();
    expect(hitView(f.layout, x, t.box.y + t.box.height + 3)).toBeUndefined();
    expect(hitView(f.layout, x, t.box.y + 3)?.position).toBeDefined();
  });

  it('the render config’s `trackPadding` is the default; a track’s own wins', async () => {
    const f = await frameOf(padded(undefined, 2), undefined, { trackPadding: { top: 4, bottom: 4 } });
    const [a, b] = f.layout.sections[0].tracks;
    expect(a.box.y - a.well.y).toBe(4);
    expect(b.box.y - b.well.y).toBe(2);
  });

  it('a ring’s padding is radial: its ground the well, its band inset by it', async () => {
    const plain = await frameOf(padded(undefined, undefined, { layout: 'rings' }));
    const f = await frameOf(padded(6, undefined, { layout: 'rings' }));
    const sp = f.layout.sections[0].placement.strips;
    const s0 = plain.layout.sections[0].placement.strips;
    expect(sp.padT[0]).toBeGreaterThan(0);
    expect(sp.padT[0]).toBe(sp.padB[0]);
    // the outer ring's band moves in by its top padding, the next ring by both sides (scaled together)
    expect(sp.bandY[0]).toBeCloseTo(sp.padT[0], 6);
    expect(sp.bandY[1]).toBeGreaterThan(s0.bandY[1]);
    // its ground is the well: the band grown by its padding, a sector (a path)
    expect(sp.bandY[0] - sp.padT[0]).toBe(0);
  });

  it('text bands honour it: the band inset in its line', async () => {
    const plain = await frameOf(padded(undefined, undefined, { layout: 'text' }));
    const f = await frameOf(padded({ top: 5 }, undefined, { layout: 'text' }));
    const a = plain.layout.sections[0].placement.strips;
    const b = f.layout.sections[0].placement.strips;
    expect(b.bandY[0]).toBe(a.bandY[0] + 5);
    expect(b.padT[0]).toBe(5);
    expect(b.bandH[0]).toBe(a.bandH[0]);
  });

  it('a grid’s plot is inset by it, its well around it', async () => {
    const f = await frameOf(withParams(gridSpec({ n: 40 }), 'plot', { padding: { top: 10 } }));
    const t = f.layout.sections[1].tracks[0];
    expect(t.box.y - t.well.y).toBe(10);
    expect(f.layout.sections[1].placement.strips.padT[0]).toBe(10);
  });

  it('is validated: non-negative numbers, one or both sides, else an issue', () => {
    expect(checkSeqSpec(padded(4, { top: 2, bottom: 3 })).filter((i) => i.severity === 'error')).toEqual([]);
    expect(checkSeqSpec(padded({ top: 2 })).filter((i) => i.severity === 'error')).toEqual([]);
    for (const bad of [-1, { top: -2 }, { left: 4 }, 'wide'] as unknown as Padding[]) {
      const issues = checkSeqSpec(padded(bad)).filter((i) => i.severity === 'error');
      expect(issues.length, JSON.stringify(bad)).toBeGreaterThan(0);
      expect(
        issues.some((i) => i.address?.includes('letters')),
        JSON.stringify(bad),
      ).toBe(true);
    }
  });
});

describe('the well a squircle, its data clipped to it', () => {
  it('the well’s radius by its height: a third, capped at 8, square below 6 px (the lists and the SVG)', async () => {
    for (const [height, radius] of [
      [4, undefined],
      [12, '4'],
      [20, '6.67'],
      [36, '8'],
    ] as const) {
      const r = await renderSvg(padded(undefined, undefined, { height }), {
        size: Size,
        viewport: { start: 0, end: 30 },
      });
      const well = recordDrawLists(r.lists, ['background']).split('\n')[3];
      expect(well.startsWith('rect 128 8 '), well).toBe(true);
      if (radius) expect(well.endsWith(`#e9e6e0ff squircle ${radius}`), well).toBe(true);
      else expect(well.endsWith('#e9e6e0ff'), well).toBe(true);
      // the SVG: a path of four corners (a cubic each), or a rect
      const own = new DrawLists();
      const w = own.writer('background');
      const p = well.split(' ').map(Number);
      w.rect(p[1], p[2], p[3], p[4], 0xe9e6e0ff, 0, 0, radius ? Number(radius) : 0, 3);
      const svg = drawListsToSvg(own);
      if (radius) expect((svg.match(/<path d="[^"]*"/)?.[0].match(/C/g) ?? []).length).toBe(4);
      else expect(svg).toContain(`<rect x="128" y="8" width="${p[3]}" height="${height}"`);
    }
  });

  it('a letters track’s data clipped to its well’s squircle, one clip a track (the lists and the SVG)', async () => {
    const r = await renderSvg(p53TracksSpec(), { size: Size, viewport: { start: 0, end: 30 } });
    const clips = lines(r, 'data').filter((l) => l.startsWith('clip '));
    expect(clips).toEqual([
      'clip 128 8 664 20 squircle 6.67',
      'clip 128 32 664 36 squircle 8',
      'clip 128 72 664 32 squircle 8',
    ]);
    // the SVG's clip paths: the same curve, a cubic a corner
    const paths = [...r.svg.matchAll(/<clipPath id="[^"]*"><path d="([^"]*)"\/><\/clipPath>/g)].map((m) => m[1]);
    expect(paths.length).toBe(3);
    for (const d of paths) expect(d.match(/C/g)?.length).toBe(4);
    expect(paths[0].startsWith('M 134.67 8 L ')).toBe(true);
  });

  it('a grid track’s plot clipped to its well', async () => {
    const r = await renderSvg(gridSpec({ n: 40 }), { size: Size });
    const clips = lines(r, 'data').filter((l) => l.startsWith('clip '));
    // the letters above, then the plot: one clip each, its squircle by its height
    expect(clips).toEqual(['clip 128 8 664 20 squircle 6.67', 'clip 128 36 664 327.6 squircle 8']);
  });

  it('`trackShape: square` draws the well and its clip as boxes, as before', async () => {
    const r = await renderSvg(p53TracksSpec(), {
      size: Size,
      viewport: { start: 0, end: 30 },
      config: { trackShape: 'square' },
    });
    expect(lines(r, 'background').slice(1, 4)).toEqual([
      'rect 128 8 664 20 #e9e6e0ff',
      'rect 128 32 664 36 #e9e6e0ff',
      'rect 128 72 664 32 #e9e6e0ff',
    ]);
    expect(lines(r, 'data').filter((l) => l.startsWith('clip '))).toEqual([
      'clip 128 8 664 20',
      'clip 128 32 664 36',
      'clip 128 72 664 32',
    ]);
    expect(r.svg).not.toMatch(/<clipPath id="[^"]*"><path/);
  });

  it('a block inside the row: its clip square where it stops short of the well’s ends', async () => {
    const f = await frameOf(blocksSpec(), undefined, {}, { start: 0, end: 400 });
    const clips = lines(f, 'data').filter((l) => l.startsWith('clip '));
    expect(clips.some((l) => / squircle [\d.]+ left$/.test(l) || / squircle [\d.]+$/.test(l))).toBe(true);
    // a middle block reaches neither end: a box
    expect(clips.some((l) => !l.includes('squircle'))).toBe(true);
  });
});

describe('every emphasis a squircle', () => {
  const Letters = '/views/main/protein/letters/letters';

  it('the projected hover: round where its range really ends, square where the view cuts it', async () => {
    const f = await frameOf(p53TracksSpec(), {
      projected: { hover: rangesLocation(Protein, [4, 10]), hoverSection: '/views/other' },
    });
    const hover = lines(f, 'hover').filter((l) => l.startsWith('rect '));
    expect(hover.length).toBeGreaterThan(0);
    // over the column (the three tracks and their gaps): 96 px tall, so the radius capped (8) plus the 1.75 out
    expect(hover.every((l) => / squircle 9\.75$/.test(l))).toBe(true);
    const cut = await frameOf(
      p53TracksSpec(),
      {
        projected: { hover: rangesLocation(Protein, [0, 10]), hoverSection: '/views/other' },
      },
      {},
      { start: 2, end: 32 },
    );
    expect(
      lines(cut, 'hover')
        .filter((l) => l.startsWith('rect '))
        .every((l) => / squircle 9\.75 right$/.test(l)),
    ).toBe(true);
  });

  it('the highlights in every style: a tint, a band (its bar below the minimum: square), an outline', async () => {
    const on = (range: [number, number], style?: 'tint' | 'band' | 'outline') => ({
      location: rangesLocation(Protein, range),
      color: 0xff0000ff,
      ...(style ? { style } : {}),
    });
    const f = await frameOf(p53TracksSpec(), {
      highlights: [on([2, 5]), on([8, 11], 'band'), on([14, 17], 'outline')],
    });
    const hl = lines(f, 'highlight').filter((l) => l.startsWith('rect '));
    expect(hl[0]).toMatch(/ #ff00004d squircle 8$/);
    expect(hl[1]).toMatch(/ #ff0000ff$/);
    expect(hl[2]).toMatch(/ none stroke #ff0000ff 1 squircle 7\.5$/);
    // overlapping channels: each stripe a squircle by its own height
    const two = await frameOf(p53TracksSpec(), { highlights: [on([2, 8]), { ...on([4, 10]), color: 0x00ff00ff }] });
    const stripes = lines(two, 'highlight').filter((l) => l.startsWith('rect ') && l.split(' ')[4] === '48');

    expect(stripes.length).toBe(2);
    expect(stripes.every((l) => / squircle 8 (left|right)$/.test(l))).toBe(true);
  });

  it('a letter cell’s hover: concentric with the cell’s squircle; square on a side the view cuts and below 6 px', async () => {
    const hover = (position: number) => ({ layer: Letters, section: 'protein', system: Protein, position });
    const f = await frameOf(p53TracksSpec(), { hover: hover(4) });
    expect(lines(f, 'hover').filter((l) => l.startsWith('rect '))).toEqual([
      'rect 215.25 6.25 25.5 23.5 none stroke #2d271eff 1.5 squircle 8.42',
    ]);
    // the cell cut by the view's left edge: square there (its visible 11 px: the radius half of it)
    const cut = await frameOf(p53TracksSpec(), { hover: hover(2) }, {}, { start: 2.5, end: 32.5 });
    expect(lines(cut, 'hover').filter((l) => l.startsWith('rect '))[0]).toMatch(/ squircle 7\.25 right$/);
    // a cell narrower than the minimum (a far zoom): square
    const far = await frameOf(p53TracksSpec(), { hover: hover(40) }, {}, { start: 0, end: 393 });
    expect(lines(far, 'hover').filter((l) => l.startsWith('rect '))[0]).not.toContain('squircle');
  });

  it('a values cell’s hover and a logo cell’s: squircles too', async () => {
    const f = await frameOf(p53TracksSpec(), {
      hover: { layer: '/views/main/protein/conservation/values', section: 'protein', system: Protein, position: 4 },
    });
    expect(lines(f, 'hover').filter((l) => l.startsWith('rect '))[0]).toMatch(/ squircle 9\.75$/);
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: { label: ['a', 'b'], letters: ['ACDEFG', 'ACDEFA'] },
    });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) =>
        s.track('logo', { label: 'logo' }, (t) => t.layer('l', { type: 'logo', data: aln } as never)),
      ),
    );
    const g = await frameOf(
      b.build(),
      { hover: { layer: '/views/v/s/logo/l', section: 's', system: '/data/aln', position: 2 } },
      {},
      {
        start: 0,
        end: 6,
      },
    );
    expect(lines(g, 'hover').filter((l) => l.startsWith('rect '))[0]).toMatch(/ squircle 9\.75$/);
  });

  it('a chip’s hover: concentric with its pill', async () => {
    const f = await frameOf(
      blocksSpec(),
      {
        hover: { layer: '/views/v/complex/letters/letters', section: 'complex', system: '/data/atp', position: 0 },
      },
      {},
      { start: 0, end: 400 },
    );
    const h = lines(f, 'hover').filter((l) => l.startsWith('rect '));
    expect(h.length).toBe(1);
    // the pill is a pixel inside the 20 px cell: its radius 9, the outline 1 + 1.75 further out
    expect(h[0]).toMatch(/ squircle 11\.75$/);
  });
});
