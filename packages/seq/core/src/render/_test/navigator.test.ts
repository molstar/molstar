/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigators (_spec/README.md "Render", Navigators): the picture drawn against the data (a strip of the tracks, an
// ideogram's every band with its mark and hit, the rows summarized), each window the image of its section's viewport
// (and back within a device pixel), the floor, the parts, the reachable span, the wrap over a cyclic origin.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, miniatureNavigator, overviewNavigator, stripNavigator, theme } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { type BandsLay, NavWindowFlags, navigatorPart, renderConfig, RowsLay, type ValuesLay } from '../index.js';
import { navigatorBox } from '../navigator.js';
import { P53, p53GeneSpec } from './fixtures.js';
import { Bands, ideogramSpec, settleNavigator } from './navigator-fixtures.js';
import { rowsSpec, Synthetic } from './rows-fixtures.js';
import { readRects } from './small-helpers.js';

const Size = { width: 800, height: 200 };
const NavSize = { width: 800, height: 48 };

/** p53 as tracks with a navigator made by `form` (a cyclic system, a reversed section, as asked). */
function p53Strip(
  form: (of: string) => Record<string, unknown>,
  options: { cyclic?: boolean; reversed?: boolean } = {},
): Node {
  const b = createSeqSpec({ title: 'p53 strip' });
  const p53 = b.data.sequence('p53', {
    alphabet: 'protein',
    letters: P53,
    ...(options.cyclic ? { cyclic: true } : {}),
  });
  const d = b.data.ranges('domains', { target: p53, rows: { start: [0, 100, 300], end: [60, 290, 393] } });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p53, ...(options.reversed ? { direction: 'reversed' as const } : {}) });
    s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('domains').layer('ranges', { type: 'ranges', data: d });
    v.navigator('nav', form(s.address) as never);
  });
  return b.build();
}

describe('the navigator node', () => {
  it('resolves its bound sections, its system and its paths, its form’s defaults', () => {
    const model = resolveSpec(p53Strip((of) => stripNavigator(of, { content: 'tracks' })));
    const [view] = resolveViews(model);
    const n = view.navigators[0];
    expect(n.address).toBe('/views/v/navigators/nav');
    expect(n.of).toEqual(['/views/v/protein']);
    expect(n.missing).toEqual([]);
    expect(n.system).toBe('/data/p53');
    expect(n.paths['/views/v/protein'].edges).toEqual([]);
    expect([n.axes, n.content, n.place, n.hidden]).toEqual(['positions', 'tracks', 'top', false]);
    expect(navigatorBox(n, 800, 300)).toEqual({ width: 800, height: 48 });
  });

  it('a section it names that the view lacks is said, and no window is drawn for it', async () => {
    const spec = p53Strip((of) => stripNavigator([of, '/views/v/gone']));
    const s = await settleNavigator(spec);
    expect(s.view.navigators[0].missing).toEqual(['/views/v/gone']);
    const { nav } = await s.frames({ size: Size }, { size: NavSize });
    expect(nav.issues.map((i) => i.code)).toContain('seq.render.navigator.of');
    expect(nav.windows.length).toBe(1);
  });

  it('the forms’ sizes', () => {
    const model = resolveSpec(rowsSpec(100, false));
    const [view] = resolveViews(model);
    const base = { ...view.navigators[0], of: [], missing: [], paths: {}, hidden: false } as never;
    expect(
      navigatorBox({ ...(base as object), axes: 'rows', content: 'rows', place: 'right' } as never, 800, 300),
    ).toEqual({ width: 96, height: 300 });
    expect(
      navigatorBox({ ...(base as object), axes: 'both', content: 'tracks', place: 'top' } as never, 800, 300),
    ).toEqual({ width: 160, height: 96 });
    expect(
      navigatorBox({ ...(base as object), axes: 'positions', content: 'bands', place: 'top' } as never, 800, 300),
    ).toEqual({ width: 800, height: 24 });
  });
});

describe('a strip', () => {
  it('draws the section’s tracks small, but its ruler, and its window is the viewport’s image', async () => {
    const s = await settleNavigator(p53Strip((of) => stripNavigator(of, { content: 'tracks' })));
    const viewport = { start: 100, end: 200 };
    const { view, nav } = await s.frames({ size: Size, viewport }, { size: NavSize, moves: { '/views/v/protein': 1 } });
    const pic = s.nav.pictureSection!;
    expect(pic.tracks.map((t) => t.track.address.split('/').pop())).toEqual(['letters', 'domains']);
    // each track a quarter of its height, at least 2 px
    expect(pic.tracks[1].box.height).toBe(36 / 4);
    const w = nav.windows[0];
    expect(w.shown).toBe(true);
    const axis = pic.axis!;
    expect(w.x0).toBeCloseTo(axis.xOfUnit(100), 9);
    expect(w.x1).toBeCloseTo(axis.xOfUnit(200), 9);
    // and back: the window's ends carried to the section's units within a device pixel
    expect(Math.abs(s.nav.sectionUnitAt('/views/v/protein', w.x0)! - 100)).toBeLessThan(1 / axis.scale);
    expect(Math.abs(s.nav.sectionUnitAt('/views/v/protein', w.x1)! - 200)).toBeLessThan(1 / axis.scale);
    expect(view.layout.sections[0].axis!.u0).toBe(100);
    // the window's outline in the window pass, grips each side
    const rects = readRects(nav.lists, 'window');
    expect(rects.some((r) => Math.abs(r.x - w.dx0) < 1e-3 && Math.abs(r.w - (w.dx1 - w.dx0)) < 1e-3)).toBe(true);
  });

  it('every domain in the picture has its mark (drawing against the data at the fit)', async () => {
    const s = await settleNavigator(p53Strip((of) => stripNavigator(of, { content: 'tracks' })));
    await s.frames({ size: Size, viewport: { start: 0, end: 50 } }, { size: NavSize });
    const pic = s.nav.pictureSection!;
    const domains = pic.tracks[1].layers[0].lay as { m: number; rows: { n: number } };
    expect(domains.rows.n).toBe(3);
    expect(domains.m).toBeGreaterThanOrEqual(3);
  });

  it('a floored window is drawn at least 16 px, centred on the truth, with a hairline at its true extent', async () => {
    const s = await settleNavigator(
      p53Strip((of) => stripNavigator(of)),
      0,
      [],
      {},
    );
    const { nav } = await s.frames(
      { size: Size, viewport: { start: 196, end: 197 } },
      { size: { width: 800, height: 20 }, moves: { '/views/v/protein': 1 } },
    );
    const w = nav.windows[0];
    expect(w.dx1 - w.dx0).toBeCloseTo(16, 9);
    expect((w.dx0 + w.dx1) / 2).toBeCloseTo((w.x0 + w.x1) / 2, 9);
    // the true extent covered by the drawn one
    expect(w.dx0).toBeLessThanOrEqual(w.x0);
    expect(w.dx1).toBeGreaterThanOrEqual(w.x1);
  });

  it('idle is a quiet outline: no wash, no grips', async () => {
    const s = await settleNavigator(p53Strip((of) => stripNavigator(of)));
    const busy = await s.frames(
      { size: Size, viewport: { start: 0, end: 100 } },
      { size: NavSize, moves: { '/views/v/protein': 1 } },
    );
    const n = readRects(busy.nav.lists, 'window').length;
    const idle = s.nav.frame({ size: NavSize, view: busy.view.layout, moves: { '/views/v/protein': 0 } });
    expect(idle.rebuilt).toBe(false);
    expect(readRects(idle.lists, 'window').length).toBe(1);
    expect(n).toBeGreaterThan(1);
  });
});

describe('the window’s parts', () => {
  const c = renderConfig().navigator.window;
  it('ends as tabs outside, at most a quarter inside; the body; outside; past the span nothing', () => {
    expect(navigatorPart(97, 100, 200, 0, 800, c)).toBe('start');
    expect(navigatorPart(103, 100, 200, 0, 800, c)).toBe('start');
    expect(navigatorPart(150, 100, 200, 0, 800, c)).toBe('body');
    expect(navigatorPart(203, 100, 200, 0, 800, c)).toBe('end');
    expect(navigatorPart(400, 100, 200, 0, 800, c)).toBe('outside');
    expect(navigatorPart(900, 100, 200, 0, 800, c)).toBeUndefined();
    // a 16 px window: its middle still pans
    expect(navigatorPart(108, 100, 116, 0, 800, c)).toBe('body');
    // a grip hits over half the hit width outside its end (the ring's window alike), by touch half the touch hit
    expect(navigatorPart(94, 100, 200, 0, 800, c)).toBe('start');
    expect(navigatorPart(93.5, 100, 200, 0, 800, c)).toBe('outside');
    expect(navigatorPart(88.5, 100, 200, 0, 800, c, true)).toBe('start');
  });

  it('a hit on the navigator names its section, part and position', async () => {
    const s = await settleNavigator(p53Strip((of) => stripNavigator(of)));
    const { nav } = await s.frames({ size: Size, viewport: { start: 100, end: 200 } }, { size: NavSize });
    const w = nav.windows[0];
    const mid = (w.dx0 + w.dx1) / 2;
    expect(s.nav.hit(mid, 10)).toMatchObject({ section: '/views/v/protein', part: 'body', axis: 'positions' });
    expect(s.nav.hit(w.dx0 - 2, 10)?.part).toBe('start');
    const out = s.nav.hit(700, 10)!;
    expect(out.part).toBe('outside');
    expect(out.position).toBe(s.nav.pictureSection!.axis!.positionAt(700));
  });
});

describe('an ideogram', () => {
  it('every band of the chromosome has its mark and its hit at the fit; the centromere two triangles', async () => {
    const s = await settleNavigator(ideogramSpec());
    const { nav } = await s.frames({ size: Size }, { size: { width: 800, height: 24 } });
    const pic = s.nav.pictureSection!;
    const lay = pic.tracks[0].layers[0].lay as BandsLay;
    expect(lay.rows.n).toBe(Bands.start.length);
    for (let i = 0; i < Bands.start.length; i++) {
      const x = pic.axis!.x(Bands.start[i]) + 1;
      const h = s.nav.hit(x, 12)!;
      expect(h.band).toBe(Bands.id[i]);
    }
    const acen = [...lay.tip.subarray(0, lay.rows.n)];
    expect(acen.filter((t) => t).length).toBe(2);
    // the p arm's points to its end (right), the q arm's to its start (left)
    expect(acen[3]).toBe(1);
    expect(acen[4]).toBe(2);
    expect(nav.issues.filter((i) => i.severity === 'error')).toEqual([]);
  });

  it('the region’s window within the chromosome, the rest of it washed as unreachable', async () => {
    const s = await settleNavigator(ideogramSpec());
    const { nav } = await s.frames(
      { size: Size },
      { size: { width: 800, height: 24 }, moves: { '/views/v/region': 1 } },
    );
    const axis = s.nav.pictureSection!.axis!;
    const w = nav.windows[0];
    expect(w.x0).toBeCloseTo(axis.x(300), 6);
    expect(w.x1).toBeCloseTo(axis.x(360), 6);
    expect(w.lo).toBeCloseTo(axis.x(300), 6);
    expect(w.hi).toBeCloseTo(axis.x(360), 6);
    // a drag maps back into the region's units
    expect(s.nav.sectionUnitAt('/views/v/region', axis.x(330))).toBeCloseTo(30, 6);
  });

  it('on a reversed section the ideogram runs reversed: its tips turn', async () => {
    const s = await settleNavigator(ideogramSpec({ reversed: true }));
    await s.frames({ size: Size }, { size: { width: 800, height: 24 } });
    const pic = s.nav.pictureSection!;
    expect(pic.axis!.dir).toBe(-1);
    const lay = pic.tracks[0].layers[0].lay as BandsLay;
    expect([lay.tip[3], lay.tip[4]]).toEqual([2, 1]);
    expect(pic.axis!.x(100)).toBeGreaterThan(pic.axis!.x(900));
  });
});

describe('a window over a cyclic origin', () => {
  const P = '/views/v/protein';
  it('two pieces, a grip at each true end; where a piece meets the axis’s end it is body; a drag maps across it', async () => {
    const s = await settleNavigator(p53Strip((of) => stripNavigator(of), { cyclic: true }));
    const { nav } = await s.frames(
      { size: Size, viewport: { start: 380, end: 410 } },
      { size: NavSize, moves: { [P]: NavWindowFlags.x } },
    );
    const w = nav.windows[0];
    const axis = s.nav.pictureSection!.axis!;
    expect(w.wrap).toBe(true);
    const right = axis.x0 + axis.width;
    expect(w.x0).toBeCloseTo(axis.xOfUnit(380), 6);
    expect(w.w1).toBeCloseTo(axis.xOfUnit(17), 6);
    // the parts: at the axis's two ends body, the grips at the window's true ends
    expect(s.nav.hit(right - 2, 10)?.part).toBe('body');
    expect(s.nav.hit(axis.x0 + 1, 10)?.part).toBe('body');
    expect(s.nav.hit(w.x0 - 1, 10)?.part).toBe('start');
    expect(s.nav.hit(w.w1 + 1, 10)?.part).toBe('end');
    // the grips drawn: one before the start, one after the true end, none at the axis's ends
    const c = renderConfig().navigator.window;
    const grips = readRects(nav.lists, 'window').filter((r) => Math.abs(r.w - c.handlePx) < 1e-6);
    expect(grips.map((r) => r.x).sort((a, b) => a - b)).toEqual(
      [w.w1, w.x0 - c.handlePx].map((x) => expect.closeTo(x, 3)),
    );
    // a drag across the origin: past the axis's end the section's units keep going
    expect(s.nav.sectionUnitAt(P, axis.xOfUnit(390))).toBeCloseTo(390, 0);
    expect(s.nav.sectionUnitAt(P, axis.xOfUnit(5))).toBeCloseTo(5, 0);
  });

  it('`direction: genomic` on a reversed cyclic section: a window not over the origin is one piece', async () => {
    const s = await settleNavigator(
      p53Strip((of) => stripNavigator(of), { cyclic: true, reversed: true }),
      0,
      [],
      {},
      () => ({ direction: 'genomic' }),
    );
    const { nav } = await s.frames(
      { size: Size, viewport: { start: 100, end: 150 } },
      { size: NavSize, moves: { [P]: NavWindowFlags.x } },
    );
    const w = nav.windows[0];
    const axis = s.nav.pictureSection!.axis!;
    expect(axis.dir).toBe(1);
    expect(w.wrap).toBe(false);
    expect((w.x1 - w.x0) / axis.scale).toBeCloseTo(50, 0);
  });
});

describe('the picture’s colors', () => {
  it('are the view’s: an open scale’s held domain, never one of its own over the whole axis', async () => {
    const b = createSeqSpec();
    const seq = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(100) });
    const values = b.data.values('a', { target: seq, values: Float32Array.from({ length: 100 }, (_, i) => i) });
    const sc = b.views.scale('cons', { label: 'Conservation' });
    b.views.view('v', (v) => {
      const sec = v.section('p', { of: seq });
      sec
        .track('a')
        .layer('a', { type: 'values', data: values, scale: sc, color: theme.field('value', { scale: sc }) });
      v.navigator('strip', stripNavigator(sec, { content: 'tracks' }));
    });
    const s = await settleNavigator(b.build(), 0, [], {}, (v) => ({ domains: v.renderer.layoutState }));
    const { view } = await s.frames(
      { size: Size, viewport: { start: 0, end: 20 }, settled: true },
      { size: NavSize, settled: true },
    );
    const viewLay = view.layout.sections[0].tracks[0].layers[0].lay as ValuesLay;
    const picLay = s.nav.pictureSection!.tracks[0].layers[0].lay as ValuesLay;
    expect(viewLay.domain).toEqual([0, 19]);
    expect(picLay.domain).toEqual(viewLay.domain);
  });
});

describe('a follower’s window', () => {
  it('a protein’s strip draws its gene follower’s window, carried along the codon path', async () => {
    const spec = p53GeneSpec();
    // a strip over the protein bound to both sections
    (spec.children![1].children![0].children as Node[]).push({
      kind: 'navigators',
      id: 'navigators',
      children: [
        { kind: 'navigator', id: 'nav', params: { of: ['/views/main/protein', '/views/main/gene'], content: 'axis' } },
      ],
    });
    const s = await settleNavigator(spec);
    const { nav } = await s.frames(
      { size: Size, viewport: { start: 50, end: 150 } },
      { size: NavSize, moves: { '/views/main/protein': 1, '/views/main/gene': 1 } },
    );
    const [p, g] = nav.windows;
    expect(p.shown && g.shown).toBe(true);
    // the gene follows the protein: its window is the protein's, a codon's edge apart at most
    const px = 800 / 393;
    expect(Math.abs(g.x0 - p.x0)).toBeLessThan(px + 1e-6);
    expect(Math.abs(g.x1 - p.x1)).toBeLessThan(px + 1e-6);
    // past the picture's end a drag keeps going as the path carries the edge: a protein position three nucleotides
    const axis = s.nav.pictureSection!.axis!;
    const end = axis.x0 + axis.width;
    const at = s.nav.sectionUnitAt('/views/main/gene', end - 1e-6)!;
    const past = s.nav.sectionUnitAt('/views/main/gene', end + 10 * axis.scale)!;
    expect(Math.abs(past - at)).toBeCloseTo(30, 0);
  });
});

describe('a vertical overview and a miniature', () => {
  it('the rows summarized down the column, the window the rows in view', async () => {
    const spec = rowsSpec(400, false);
    (spec.children![1].children![0].children as Node[]).push({
      kind: 'navigators',
      id: 'navigators',
      children: [{ kind: 'navigator', id: 'nav', params: overviewNavigator('/views/v/s') as never }],
    });
    const s = await settleNavigator(spec, 0, [Synthetic]);
    const { view, nav } = await s.frames(
      { size: { width: 800, height: 600 } },
      { size: { width: 96, height: 300 }, moves: { '/views/v/s': NavWindowFlags.y } },
    );
    const pic = s.nav.pictureSection!;
    const lay = pic.tracks[0].layers[0].lay as RowsLay;
    expect(lay.display.count).toBe(400);
    // every row in the band
    expect(lay.rowScale * 400).toBeCloseTo(lay.box.height, 6);
    const vr = view.layout.sections[0].lays[0][0].lay as RowsLay;
    const w = nav.windows[0];
    const shown = (vr.box.height - vr.pinnedHeight) / vr.rowScale;
    expect((w.y1 - w.y0) / lay.rowScale).toBeCloseTo(Math.min(400, shown), 6);
    expect(w.x0).toBe(pic.placement.box.x);
    expect(s.nav.rowAt(w.y0)).toBeCloseTo(vr.top, 6);
  });

  it('a lazy alignment’s rows are left out, said', async () => {
    const spec = rowsSpec(200_000, false);
    (spec.children![1].children![0].children as Node[]).push({
      kind: 'navigators',
      id: 'navigators',
      children: [{ kind: 'navigator', id: 'nav', params: overviewNavigator('/views/v/s') as never }],
    });
    const s = await settleNavigator(spec, 0, [Synthetic]);
    const { nav } = await s.frames({ size: { width: 800, height: 600 } }, { size: { width: 96, height: 300 } });
    expect(nav.issues.map((i) => i.code)).toContain('seq.render.navigator.rows');
    expect(s.nav.pictureSection!.tracks.length).toBe(0);
  });

  it('a miniature covers both axes', async () => {
    const spec = rowsSpec(400, false);
    (spec.children![1].children![0].children as Node[]).push({
      kind: 'navigators',
      id: 'navigators',
      children: [{ kind: 'navigator', id: 'nav', params: miniatureNavigator('/views/v/s') as never }],
    });
    const s = await settleNavigator(spec, 0, [Synthetic]);
    const { nav } = await s.frames(
      { size: { width: 800, height: 600 }, viewport: { start: 0, end: 100 } },
      { size: { width: 160, height: 96 }, moves: { '/views/v/s': 3 } },
    );
    const w = nav.windows[0];
    expect(w.axes).toBe(3);
    expect(w.x1 - w.x0).toBeCloseTo((160 * 100) / 300, 6);
    expect(w.y1).toBeLessThanOrEqual(96);
  });
});
