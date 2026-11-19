/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text layout: lines of whole groups, a number a line, group labels above, a band no layer used on a line
// collapsed there, the placement as the inverse of hit-testing over random widths and scales, drawing checked against
// the data across a line break, a reflow holding the anchor's line, an export drawing every line, the navigation's
// arithmetic, the Canvas2D and SVG lists the same lists, too many lines said.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import {
  Cut,
  type FrameLay,
  hitView,
  lineOf,
  lineStart,
  type NavContext,
  PlacePoint,
  recordDrawLists,
  renderSvg,
  type SectionLayout,
  TextNav,
  type TextPlacement,
} from '../index.js';
import { P53, p53TextSpec } from './fixtures.js';
import { settle } from './helpers.js';
import { rng } from './small-helpers.js';

const Size = { width: 720, height: 400 };

function text(s: SectionLayout): TextPlacement {
  return s.placement as TextPlacement;
}

describe('the text layout', () => {
  it('lines of whole groups of ten, a number a line, a band collapsed where no layer drew', async () => {
    const s = await settle(p53TextSpec());
    const f = await s.frame({ size: Size });
    const sec = f.layout.sections[0];
    const p = text(sec);
    const g = p.geometry!;
    expect(g.unitsPerLine % 10).toBe(0);
    expect(g.lines).toBe(Math.ceil(393 / g.unitsPerLine));
    const strips = p.strips;
    // each line's transform holds its units; the lines run in reading order, each a line below the one before
    for (let i = 0; i < strips.n; i++) {
      expect(strips.axes[i].u0).toBe(i * g.unitsPerLine);
      if (i) expect(strips.top[i]).toBeGreaterThan(strips.top[i - 1]);
    }
    // the first line's number is 1, the second's the first position on it
    const starts = p.marks.filter((m) => m.anchor === 2).map((m) => m.text);
    expect(starts.slice(0, 2)).toEqual(['1', String(g.unitsPerLine + 1)]);
    // the domains' band is used where a domain is (positions 0–61 and 94–292) and collapsed past 292
    const lastLine = Math.floor(392 / g.unitsPerLine);
    const firstFree = Math.floor(292 / g.unitsPerLine) + 1;
    if (firstFree <= lastLine && firstFree < strips.n) expect(strips.bandH[firstFree * 2 + 1]).toBe(0);
    expect(strips.bandH[1]).toBeGreaterThan(0);
  });

  it('placement is the inverse of hit-testing over random widths and scales; offsetAt runs in reading order', async () => {
    const s = await settle(p53TextSpec());
    const r = rng(3);
    for (let n = 0; n < 6; n++) {
      const width = 360 + Math.floor(r() * 600);
      const scale = 11 + Math.floor(r() * 8);
      const f = await s.frame({ size: { width, height: 2_000 }, viewport: { start: 0, end: 393, scale } });
      const p = text(f.layout.sections[0]);
      const pt = new PlacePoint();
      let last = -1;
      const top = p.strips.top;
      for (let u = 0; u < p.window.end; u++) {
        expect(p.place(u + 0.5, pt)).toBe(true);
        const y = top[pt.strip] + p.strips.bandY[pt.strip * 2] + 2;
        // the last line in view may run past the box's bottom (it scrolls inside)
        if (y >= p.box.y + p.box.height) break;
        expect(p.at(pt.x, y, pt)).toBe(true);
        expect(pt.kind).toBe('data');
        expect(Math.floor(pt.offset)).toBe(u);
        const o = p.offsetAt(p.strips.axes[pt.strip].xOfUnit(u + 0.5), y);
        expect(o).toBeGreaterThan(last);
        last = o;
      }
    }
  });

  it('a selection across a line break is cut there: square at the break, round where it really ends', async () => {
    const s = await settle(p53TextSpec());
    const f0 = await s.frame({ size: Size });
    const upl = text(f0.layout.sections[0]).geometry!.unitsPerLine;
    const selection = rangesLocation('/data/p53', [upl - 5, upl + 5]);
    const f = await s.frame({ size: Size, emphasis: { selection } });
    const l = f.lists.passes.selection.rects;
    const tints = [];
    for (let i = 0; i < l.n; i++) if (l.fill[i] !== 0) tints.push({ y: l.y[i], round: l.round[i] });
    // one tint a line it touches: round at its start on the first line, at its end on the next, square at the break
    expect(tints.length).toBe(2);
    expect(tints[0].round).toBe(1);
    expect(tints[1].round).toBe(2);
    expect(tints[1].y).toBeGreaterThan(tints[0].y);
  });

  it('a reflow keeps the anchor’s line at the top; a narrower width, more lines', async () => {
    const s = await settle(p53TextSpec());
    // a placement is kept frame to frame (its buffers): read what each frame placed before the next
    const lines: number[] = [];
    for (const width of [900, 500]) {
      const f = await s.frame({ size: { width, height: 300 }, viewport: { start: 120, end: 393 } });
      const p = text(f.layout.sections[0]);
      lines.push(p.geometry!.lines);
      const first = p.strips.axes[0];
      expect(first.u0 <= 120 && 120 < first.u1).toBe(true);
    }
    expect(lines[1]).toBeGreaterThan(lines[0]);
  });

  it('an export draws every line (its height auto), the same lists the canvas draws', async () => {
    const r = await renderSvg(p53TextSpec(), { size: { width: 720, height: 'auto' } });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    const glyphs = recordDrawLists(r.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('glyphs'));
    // every letter of p53 drawn once, in order
    const letters = glyphs
      .map((l) =>
        l
          .replace(/^glyphs "[^"]*" [\d.]+ \S+ /, '')
          .split(', ')
          .map((x) => x[0])
          .join(''),
      )
      .join('');
    expect(letters).toBe(P53);
    expect(r.lists.height).toBeGreaterThan(300);
  });

  it('the numbers print a layer’s numbering (a ruler’s `numbering`), else ordinals', async () => {
    const b = createSeqSpec({ title: 'numbered' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
    b.views.view('main', { layout: 'text' }, (v) => {
      const sec = v.section('protein', { of: p });
      sec.track('letters').layer('letters', { type: 'letters' });
      sec.track('ruler').layer('ruler', { type: 'ruler', numbering: 'index' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: Size });
    const t = text(f.layout.sections[0]);
    const starts = t.marks.filter((m) => m.anchor === 2).map((m) => m.text);
    expect(starts.slice(0, 2)).toEqual(['0', String(t.geometry!.unitsPerLine)]);
  });

  it('`text.pieceBreak: line`: a piece edge starts a line; steps, focus and clamp follow the lines', async () => {
    const b = createSeqSpec({ title: 'pieces' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
    const d = b.data.ranges('d', {
      target: p,
      rows: { start: Int32Array.from([0, 94]), end: Int32Array.from([61, 292]), id: ['TAD', 'DBD'] },
    });
    b.views.view('main', { layout: 'text' }, (v) => {
      v.section('protein', { of: p, pieces: d }).track('letters').layer('letters', { type: 'letters' });
    });
    const s = await settle(b.build(), [], { text: { pieceBreak: 'line', maxHeight: 2000 } });
    const f = await s.frame({ size: Size });
    const t = text(f.layout.sections[0]);
    const g = t.geometry!;
    const upl = g.unitsPerLine;
    // the first piece (61 units) in lines of upl, the second from a line of its own
    const firsts = [];
    for (let i = 0; i < t.strips.n; i++) firsts.push(t.strips.axes[i].u0);
    const second = Math.ceil(61 / upl);
    expect(firsts[second]).toBe(61);
    expect(t.strips.axes[second - 1].u1).toBe(61);
    expect(g.lines).toBe(second + Math.ceil((292 - 94) / upl));
    // every unit is placed once, `at` its inverse
    const pt = new PlacePoint();
    for (let u = 0; u < 61 + 198; u += 7) {
      expect(t.place(u + 0.5, pt)).toBe(true);
      const page = { x: pt.x, y: pt.y };
      expect(t.at(page.x, page.y, pt)).toBe(true);
      expect(Math.floor(pt.offset)).toBe(u);
    }
    const c: NavContext = { total: 259, width: 600, height: 400, cyclic: false, maxScale: 32, placement: t };
    // down from the line above the short one: its last unit; then the next piece's line, the same column
    expect(TextNav.step(c, (second - 1) * upl - 2, 'down')).toBe(60);
    expect(TextNav.step(c, (second - 1) * upl + 0, 'down')).toBe(61);
    expect(TextNav.step(c, 61 + 3, 'up')).toBe((second - 1) * upl + 3);
    // a unit's line and its line's start, as a focus reads them
    expect(lineStart(g, lineOf(g, 70))).toBe(61);
    // groups start at each line's first unit: every line at the cell's scale, each group label over the letter before
    // its gap (the ordinal of the group's last unit)
    for (let i = 0; i < t.strips.n; i++) {
      const ax = t.strips.axes[i];
      if (ax.u1 - ax.u0 >= upl) expect(ax.scale).toBeCloseTo(g.cell, 9);
      const labels = t.marks.filter((m) => m.strip === i && m.anchor === 1);
      let k = 0;
      for (let e = ax.u0 + 10; e <= ax.u1; e += 10, k++) {
        ax.cell(ax.elementAtUnit(e - 1));
        expect(labels[k].text).toBe(String(ax.elementAtUnit(e - 1) + 1));
        expect(labels[k].x).toBeCloseTo((ax.sx0 + ax.sx1) / 2, 6);
      }
    }
  });

  it('frames across tracks: an outline over the named tracks’ bands, open where a line break cuts it, titled where it starts', async () => {
    const b = createSeqSpec({ title: 'frames' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
    const regions = b.data.ranges('regions', {
      target: p,
      rows: { start: Int32Array.from([95]), end: Int32Array.from([125]), id: ['A'], label: ['Frame A'] },
    });
    b.views.view('main', { layout: 'text' }, (v) => {
      const sec = v.section('protein', { of: p, group: 10 });
      sec.track('frames').layer('f', { type: 'frame', data: regions, across: ['letters', 'ruler'] });
      sec.track('letters').layer('letters', { type: 'letters' });
      sec.track('ruler').layer('ruler', { type: 'ruler' });
    });
    const s = await settle(b.build(), [], { text: { maxHeight: 4000 } });
    const f = await s.frame({ size: { width: 520, height: 400 } });
    const sec = f.layout.sections[0];
    const t = text(sec);
    const upl = t.geometry!.unitsPerLine;
    const first = Math.floor(95 / upl);
    const last = Math.floor(124 / upl);
    expect(last).toBeGreaterThan(first);
    const nb = t.strips.bands;
    const frameOn = (strip: number) => sec.lays[strip * nb][0].lay as FrameLay;
    const a = frameOn(first);
    const z = frameOn(last);
    // one mark a line: open at the break (its right end on the first line, its left on the last), titled where it starts
    expect(a.m).toBe(1);
    expect(a.marks[2] & Cut.right).toBe(Cut.right);
    expect(a.marks[2] & Cut.left).toBe(0);
    expect(a.titled[0]).toBe(1);
    expect(z.marks[2] & Cut.left).toBe(Cut.left);
    expect(z.titled[0]).toBe(0);
    // over the letters and the ruler: from the letters' band's top to the ruler's bottom on its strip
    const k = first * nb;
    expect(a.span.y).toBe(t.strips.bandY[k + 1]);
    expect(a.span.y + a.span.height).toBe(t.strips.bandY[k + 2] + t.strips.bandH[k + 2]);
    // a line without the region has no outline, its title band collapsed there
    const away = Math.min(t.strips.n - 1, last + 2);
    expect(t.strips.bandH[away * nb]).toBe(0);
    // its outline hits as the region's item
    const x = (a.marks[0] + a.marks[1]) / 2;
    const h = hitView(f.layout, x, t.strips.dy[first] + a.span.y);
    expect(h?.element).toEqual({ owner: '/data/regions', key: 'A' });
    const record = recordDrawLists(f.lists, ['labels']);
    expect(record).toContain('Frame A');
  });

  it('too many lines are said, and the lines past the bound are not laid out', async () => {
    const s = await settle(p53TextSpec(), [], { text: { maxLines: 2 } });
    const f = await s.frame({ size: Size });
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.text.lines');
    expect(text(f.layout.sections[0]).strips.n).toBeLessThanOrEqual(2);
  });
});

describe('the text layout’s navigation', () => {
  it('a pan scrolls lines, a step goes a line up and down, a reveal brings a line to the top, clamp holds the end', async () => {
    // a box of 120 px: the text scrolls inside it
    const s = await settle(p53TextSpec(), [], { text: { maxHeight: 120 } });
    const f = await s.frame({ size: { width: 720, height: 400 }, viewport: { start: 0, end: 393 } });
    const p = text(f.layout.sections[0]);
    const g = p.geometry!;
    const c: NavContext = {
      total: 393,
      width: 600,
      height: 120,
      cyclic: false,
      maxScale: 32,
      minScale: 10,
      placement: p,
    };
    const pitch = g.lineHeight + g.lineGap;
    const v = { start: 0, end: 393 };
    expect(TextNav.pan(c, v, 0, pitch).start).toBeCloseTo(g.unitsPerLine);
    expect(TextNav.pan(c, v, 0, pitch / 2).start).toBeCloseTo(g.unitsPerLine / 2);
    expect(TextNav.step(c, 5, 'down')).toBe(5 + g.unitsPerLine);
    expect(TextNav.step(c, 5, 'up')).toBe(-1);
    expect(TextNav.step(c, 392, 'down')).toBe(-1);
    expect(TextNav.reveal(c, v, 2, 4)).toBe(v);
    expect(TextNav.reveal(c, v, 380, 390).start).toBeLessThanOrEqual(380);
    // never past where the last line meets the box's bottom
    const end = TextNav.clamp(c, { start: 393, end: 393 });
    expect(end.start).toBeLessThan(393);
    expect(TextNav.zoom(c, { ...v, scale: 12 }, 2, 100).scale).toBe(24);
    expect(TextNav.zoom(c, { ...v, scale: 12 }, 0.1, 100).scale).toBe(10);
    expect(TextNav.fit(c, { start: 100, end: 200, scale: 20 })).toEqual({ start: 0, end: 100 });
  });

  it('a pan reads the viewport it is given (two pans before a layout add up); a zoom keeps the pointer’s line at its y', async () => {
    const s = await settle(p53TextSpec(), [], { text: { maxHeight: 120 } });
    const f = await s.frame({ size: { width: 720, height: 400 }, viewport: { start: 0, end: 393 } });
    const p = text(f.layout.sections[0]);
    const g = p.geometry!;
    const c: NavContext = {
      total: 393,
      width: 600,
      height: 120,
      cyclic: false,
      maxScale: 32,
      minScale: 10,
      placement: p,
    };
    const pitch = g.lineHeight + g.lineGap;
    const once = TextNav.pan(c, { start: 0, end: 393 }, 0, pitch);
    expect(TextNav.pan(c, once, 0, pitch).start).toBeCloseTo(2 * g.unitsPerLine);
    // the scroll's sign: a positive dy shows what is below
    expect(TextNav.pan(c, once, 0, -pitch).start).toBeCloseTo(0);
    // zoomed in about a unit two lines below the anchor: it stays two lines below, on the new breaking
    const at = 2 * g.unitsPerLine + 3;
    const z = TextNav.zoom(c, { start: 0, end: 393, scale: g.cell }, 1.5, at);
    const laid = await s.frame({
      size: { width: 720, height: 400 },
      viewport: { start: z.start, end: 393, scale: z.scale },
    });
    const q = text(laid.layout.sections[0]);
    const lineOf = (u: number, upl: number) => Math.floor(u / upl);
    expect(lineOf(at, q.geometry!.unitsPerLine) - lineOf(q.anchor, q.geometry!.unitsPerLine)).toBe(2);
    expect(q.geometry!.unitsPerLine).toBeLessThan(g.unitsPerLine);
  });
});
