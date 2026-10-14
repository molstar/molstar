/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Selecting in a wrapped section (seqspec:interaction#selection), through input records on the built-in text preset: a
// click on a residue selects it, again clears it (a double or triple click is a click), Shift extends from the anchor;
// a plain drag selects positions line by line, past a line's end, between
// groups and lines, the press's side by the drag's direction, inline labels, into another block's lines (clamped to its
// own), on cyclic and reversed blocks and in another text section; the scroll past a box section's edges (from the
// first drag event, per millisecond, stopping at the end, on blur, a cancel and the behavior's removal); Shift- and
// Alt-drag; an underline under the letters is its item, by lane, and the item drawn there for ranges, points and
// secondary; the keyboard; a double click a zoom takes, which puts back its first click's selection (linear); fast
// pointers and releases the page lost, linear drags too; and the drag drawn as the
// selection it becomes (seqspec:interaction#input-events for the releases).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { OrderedSet, OrderedSet as OS } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { Pass, SkinToken } from '../src/draw/index.js';
import {
  DefaultPresets,
  Hit,
  HoverKind,
  Inputs,
  KernelRenderers,
  type LayerRenderer,
  PointerType,
  type PresetEntry,
  type Selection,
  type SeqViewerSpec,
} from '../src/index.js';
import type { Layouts } from '../src/layout.js';
import { runX } from '../src/runs.js';
import { setPiece } from './axis-fixtures.js';
import { type Driven, driven } from './input-fixtures.js';
import { lit, passRects, seqLoci } from './overlay-fixtures.js';
import { allTypes } from './prepare-fixtures.js';

const S = 'v.main';
const WRAPPED = DefaultPresets.find((p) => p.name === 'wrapped') as PresetEntry;
const TEXT: PresetEntry = { ...WRAPPED, name: 'lab:text' };
const INLINE: PresetEntry = {
  ...WRAPPED,
  name: 'lab:inline',
  render: { ...WRAPPED.render, wrappedGroupLabels: 'inlineStart' },
};

/** 400 residues, a second block of 60, the letters and two overlapping features in two lanes. */
function textDoc() {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(40) });
  const c = b.data.sequence({ ref: 'c', value: 'MKVLATGHEK'.repeat(6) });
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: [
      { start: 2, stop: 30, label: 'one' },
      { start: 20, stop: 40, label: 'two' },
    ],
  });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ sequence: a });
  v.block({ sequence: c });
  v.track({ ref: 'v.letters' }).letters();
  v.track({ ref: 'v.features' }).ranges({ data: f });
  return b.getDocument();
}

interface TextOptions {
  height?: number;
  doc?: SeqViewSpecDoc;
  /** The sections styled as text (default `v.main`). */
  sections?: readonly string[];
  preset?: string;
}
async function text(height = 900, o: TextOptions = {}) {
  const sections = o.sections ?? [S];
  const t = await driven(o.doc ?? textDoc(), {
    size: Object.fromEntries(sections.map((s) => [s, [400, o.height ?? height] as const])),
    registry: [{ presets: [TEXT, INLINE] }],
    renderers: KernelRenderers as LayerRenderer[],
    labels: false,
  });
  for (const s of sections) t.ctx.ui.setStyle(s, o.preset ?? 'lab:text');
  await t.ctx.graph.update();
  const l = t.layout(sections[0]);
  /** Shown position k of block b: its line, the center x of its letter and the center y of `track`'s band. */
  const at = (b: number, k: number, track = 'v.letters', section = sections[0]) => {
    const ls = t.layout(section);
    let line = -1;
    for (let i = 0; i < ls.lines.count; i++) if (ls.lines.block[i] === b && ls.lines.start[i] <= k) line = i;
    const v = t.ctx.layout.axis(section);
    const run = v?.runs[b];
    if (!v || !run || line < 0) throw new Error(`no line for ${b}:${k}`);
    const x = ls.dataLeft + runX(run, ls.lines.start[line], k, v.gap) + run.advance / 2;
    const tr = t.trackIndex(track);
    for (let e = 0; e < ls.entryCount; e++)
      if (ls.entryTrack[e] === tr && ls.entryRow[e] === line) {
        const top = ls.scrolledTop + ls.entryTop[e] - t.vp(section).live.current.scrollY;
        return { x, y: top + ls.entryHeight[e] / 2, top, height: ls.entryHeight[e], line, advance: run.advance };
      }
    throw new Error(`no ${track} band on line ${line}`);
  };
  /** The line holding shown k of block b. */
  const lineOf = (b: number, k: number) => {
    let line = -1;
    for (let i = 0; i < l.lines.count; i++) if (l.lines.block[i] === b && l.lines.start[i] <= k) line = i;
    return line;
  };
  return { ...t, l, at, lineOf };
}

/** The selection's positions per unit, as [unit, indices]. */
const positions = (s: Selection) => {
  const out: [number, number[]][] = [];
  for (const e of s.seq?.elements ?? []) {
    const xs: number[] = [];
    OrderedSet.forEach(e.indices, (i: number) => xs.push(i), null);
    out.push([e.unit, xs]);
  }
  return out;
};
const range = (a: number, b: number) => Array.from({ length: b - a }, (_, i) => a + i);
const sel = (t: Driven) => t.ctx.selection.current.value;

/** A y in [y0, y1) of a section at x whose hit is item `item` (by the hit, so a line's lanes need not be known). */
function itemY(t: Driven, section: string, x: number, y0: number, y1: number, item: number) {
  const h = Hit.create();
  for (let y = y0; y < y1; y += 0.5)
    if (t.ctx.layout.hitTest(section, x, y, h) && h.kind === HoverKind.item && h.item === item) return y;
  throw new Error(`no item ${item} at ${x} in [${y0}, ${y1})`);
}
/** A press, then its moves (no release). */
function press(t: Driven, x: number, y: number, moves: [number, number][], opts: { shift?: boolean } = {}) {
  t.pointer('down', S, x, y, opts);
  for (const [mx, my] of moves) t.pointer('move', S, mx, my, { ...opts, onDocument: true });
}
/** A pointer record with no button down: a move after a release the page lost. */
function unpressedMove(t: Driven, x: number, y: number) {
  const p = Inputs.pointer();
  Object.assign(p, { type: PointerType.move, section: t.ctx.layout.sectionIndex(S), x, y, buttons: 0, button: -1 });
  t.ctx.input.pointer(p);
}

describe('selecting in text', () => {
  it('the text preset selects: a click selects, again clears, a double click keeps; Shift extends; nothing clears', async () => {
    const t = await text();
    const s = t.ctx.layout.sectionIndex(S);
    expect(t.ctx.behaviors.get('select')?.forSection(s)).not.toBeNull();
    // The text's inset from the canvas's edges: left, top, and below its last line.
    const last = t.l.entryCount - 1;
    expect([t.l.dataLeft, t.l.lines.top[0], t.l.scrolledHeight - t.l.entryTop[last] - t.l.entryHeight[last]]).toEqual([
      8, 8, 8,
    ]);
    const p = t.at(0, 15);
    t.pointer('move', S, p.x, p.y);
    expect(t.ctx.ui.cursor(S).value).toBe('crosshair');
    t.click(S, p.x, p.y);
    expect(positions(sel(t))).toEqual([[0, [15]]]);
    t.advance(500);
    t.click(S, p.x, p.y);
    expect(sel(t).seq).toBeNull();
    t.advance(500);
    t.click(S, p.x, p.y);
    t.advance(500);
    // Shift extends from the anchor (the last click's position), forward and back, replacing the selection.
    const k = t.l.lines.start[1] + 2;
    const q = t.at(0, k);
    t.click(S, q.x, q.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(15, k + 1)]]);
    t.advance(500);
    // Before the selection: from the click to its far end.
    const r = t.at(0, 3);
    t.click(S, r.x, r.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(3, k + 1)]]);
    t.advance(500);
    // Between two groups: no position, so nothing, and the selection clears.
    const g = t.at(0, 9);
    t.click(S, g.x + g.advance / 2 + 1, g.y);
    expect(sel(t).seq).toBeNull();
    expect(sel(t).origin).toBe('pointer');
  });

  it('a double and a triple click are a click: the position stays selected, never cleared', async () => {
    const t = await text();
    const p = t.at(0, 23);
    t.click(S, p.x, p.y);
    t.click(S, p.x, p.y);
    expect(positions(sel(t))).toEqual([[0, [23]]]);
    t.click(S, p.x, p.y);
    expect(positions(sel(t))).toEqual([[0, [23]]]);
    t.advance(500);
    // A click again later clears it, as a single click does.
    t.click(S, p.x, p.y);
    expect(sel(t).seq).toBeNull();
  });

  for (const preset of ['lab:text', 'lab:inline'] as const)
    it(`Shift-click after a drag across lines (${preset}): before it, after it and inside it, both ends exact`, async () => {
      const t = await text(900, { preset });
      const L = t.l.lines;
      const [s, e] = [30, L.start[2] + 5];
      const drag = (from: number, to: number) => {
        const [a, z] = [t.at(0, from), t.at(0, to)];
        t.drag(S, a.x, a.y, z.x, z.y);
        t.advance(500);
      };
      const shift = (k: number) => {
        const p = t.at(0, k);
        t.click(S, p.x, p.y, { shift: true });
        t.advance(500);
        return positions(sel(t));
      };
      drag(s, e);
      // Before: from the click to the old end; after: from the old start to the click.
      expect(shift(10)).toEqual([[0, range(10, e + 1)]]);
      drag(s, e);
      expect(shift(e + 20)).toEqual([[0, range(s, e + 21)]]);
      // Inside: the anchor (the drag's start) stays, the other end moves; dragged backward, the anchor is the end.
      drag(s, e);
      expect(shift(40)).toEqual([[0, range(s, 41)]]);
      drag(e, s);
      expect(shift(40)).toEqual([[0, range(40, e + 1)]]);
    });

  it("Shift-click on an item's selection (an underline clicked): its far end stays exactly the item's", async () => {
    const t = await text();
    // Item "two" covers [20, 40): a point on its underline at 33, found by the hit.
    const u = t.at(0, 33, 'v.features');
    const y = itemY(t, S, u.x, u.top, u.top + u.height, 1);
    t.click(S, u.x, y);
    expect(sel(t).items).not.toBeNull();
    t.advance(500);
    const before = t.at(0, 8);
    t.click(S, before.x, before.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(8, 40)]]);
    t.advance(500);
    t.click(S, u.x, y);
    t.advance(500);
    const after = t.at(0, 52);
    t.click(S, after.x, after.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(20, 53)]]);
  });

  it('Shift-click in text extends a selection made in another view (an item, a drag): from the model, ends exact', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(40) });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: a,
      value: [
        { start: 2, stop: 30, label: 'one' },
        { start: 20, stop: 40, label: 'two' },
      ],
    });
    const tr = b.view({ ref: 't', label: 'Tracks' });
    tr.block({ sequence: a });
    tr.track({ ref: 't.letters' }).letters();
    tr.track({ ref: 't.features' }).ranges({ data: f });
    const v = b.view({ ref: 'v', label: 'Text', layout: 'text', height: 'box' });
    v.block({ sequence: a });
    v.track({ ref: 'v.letters' }).letters();
    v.track({ ref: 'v.features' }).ranges({ data: f });
    const T = 't.main';
    const t = await text(900, { doc: b.getDocument() });
    const show = async (view: string) => {
      t.ctx.views.setCurrent(view);
      await t.ctx.graph.update();
      for (const ref of t.ctx.viewports.sections.value) t.vp(ref).setSize(ref === T ? 1000 : 400, 900, 1);
      await t.ctx.graph.update();
    };
    // In Tracks: a click on item "two" where it alone lies (35), at a y its bar hits.
    await show('t');
    const lt = t.layout(T);
    const h = Hit.create();
    const x35 = t.xOf(T, 0, 35);
    let y = -1;
    for (let yy = lt.scrolledTop; yy < lt.height && y < 0; yy += 0.5)
      if (t.ctx.layout.hitTest(T, x35, yy, h) && h.kind === HoverKind.item && h.item === 1) y = yy;
    expect(y).toBeGreaterThan(0);
    t.click(T, x35, y);
    // Its double click zooms there: the click's selection waits for the double-click time.
    expect(sel(t).items).toBeNull();
    t.advance(500);
    expect(sel(t).items).not.toBeNull();
    // In Text: Shift-click before it, after it, inside it.
    await show('v');
    const shift = (k: number) => {
      const p = t.at(0, k);
      t.click(S, p.x, p.y, { shift: true });
      t.advance(500);
      return positions(sel(t));
    };
    expect(shift(8)).toEqual([[0, range(8, 40)]]);
    await show('t');
    t.click(T, x35, y);
    t.advance(500);
    await show('v');
    expect(shift(52)).toEqual([[0, range(20, 53)]]);
    await show('t');
    t.click(T, x35, y);
    t.advance(500);
    await show('v');
    // Inside, with the anchor (35) at neither end: the far end stays.
    expect(shift(25)).toEqual([[0, range(25, 40)]]);
    // A drag in Tracks (Alt: positions), then Shift-click in Text inside it and before it.
    await show('t');
    const ly = t.yOf(T, 't.letters');
    t.drag(T, t.xOf(T, 0, 10), ly, t.xOf(T, 0, 30), ly, { alt: true });
    expect(positions(sel(t))).toEqual([[0, range(10, 31)]]);
    t.advance(500);
    await show('v');
    expect(shift(20)).toEqual([[0, range(10, 21)]]);
    expect(shift(5)).toEqual([[0, range(5, 21)]]);
  });

  it('Shift-click on a cyclic block across its origin extends the run as shown; after a Shift-drag added a run, from the anchor', async () => {
    const b = createSeqViewSpec();
    const c = b.data.sequence({ ref: 's', value: 'ACGTACGTTA'.repeat(4), alphabet: 'dna', cyclic: true });
    b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' })
      .block({ sequence: c, origin: 2 })
      .track({ ref: 'v.letters' })
      .letters();
    const t = await text(900, { doc: b.getDocument() });
    const [a, z] = [t.at(0, 36), t.at(0, 39)];
    t.drag(S, a.x, a.y, z.x, z.y);
    expect(positions(sel(t))).toEqual([[0, [0, 1, 38, 39]]]);
    t.advance(500);
    // Shown 30 (index 32), before the run shown 36..39: from it to the run's far end (shown 39, index 1).
    const p = t.at(0, 30);
    t.click(S, p.x, p.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, [0, 1, ...range(32, 40)]]]);

    const u = await text();
    u.ctx.behaviors.get('select')?.setParams({ shiftPositions: 'toggle' });
    const at = (k: number) => u.at(0, k);
    u.drag(S, at(2).x, at(2).y, at(4).x, at(4).y);
    u.advance(500);
    u.drag(S, at(30).x, at(30).y, at(33).x, at(33).y, { shift: true });
    u.advance(500);
    expect(positions(sel(u))).toEqual([[0, [2, 3, 4, 30, 31, 32, 33]]]);
    // Two runs, the anchor on the block (the added drag's start, 30): extending from it, not across the gap.
    u.ctx.behaviors.get('select')?.setParams({ shiftPositions: 'extend' });
    u.click(S, at(40).x, at(40).y, { shift: true });
    expect(positions(sel(u))).toEqual([[0, range(30, 41)]]);
  });

  it('linear: Shift-click on positions extends from the anchor; on an item adds or removes it; Ctrl/⌘ adds or removes', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(30) });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: a,
      value: [
        { start: 40, stop: 60, label: 'D' },
        { start: 100, stop: 120, label: 'E' },
      ],
    });
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ sequence: a });
    v.track({ ref: 'v.letters' }).letters();
    v.track({ ref: 'v.f' }).ranges({ data: f });
    const L = 'v.main';
    const t = await text(300, { doc: b.getDocument(), sections: [L] });
    t.ctx.ui.setStyle(L, 'tracks');
    await t.ctx.graph.update();
    t.vp(L).zoomTo(0, 0, 150);
    const ly = t.yOf(L, 'v.letters');
    const fy = t.yOf(L, 'v.f');
    const x = (i: number) => t.xOf(L, 0, i);
    const click = (i: number, y: number, o: { shift?: boolean; ctrl?: boolean } = {}) => {
      t.click(L, x(i), y, o);
      t.advance(500);
    };
    click(10, ly);
    expect(positions(sel(t))).toEqual([[0, [10]]]);
    click(20, ly, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(10, 21)]]);
    click(5, ly, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(5, 21)]]);
    click(25, ly, { ctrl: true });
    expect(positions(sel(t))).toEqual([[0, [...range(5, 21), 25]]]);
    click(25, ly, { ctrl: true });
    expect(positions(sel(t))).toEqual([[0, range(5, 21)]]);
    // Items: Shift adds or removes the clicked one.
    const items = () => {
      const it = sel(t).items;
      const out: number[] = [];
      if (it) for (const e of it.elements) OrderedSet.forEach(e.items, (i: number) => out.push(i), null);
      return out;
    };
    click(50, fy);
    click(110, fy, { shift: true });
    expect(items()).toEqual([0, 1]);
    click(50, fy, { shift: true });
    expect(items()).toEqual([1]);
  });

  it("Shift-click extends from a drag's start; on another block, from that block's end nearest the anchor", async () => {
    const t = await text();
    const a = t.at(0, 5);
    const b = t.at(0, 12);
    t.drag(S, a.x, a.y, b.x, b.y);
    t.advance(500);
    const c = t.at(0, 40);
    t.click(S, c.x, c.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(5, 41)]]);
    t.advance(500);
    const d = t.at(1, 3);
    t.click(S, d.x, d.y, { shift: true });
    expect(positions(sel(t))).toEqual([[1, range(0, 4)]]);
  });

  it('a plain drag selects positions line by line, and back; never pans', async () => {
    const t = await text();
    const L = t.l.lines;
    const a = t.at(0, 5);
    const pend: number[] = [];
    t.ctx.selection.pending.version.subscribe(() => pend.push(t.ctx.selection.pending.current.endShown));
    const z = t.at(0, L.start[2] + 3);
    const sx = t.vp(S).live.current.scrollX;
    const sy = t.vp(S).live.current.scrollY;
    t.drag(S, a.x, a.y, z.x, z.y, {}, 6);
    expect(t.ctx.ui.drag.value).toBeNull();
    expect(positions(sel(t))).toEqual([[0, range(5, L.start[2] + 4)]]);
    // The end followed the pointer onto each line it crossed (the press's, then each move's), then the release.
    expect(pend.length).toBeGreaterThan(3);
    expect(Math.max(...pend)).toBe(L.start[2] + 3);
    expect([t.vp(S).live.current.scrollX, t.vp(S).live.current.scrollY]).toEqual([sx, sy]);
    // Backward, from a later line to an earlier one.
    t.drag(S, z.x, z.y, a.x, a.y);
    expect(positions(sel(t))).toEqual([[0, range(5, L.start[2] + 4)]]);
  });

  it("past a line's end, between groups, between lines, above and below the section (by x), into another block", async () => {
    const t = await text();
    const L = t.l.lines;
    const a = t.at(0, 5);
    // Past the end of line 1 (in the inset): its last position.
    t.drag(S, a.x, a.y, t.l.dataRight + 3, t.at(0, L.start[1]).y);
    expect(positions(sel(t))).toEqual([[0, range(5, L.stop[1])]]);
    // In the gap after position 19, forward: 19.
    const g = t.at(0, 19);
    t.drag(S, a.x, a.y, g.x + g.advance / 2 + 1, g.y);
    expect(positions(sel(t))).toEqual([[0, range(5, 20)]]);
    // On the labels band above line 2 (its group numbers): line 2's position there.
    const two = t.at(0, L.start[2] + 4);
    t.drag(S, a.x, a.y, two.x, two.top - 3);
    expect(positions(sel(t))).toEqual([[0, range(5, L.start[2] + 5)]]);
    // A press on the labels band starts there too.
    t.drag(S, two.x, two.top - 3, a.x, a.y);
    expect(positions(sel(t))).toEqual([[0, range(5, L.start[2] + 5)]]);
    // Onto the second block's lines: to the first block's last position.
    const c = t.at(1, 3);
    t.drag(S, a.x, a.y, c.x, c.y);
    expect(positions(sel(t))).toEqual([[0, range(5, 400)]]);
    // Above the section: the point on its top edge, the first line, by x (position 2 there).
    const two2 = t.at(0, 2);
    t.drag(S, a.x, a.y, two2.x, -40, { onDocument: true });
    expect(positions(sel(t))).toEqual([[0, range(2, 6)]]);
    // Starting in the second block, onto the first block's lines: its first position; below the section (its bottom
    // edge, under the last line here): the last line by x.
    t.drag(S, c.x, c.y, a.x, a.y);
    expect(positions(sel(t))).toEqual([[1, range(0, 4)]]);
    const c40 = t.at(1, 40);
    t.drag(S, c.x, c.y, c40.x, t.l.height + 40, { onDocument: true });
    expect(positions(sel(t))).toEqual([[1, range(3, 41)]]);
  });

  it("a press between positions starts on the drag's side: forward the next, backward the one before", async () => {
    const t = await text();
    const L = t.l.lines;
    const g = t.at(0, 19);
    const gap = g.x + g.advance / 2 + 1;
    const f = t.at(0, 25);
    const b = t.at(0, 12);
    t.drag(S, gap, g.y, f.x, f.y);
    expect(positions(sel(t))).toEqual([[0, range(20, 26)]]);
    t.drag(S, gap, g.y, b.x, b.y);
    expect(positions(sel(t))).toEqual([[0, range(12, 20)]]);
    // Past a line's end, within the text: forward from the next line's first position, backward from the line's last.
    let line = -1;
    for (let i = 0; i + 1 < L.count && line < 0; i++) {
      const e = t.at(0, L.stop[i] - 1);
      if (L.block[i + 1] === 0 && t.l.dataRight - (e.x + e.advance / 2) > 4) line = i;
    }
    expect(line).toBeGreaterThanOrEqual(0);
    const end = t.at(0, L.stop[line] - 1);
    const past = end.x + end.advance / 2 + 2;
    const n1 = t.at(0, L.start[line + 1] + 3);
    t.drag(S, past, end.y, n1.x, n1.y);
    expect(positions(sel(t))).toEqual([[0, range(L.start[line + 1], L.start[line + 1] + 4)]]);
    const back = t.at(0, L.start[line] + 3);
    t.drag(S, past, end.y, back.x, back.y);
    expect(positions(sel(t))).toEqual([[0, range(L.start[line] + 3, L.stop[line])]]);
    // In the text's inset or its right strip a drag scrolls the rows (a touch reader's grip), never selects.
    const before = sel(t);
    t.drag(S, t.l.dataRight + 6, end.y, t.l.dataRight + 6, end.y - 40);
    t.drag(S, 3, n1.y, n1.x, n1.y);
    expect(sel(t)).toBe(before);
    expect(t.l.width - t.l.dataRight).toBeGreaterThanOrEqual(24);
    // The end between positions, backward: the one after the gap.
    const h = t.at(0, 9);
    t.drag(S, f.x, f.y, h.x + h.advance / 2 + 1, h.y);
    expect(positions(sel(t))).toEqual([[0, range(10, 26)]]);
  });

  it("inline group labels: a label belongs between its group and the one before; a press there starts on the drag's side", async () => {
    const t = await text(900, { preset: 'lab:inline' });
    const p = t.at(0, 20);
    const label = p.x - p.advance / 2 - 3;
    const f = t.at(0, 25);
    const b = t.at(0, 14);
    t.drag(S, label, p.y, f.x, f.y);
    expect(positions(sel(t))).toEqual([[0, range(20, 26)]]);
    t.drag(S, label, p.y, b.x, b.y);
    expect(positions(sel(t))).toEqual([[0, range(14, 20)]]);
    t.advance(500);
    // A click on the label: no position.
    t.click(S, label, p.y);
    expect(sel(t).seq).toBeNull();
  });

  it('a cyclic block across its origin, a reversed piece: the block indices of the shown range', async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 's', value: 'ACGTACGTTA'.repeat(4), alphabet: 'dna', cyclic: true });
    b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' })
      .block({ sequence: s, origin: 2 })
      .track({ ref: 'v.letters' })
      .letters();
    const t = await text(900, { doc: b.getDocument() });
    const a = t.at(0, 36);
    const z = t.at(0, 39);
    t.drag(S, a.x, a.y, z.x, z.y);
    expect(positions(sel(t))).toEqual([[0, [0, 1, 38, 39]]]);

    const r = createSeqViewSpec();
    const q = r.data.sequence({ ref: 'q', value: 'MKVLAAGIVGLLLAQRSTVW' });
    r.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' })
      .block({ sequence: q })
      .track({ ref: 'v.letters' })
      .letters();
    const u = await text(900, { doc: r.getDocument() });
    const m = u.ctx.graph.assembled.value?.view;
    if (!m) throw new Error('no view');
    setPiece(m.resolved, 0, 0, 20, -1, 20);
    (u.layouts as Layouts).update(m);
    const c = u.at(0, 2);
    const d = u.at(0, 5);
    u.drag(S, c.x, c.y, d.x, d.y);
    expect(positions(sel(u))).toEqual([[0, [14, 15, 16, 17]]]);
  });

  it('several text sections: each selects its own; a drag from one over another stays on its own block', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(20) });
    const c = b.data.sequence({ ref: 'c', value: 'MKVLATGHEK'.repeat(20) });
    const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
    v.section({ ref: 'v.one', height: 'box' }).block({ sequence: a }).track({ ref: 'v.letters' }).letters();
    v.section({ ref: 'v.two', height: 'box' }).block({ sequence: c }).track({ ref: 'v.letters2' }).letters();
    const t = await text(300, { doc: b.getDocument(), sections: ['v.one', 'v.two'] });
    const two = (k: number) => t.at(0, k, 'v.letters2', 'v.two');
    const [p, q] = [two(4), two(30)];
    t.drag('v.two', p.x, p.y, q.x, q.y);
    const unitOf = (section: string) => t.vp(section).axis?.unit[0] ?? -1;
    expect(positions(sel(t))).toEqual([[unitOf('v.two'), range(4, 31)]]);
    // From the first section down past its bottom (the press's surface keeps the drag, in its own px).
    const o = t.at(0, 10, 'v.letters', 'v.one');
    t.pointer('down', 'v.one', o.x, o.y);
    t.pointer('move', 'v.one', o.x, o.y + 10);
    t.pointer('up', 'v.one', o.x, t.layout('v.one').height + 80);
    const [[unit, xs]] = positions(sel(t));
    expect([unit, xs[0]]).toEqual([unitOf('v.one'), 10]);
  });
});

describe("the scroll past the rows' edges", () => {
  it('each frame while the drag lasts, the range following; back inside it stops; above, up again', async () => {
    const t = await text(160);
    const vp = t.vp(S);
    expect(t.l.vscroll).toBe(true);
    const a = t.at(0, 5);
    press(t, a.x, a.y, [
      [a.x, a.y + 10],
      [a.x, t.l.height + 30],
    ]);
    expect(t.ctx.ui.drag.value?.owner).toBe('select');
    const y0 = vp.live.current.scrollY;
    const end0 = t.ctx.selection.pending.current.endShown;
    for (let i = 0; i < 20; i++) t.advance(16);
    const y1 = vp.live.current.scrollY;
    const end1 = t.ctx.selection.pending.current.endShown;
    expect(y1).toBeGreaterThan(y0 + 100);
    expect(end1).toBeGreaterThan(end0 + 20);
    expect(vp.live.current.interacting).toBe(true);
    // Back inside: the scroll stops.
    t.pointer('move', S, a.x, t.l.scrolledTop + 40);
    t.advance(16);
    const y2 = vp.live.current.scrollY;
    t.advance(16);
    expect([vp.live.current.scrollY, t.scheduler.pendingFrames]).toEqual([y2, 0]);
    // Above the top: up again.
    t.pointer('move', S, a.x, -20, { onDocument: true });
    for (let i = 0; i < 5; i++) t.advance(16);
    expect(vp.live.current.scrollY).toBeLessThan(y2);
    t.pointer('up', S, a.x, -20, { onDocument: true });
    expect(vp.live.current.interacting).toBe(false);
    expect(t.scheduler.pendingFrames).toBe(0);
    // From the press to where the scroll left the point: one range.
    const [[unit, xs]] = positions(sel(t));
    expect([unit, xs[0], xs.length]).toEqual([0, 5, (xs.at(-1) ?? 0) - 4]);
    expect(xs.at(-1)).toBeGreaterThan(t.l.lines.stop[0]);
  });

  it('starts at the first drag event when it is already past the edge', async () => {
    const t = await text(160);
    const a = t.at(0, 5);
    press(t, a.x, a.y, [[a.x, t.l.height + 30]]);
    expect(t.ctx.ui.drag.value?.owner).toBe('select');
    const y0 = t.vp(S).live.current.scrollY;
    for (let i = 0; i < 10; i++) t.advance(16);
    expect(t.vp(S).live.current.scrollY).toBeGreaterThan(y0 + 40);
    t.pointer('up', S, a.x, t.l.height + 30, { onDocument: true });
  });

  it('goes as far in the same time at 60 and 120 frames a second', async () => {
    const run = async (ms: number, frames: number) => {
      const t = await text(160);
      const a = t.at(0, 5);
      press(t, a.x, a.y, [[a.x, t.l.height + 30]]);
      const y0 = t.vp(S).live.current.scrollY;
      for (let i = 0; i < frames; i++) t.advance(ms);
      const y = t.vp(S).live.current.scrollY - y0;
      t.pointer('up', S, a.x, t.l.height + 30, { onDocument: true });
      return y;
    };
    const [slow, fast] = [await run(16, 20), await run(8, 40)];
    expect(slow).toBeGreaterThan(60);
    expect(Math.abs(fast - slow)).toBeLessThan(8);
  });

  it("stops at the rows' end, on blur, on a cancel and when the behavior is removed; nothing pending after", async () => {
    const t = await text(160);
    const vp = t.vp(S);
    const a = t.at(0, 5);
    press(t, a.x, a.y, [[a.x, t.l.height + 300]]);
    for (let i = 0; i < 400 && t.scheduler.pendingFrames > 0; i++) t.advance(50);
    // At the end: no more frames, the range at the block's last position.
    expect(t.scheduler.pendingFrames).toBe(0);
    expect(vp.live.current.scrollY).toBeGreaterThan(0);
    const end = vp.live.current.scrollY;
    t.advance(50);
    expect(vp.live.current.scrollY).toBe(end);
    t.pointer('up', S, a.x, t.l.height + 300, { onDocument: true });
    expect(positions(sel(t))[0][1].at(-1)).toBe(399);

    // Blur pauses it; a move starts it again.
    vp.scrollTo(0, 0, 'api');
    press(t, a.x, a.y, [[a.x, t.l.height + 30]]);
    t.advance(16);
    t.focus(S, false);
    expect(t.scheduler.pendingFrames).toBe(0);
    const y = vp.live.current.scrollY;
    t.advance(16);
    expect(vp.live.current.scrollY).toBe(y);
    t.pointer('move', S, a.x, t.l.height + 31, { onDocument: true });
    t.advance(16);
    expect(vp.live.current.scrollY).toBeGreaterThan(y);

    // A cancel: nothing selected, nothing pending, the view committed.
    const before = sel(t);
    t.ctx.input.cancel();
    expect([t.scheduler.pendingFrames, t.ctx.selection.pending.current.active, vp.live.current.interacting]).toEqual([
      0,
      false,
      false,
    ]);
    expect(sel(t)).toBe(before);

    // The behavior removed mid-drag.
    press(t, a.x, a.y, [[a.x, t.l.height + 30]]);
    t.advance(16);
    t.ctx.behaviors.get('select')?.remove();
    expect([t.scheduler.pendingFrames, t.ctx.selection.pending.current.active]).toEqual([0, false]);
    t.pointer('up', S, a.x, t.l.height + 30, { onDocument: true });
  });
});

describe('modifiers, items and the keyboard', () => {
  it('Shift-drag adds the dragged positions (`shiftPositions: toggle`); Alt-drag selects as a plain drag', async () => {
    const t = await text();
    t.ctx.behaviors.get('select')?.setParams({ shiftPositions: 'toggle' });
    const a = t.at(0, 2);
    const b = t.at(0, 4);
    t.drag(S, a.x, a.y, b.x, b.y);
    const c = t.at(0, 30);
    const d = t.at(0, 33);
    t.drag(S, c.x, c.y, d.x, d.y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, [2, 3, 4, 30, 31, 32, 33]]]);
    t.drag(S, c.x, c.y, d.x, d.y, { alt: true });
    expect(positions(sel(t))).toEqual([[0, [30, 31, 32, 33]]]);
  });

  it('Shift with a click or a drag means the same: `extend` from the anchor; `toggle` adds, or removes from a selected press', async () => {
    const t = await text();
    const L = t.l.lines;
    const at = (k: number) => t.at(0, k);
    const drag = (a: number, z: number) => {
      t.drag(S, at(a).x, at(a).y, at(z).x, at(z).y, { shift: true });
      t.advance(500);
    };
    // Extend (the default): a click's anchor, then a Shift-drag from elsewhere extends from it to the drag's end.
    t.click(S, at(10).x, at(10).y);
    t.advance(500);
    drag(20, L.start[1] + 5);
    expect(positions(sel(t))).toEqual([[0, range(10, L.start[1] + 6)]]);
    drag(25, 4);
    expect(positions(sel(t))).toEqual([[0, range(4, 11)]]);
    // Toggle: a drag pressed on an unselected position adds, on a selected one removes, across lines; its preview says.
    t.ctx.behaviors.get('select')?.setParams({ shiftPositions: 'toggle' });
    t.ctx.selection.set(null, 'api');
    t.drag(S, at(5).x, at(5).y, at(L.start[2] + 5).x, at(L.start[2] + 5).y);
    t.advance(500);
    drag(L.start[2] + 20, L.start[2] + 22);
    expect(positions(sel(t))).toEqual([[0, [...range(5, L.start[2] + 6), ...range(L.start[2] + 20, L.start[2] + 23)]]]);
    const pd = t.ctx.selection.pending.current;
    t.pointer('down', S, at(12).x, at(12).y, { shift: true });
    t.pointer('move', S, at(L.start[1] + 3).x, at(L.start[1] + 3).y, { shift: true });
    expect([pd.active, pd.remove]).toEqual([true, true]);
    t.pointer('up', S, at(L.start[1] + 3).x, at(L.start[1] + 3).y, { shift: true });
    expect([pd.active, pd.remove]).toEqual([false, false]);
    expect(positions(sel(t))).toEqual([
      [0, [...range(5, 12), ...range(L.start[1] + 4, L.start[2] + 6), ...range(L.start[2] + 20, L.start[2] + 23)]],
    ]);
    // And a click: Shift adds or removes one position.
    t.advance(500);
    t.click(S, at(6).x, at(6).y, { shift: true });
    expect(positions(sel(t))[0][1].includes(6)).toBe(false);
  });

  it('an underline under the letters is its item, by lane; between lanes and under none, the position', async () => {
    const t = await text();
    const at = (k: number) => {
      const p = t.at(0, k, 'v.features');
      return { x: p.x, lane: (n: number) => p.top + p.height - 5 - n * 5 + 1.5, top: p.top };
    };
    const items = () => {
      const it = sel(t).items;
      const out: number[] = [];
      if (it) for (const e of it.elements) OrderedSet.forEach(e.items, (i: number) => out.push(i), null);
      return out;
    };
    // Position 25 is under both: "one" packed in lane 0, "two" in lane 1.
    const p = at(25);
    t.click(S, p.x, p.lane(0));
    expect([items(), sel(t).seq]).toEqual([[0], null]);
    t.advance(500);
    t.click(S, p.x, p.lane(1));
    expect(items()).toEqual([1]);
    t.advance(500);
    // On the band's top padding: the position.
    t.click(S, p.x, p.top + 1);
    expect([sel(t).items, positions(sel(t))]).toEqual([null, [[0, [25]]]]);
    t.advance(500);
    // Lane 1 at position 5: no feature there, the position.
    const q = at(5);
    t.click(S, q.x, q.lane(1));
    expect([sel(t).items, positions(sel(t))]).toEqual([null, [[0, [5]]]]);
  });

  it('ranges, points and secondary: every underline drawn hits the item the hover then outlines there', async () => {
    const t = await lit(allTypes({ layout: 'text' }), { width: 600, height: 4000 });
    const l = t.layout(t.section);
    const main = passRects(t.frame(), Pass.main);
    const failures: string[] = [];
    for (const track of ['Features', 'Sites', 'Secondary']) {
      const m = t.ctx.graph.assembled.value?.view?.resolved;
      const set = t.ctx.graph.assembled.value?.set;
      let tr = -1;
      for (let k = 0; m && set && k < m.tracks.count; k++) if (set.strings.get(m.tracks.label[k]) === track) tr = k;
      // The track's bands, as device rows.
      const bands: [number, number][] = [];
      for (let e = 0; e < l.entryCount; e++)
        if (l.entryTrack[e] === tr)
          bands.push([l.scrolledTop + l.entryTop[e], l.scrolledTop + l.entryTop[e] + l.entryHeight[e]]);
      const under = main.filter((r) => r[3] === 3 && bands.some(([y0, y1]) => r[1] >= y0 && r[1] + r[3] <= y1));
      let checked = 0;
      for (const r of under.slice(0, 40)) {
        const [x, y] = [r[0] + r[2] / 2, r[1] + r[3] / 2];
        t.hover(x, y);
        const h = t.ctx.hover.live.current;
        if (h.kind !== HoverKind.item) {
          failures.push(`${track} ${x},${y}: kind ${h.kind}`);
          continue;
        }
        const outlined = passRects(t.overlay(), Pass.overlay).some(
          (o) => o[0] === r[0] && o[1] === r[1] && o[2] === r[2] && o[3] === r[3],
        );
        if (!outlined) failures.push(`${track} ${x},${y}: item ${h.item} is not drawn there`);
        checked++;
      }
      expect([track, checked > 3]).toEqual([track, true]);
    }
    expect(failures).toEqual([]);
  });

  it("the keyboard: Enter selects the cursor's position, Shift+Enter extends to the next, Escape clears", async () => {
    const t = await text();
    t.focus(S);
    for (let i = 0; i < 3; i++) t.key(S, 'ArrowRight');
    const k = t.ctx.ui.keyCursor.value?.index ?? -1;
    t.key(S, 'Enter');
    expect(positions(sel(t))).toEqual([[0, [k]]]);
    for (let i = 0; i < 12; i++) t.key(S, 'ArrowRight');
    const k2 = t.ctx.ui.keyCursor.value?.index ?? -1;
    t.key(S, 'Enter', { shift: true });
    expect(positions(sel(t))).toEqual([[0, range(k, k2 + 1)]]);
    t.key(S, 'Escape');
    expect(sel(t).seq).toBeNull();
  });
});

describe('clicks where a double click means something else', () => {
  /** A linear section with item "D" [40, 90), and a text view of the same sequence; a host counting selection events. */
  async function both() {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(30) });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: a,
      value: [{ start: 40, stop: 90, label: 'D' }],
    });
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ sequence: a });
    v.track({ ref: 'v.f' }).ranges({ data: f });
    const w = b.view({ ref: 'w', label: 'Text', layout: 'text', height: 'box' });
    w.block({ sequence: a });
    w.track({ ref: 'v.letters' }).letters();
    const t = await text(400, { doc: b.getDocument(), sections: ['v.main', 'w.main'] });
    t.ctx.ui.setStyle('v.main', 'tracks');
    await t.ctx.graph.update();
    const events: string[] = [];
    t.ctx.selection.current.subscribe((s) => events.push(s.origin));
    events.length = 0;
    return { t, events };
  }

  it('a single click on an item selects once, after the double-click time; a double click zooms and never selects', async () => {
    const { t, events } = await both();
    const L = 'v.main';
    const vp = t.vp(L);
    const [x, y] = [t.xOf(L, 0, 60), t.yOf(L, 'v.f')];
    t.click(L, x, y);
    expect([events, sel(t).items]).toEqual([[], null]);
    t.advance(320);
    expect(events).toEqual(['pointer']);
    expect(sel(t).items).not.toBeNull();
    t.advance(500);
    // Cleared, then a double click: the zoom, and no selection event at all, ever.
    t.ctx.selection.set(null, 'api');
    events.length = 0;
    const cell = vp.live.current.cell;
    t.click(L, x, y);
    t.click(L, x, y);
    t.advance(1000);
    expect(vp.live.current.cell).toBeGreaterThan(cell);
    expect([events, sel(t).items]).toEqual([[], null]);
    expect(t.ctx.history.canUndo.value).toBe(false);
    // A held click, then Escape, or a host's selection: it never applies.
    t.click(L, x, y);
    t.key(L, 'Escape');
    t.advance(500);
    expect(events).toEqual([]);
    t.click(L, x, y);
    t.ctx.selection.set(seqLoci(t.ctx.graph.assembled.value?.set as never, 'a', 3, 4), 'host:x');
    t.advance(500);
    expect(events).toEqual(['host:x']);
  });

  it('a held click still applies after a patch carries the selection to a new set meanwhile (VS-06)', async () => {
    const { t, events } = await both();
    const L = 'v.main';
    t.ctx.selection.set(seqLoci(t.ctx.graph.assembled.value?.set as never, 'a', 3, 4), 'api');
    expect(events).toEqual(['api']);
    const before = sel(t).version;
    // A click held on the item, then a data patch (a new set: the selection is carried with origin `graph`).
    t.click(L, t.xOf(L, 0, 60), t.yOf(L, 'v.f'));
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'f', params: { label: 'Domains' } }],
    });
    expect(events).toEqual(['api', 'graph']);
    t.advance(500);
    expect(events).toEqual(['api', 'graph', 'pointer']);
    expect(sel(t).items).not.toBeNull();
    expect(sel(t).version).toBeGreaterThan(before);
    // The held click's loci were carried to the new set too (LB-06): the selection holds none of the old one's.
    expect(sel(t).items?.set).toBe(t.ctx.graph.assembled.value?.set);
  });

  it("in text a click selects at once (nothing takes its double click); one event a click, none for a double click's second", async () => {
    const { t, events } = await both();
    t.ctx.views.setCurrent('w');
    await t.ctx.graph.update();
    t.vp('w.main').setSize(400, 400, 1);
    await t.ctx.graph.update();
    const p = t.at(0, 15, 'v.letters', 'w.main');
    t.click('w.main', p.x, p.y);
    expect(events).toEqual(['pointer']);
    t.click('w.main', p.x, p.y);
    expect(events).toEqual(['pointer']);
    expect(positions(sel(t))).toEqual([[0, [15]]]);
  });
});

describe("the rules are the select behavior's params, in the viewer spec", () => {
  /** A linear section of letters and an item [40, 90); `behaviors` composes the spec's. */
  async function linear(behaviors?: (l: SeqViewerSpec['behaviors']) => SeqViewerSpec['behaviors']) {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(30) });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: a,
      value: [{ start: 40, stop: 90, label: 'D' }],
    });
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ sequence: a });
    v.track({ ref: 'v.letters' }).letters();
    v.track({ ref: 'v.f' }).ranges({ data: f });
    const L = 'v.main';
    const t = await driven(b.getDocument(), {
      size: { [L]: [1000, 200] },
      renderers: KernelRenderers as LayerRenderer[],
      labels: false,
      behaviors,
    });
    t.vp(L).zoomTo(0, 0, 150);
    return { t, L, x: (i: number) => t.xOf(L, 0, i), y: t.yOf(L, 'v.letters') };
  }

  it('a spec without the select behavior selects nothing on a click, a drag or a key', async () => {
    const { t, L, x, y } = await linear((l) => l.filter((b) => b.behavior.name !== 'select'));
    expect(t.ctx.behaviors.get('select')).toBeNull();
    t.click(L, x(10), y);
    t.advance(500);
    t.drag(L, x(10), y, x(20), y, { alt: true });
    t.focus(L);
    t.key(L, 'ArrowRight');
    t.key(L, 'Enter');
    expect(sel(t).seq).toBeNull();
  });

  it('`shiftPositions: toggle`, `doubleClickWait: false`, `autoscroll: false`, set over the defaults', async () => {
    const { t, L, x, y } = await linear();
    t.ctx.behaviors.get('select')?.setParams({ shiftPositions: 'toggle', doubleClickWait: false, autoscroll: false });
    // No wait: a click on data selects at once (its double click still zooms).
    t.click(L, x(10), y);
    expect(positions(sel(t))).toEqual([[0, [10]]]);
    t.advance(500);
    t.click(L, x(20), y, { shift: true });
    expect(positions(sel(t))).toEqual([[0, [10, 20]]]);
    t.advance(500);
    // No scroll past an edge: the view stays, the end at the edge's column.
    const vp = t.vp(L);
    const sx = vp.live.current.scrollX;
    t.pointer('down', L, x(30), y, { alt: true });
    t.pointer('move', L, x(35), y, { alt: true });
    t.pointer('move', L, t.layout(L).width + 40, y, { alt: true, onDocument: true });
    for (let i = 0; i < 10; i++) t.advance(16);
    expect(vp.live.current.scrollX).toBe(sx);
    t.pointer('up', L, t.layout(L).width + 40, y, { alt: true, onDocument: true });
  });
});

describe('a release never seen is no click', () => {
  it('a press whose release was lost, then a hover minutes later: nothing selected or cleared', async () => {
    const t = await text();
    const p = t.at(0, 15);
    t.pointer('down', S, p.x, p.y);
    t.advance(120_000);
    unpressedMove(t, p.x + 1, p.y);
    expect(sel(t).seq).toBeNull();
    // Selected, then the same: never cleared by it either.
    t.click(S, p.x, p.y);
    t.advance(500);
    t.pointer('down', S, p.x, p.y);
    unpressedMove(t, p.x, p.y);
    expect(positions(sel(t))).toEqual([[0, [15]]]);
  });

  it('a lost release, then a real click: one click (count 1), never a double click', async () => {
    const t = await text();
    const counts: number[] = [];
    t.ctx.interaction.click.subscribe((e) => counts.push(e.count));
    const p = t.at(0, 15);
    t.pointer('down', S, p.x, p.y);
    t.advance(100);
    t.click(S, p.x, p.y);
    expect(counts).toEqual([1]);
  });
});

describe('a pressed point outside its surface', () => {
  it("is clamped to the surface's bounds on every side and corner; its overshoot kept", async () => {
    const t = await text(300);
    const l = t.l;
    const seen: number[][] = [];
    t.ctx.interaction.pointer.subscribe((e) => seen.push([e.input.x, e.input.y, e.input.overX, e.input.overY]));
    const p = t.at(0, 15);
    t.pointer('down', S, p.x, p.y);
    const pts: [number, number][] = [
      [-50, p.y],
      [l.width + 30, p.y],
      [p.x, -20],
      [p.x, l.height + 40],
      [-5, -6],
      [l.width + 7, l.height + 8],
    ];
    seen.length = 0;
    for (const [x, y] of pts) t.pointer('move', S, x, y, { onDocument: true });
    t.pointer('up', S, -50, l.height + 40, { onDocument: true });
    expect(seen).toEqual([
      [0, p.y, -50, 0],
      [l.width, p.y, 30, 0],
      [p.x, 0, 0, -20],
      [p.x, l.height, 0, 40],
      [0, 0, -5, -6],
      [l.width, l.height, 7, 8],
      [0, l.height, -50, 40],
    ]);
    // Not pressed: as given, no overshoot.
    t.pointer('move', S, l.width + 30, p.y, { onDocument: true });
    expect(seen.at(-1)).toEqual([l.width + 30, p.y, 0, 0]);
  });

  it("text: out above or below, the end follows x along the edge's line; in a scrolled box, the line at the edge", async () => {
    const t = await text(160);
    const L = t.l.lines;
    const vp = t.vp(S);
    // Scrolled so line 2 is at the top.
    vp.scrollTo(0, L.top[2] - t.l.scrolledTop + 1, 'api');
    await t.ctx.graph.update();
    const a = t.at(0, L.start[3] + 2);
    t.pointer('down', S, a.x, a.y);
    t.pointer('move', S, a.x + 8, a.y);
    const pd = t.ctx.selection.pending.current;
    for (const k of [3, 7, 1]) {
      t.pointer('move', S, t.at(0, L.start[2] + k).x, -15, { onDocument: true });
      expect(pd.endShown).toBe(L.start[2] + k);
    }
    t.pointer('up', S, a.x, t.l.scrolledTop + 40);
  });

  it('linear: out above or below, the end follows x; past the left or right it pans and extends; back inside it adjusts', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(60) });
    const v = b.view({ ref: 'v', label: 'V' });
    v.block({ sequence: a });
    v.track({ ref: 'v.letters' }).letters();
    const L = 'v.main';
    const t = await text(200, { doc: b.getDocument(), sections: [L] });
    t.ctx.ui.setStyle(L, 'tracks');
    await t.ctx.graph.update();
    const vp = t.vp(L);
    vp.zoomTo(0, 100, 160);
    const y = t.yOf(L, 'v.letters');
    const x = (i: number) => t.xOf(L, 0, i);
    const pd = t.ctx.selection.pending.current;
    // Alt: a selection drag on data.
    t.pointer('down', L, x(110), y, { alt: true });
    t.pointer('move', L, x(112), y, { alt: true });
    for (const k of [120, 130, 105]) {
      t.pointer('move', L, x(k), -30, { alt: true, onDocument: true });
      expect(pd.endShown).toBe(vp.axis?.shownOf(0, k));
      t.pointer('move', L, x(k), t.layout(L).height + 30, { alt: true, onDocument: true });
      expect(pd.endShown).toBe(vp.axis?.shownOf(0, k));
    }
    // Past the right edge: the view pans right each frame and the end goes with it, past what was in view.
    const sx = vp.live.current.scrollX;
    t.pointer('move', L, t.layout(L).width + 40, y, { alt: true, onDocument: true });
    for (let i = 0; i < 20; i++) t.advance(16);
    expect(vp.live.current.scrollX).toBeGreaterThan(sx + 50);
    expect(pd.endShown).toBeGreaterThan(160);
    // Back inside: it adjusts to the column under the point, and the pan stops.
    const xs = vp.live.current.scrollX;
    const inside = (t.layout(L).dataLeft + t.layout(L).dataRight) / 2;
    t.pointer('move', L, inside, y, { alt: true });
    t.advance(16);
    t.advance(16);
    expect(vp.live.current.scrollX).toBe(xs);
    const h = Hit.create();
    t.ctx.layout.hitTest(L, inside, y, h);
    expect(pd.endShown).toBe(h.shown);
    t.pointer('up', L, inside, y, { alt: true });
    const [[, xs2]] = positions(sel(t));
    expect([xs2[0], xs2.at(-1)]).toEqual([110, h.index]);
    // Past the left edge: it pans left.
    t.advance(500);
    t.pointer('down', L, inside, y, { alt: true });
    t.pointer('move', L, inside - 30, y, { alt: true });
    const sl = vp.live.current.scrollX;
    t.pointer('move', L, -40, y, { alt: true, onDocument: true });
    for (let i = 0; i < 20; i++) t.advance(16);
    expect(vp.live.current.scrollX).toBeLessThan(sl - 50);
    t.pointer('up', L, -40, y, { alt: true, onDocument: true });
  });

  it('rows: an Alt-drag over cells out below follows x on the bottom row in view; past the bottom the rows scroll', async () => {
    const F = 'family.main';
    const t = await driven(case3(), { size: { [F]: [1000, 150] } });
    const l = t.layout(F);
    const vp = t.vp(F);
    vp.zoomTo(0, 5, 25);
    const x = (i: number) => t.xOf(F, 0, i);
    const rowY = (d: number) => l.scrolledTop + l.entryTop[l.rows.entry] + d * l.rows.height + l.rows.height / 2;
    const pd = t.ctx.selection.pending.current;
    t.pointer('down', F, x(8), rowY(0), { alt: true });
    t.pointer('move', F, x(9), rowY(0), { alt: true });
    t.pointer('move', F, x(12), l.height + 10, { alt: true, onDocument: true });
    const bottom = pd.rowEnd;
    expect([pd.endShown, bottom > 0]).toEqual([vp.axis?.shownOf(0, 12), true]);
    t.pointer('move', F, x(15), l.height + 10, { alt: true, onDocument: true });
    expect(pd.endShown).toBe(vp.axis?.shownOf(0, 15));
    if (l.vscroll) {
      const y0 = vp.live.current.scrollY;
      for (let i = 0; i < 10; i++) t.advance(16);
      expect(vp.live.current.scrollY).toBeGreaterThan(y0);
      expect(pd.rowEnd).toBeGreaterThan(bottom);
    }
    t.pointer('up', F, x(15), l.height + 10, { alt: true, onDocument: true });
  });

  it('the scroll past an edge reads how far past: faster further out', async () => {
    const run = async (past: number) => {
      const t = await text(160);
      const a = t.at(0, 5);
      press(t, a.x, a.y, [[a.x, t.l.height + past]]);
      const y0 = t.vp(S).live.current.scrollY;
      for (let i = 0; i < 10; i++) t.advance(16);
      const y = t.vp(S).live.current.scrollY - y0;
      t.pointer('up', S, a.x, t.l.height + past, { onDocument: true });
      return y;
    };
    const [near, far] = [await run(10), await run(60)];
    expect(near).toBeGreaterThan(0);
    expect(far).toBeGreaterThan(near * 2);
  });
});

describe('fast pointers and lost releases', () => {
  it("jumps, a release right after a jump, across lines, outside the section: the release's range committed, every time", async () => {
    const t = await text(300);
    const L = t.l.lines;
    const pd = t.ctx.selection.pending.current;
    const failures: string[] = [];
    for (let k = 0; k < 24; k++) {
      const a = (k * 37) % L.stop[2];
      const z = (k * 53 + 11) % L.stop[3];
      const how = k % 3 === 0 ? 'jump' : k % 3 === 1 ? 'release' : 'outside';
      t.vp(S).scrollTo(0, 0, 'api');
      const p0 = t.at(0, a);
      const p1 = t.at(0, z);
      t.advance(1000);
      t.pointer('down', S, p0.x, p0.y);
      // The end the release must commit, by the nearest-position rule.
      let want = z;
      if (how === 'jump') {
        // One move past the slop straight to the end, then the release there.
        t.pointer('move', S, p1.x, p1.y);
        if (pd.endShown !== z) failures.push(`${k} jump: the first drag event's end ${pd.endShown}, want ${z}`);
        t.pointer('up', S, p1.x, p1.y);
      } else if (how === 'release') {
        // A first move just past the slop, then the release far away with no move between.
        t.pointer('move', S, p0.x + 6, p0.y);
        t.pointer('up', S, p1.x, p1.y);
      } else {
        // A move to the end's letter, then the release outside, right of its line.
        t.pointer('move', S, p1.x, p1.y);
        t.pointer('up', S, t.l.width + 60, p1.y, { onDocument: true });
        // Past the line's end: forward its last position, backward the one after it (the next line's first).
        want = L.stop[p1.line] - 1 < a ? L.stop[p1.line] : L.stop[p1.line] - 1;
      }
      const got = JSON.stringify(positions(sel(t)));
      const exp = JSON.stringify([[0, range(Math.min(a, want), Math.max(a, want) + 1)]]);
      if (got !== exp) failures.push(`${k} ${how} ${a}→${want}: got ${got.slice(0, 60)}`);
      expect(pd.active).toBe(false);
    }
    expect(failures).toEqual([]);
  });

  it('a release the page lost (a move with no button, then another press): the range is committed at the last point', async () => {
    const t = await text();
    const a = t.at(0, 5);
    const z = t.at(0, t.l.lines.start[1] + 3);
    t.pointer('down', S, a.x, a.y);
    t.pointer('move', S, z.x, z.y);
    unpressedMove(t, 20, 20);
    expect([t.ctx.ui.drag.value, t.ctx.selection.pending.current.active]).toEqual([null, false]);
    expect(positions(sel(t))).toEqual([[0, range(5, t.l.lines.start[1] + 4)]]);
    t.advance(500);
    // The next press while the last one was never released: that one ends first, where it last was.
    const b = t.at(0, 30);
    t.pointer('down', S, a.x, a.y);
    t.pointer('move', S, b.x, b.y);
    t.pointer('down', S, z.x, z.y);
    expect(positions(sel(t))).toEqual([[0, range(5, 31)]]);
    t.pointer('up', S, z.x, z.y);
  });

  it('linear: a select drag follows its first drag event and its release', async () => {
    const F = 'family.main';
    const t = await driven(case3(), { size: { [F]: [1000, 150] } });
    t.vp(F).zoomTo(0, 5, 25);
    const l = t.layout(F);
    const y = l.entryTop[1] + 9;
    const x = (i: number) => t.xOf(F, 0, i);
    const cols = () => {
      const c = t.ctx.selection.current.value.columns;
      const out: number[] = [];
      if (c) OS.forEach(c.columns, (i: number) => out.push(i), null);
      return out;
    };
    t.pointer('down', F, x(8), y);
    t.pointer('move', F, x(12), y);
    expect(t.ctx.selection.pending.current.endShown).toBe(t.vp(F).axis?.shownOf(0, 12));
    t.pointer('up', F, x(16), y);
    expect(cols()).toEqual(range(8, 17));
  });
});

describe('the drag drawn', () => {
  it('wrapped, the pending range is drawn on each of its lines as the selection it becomes', async () => {
    const doc = createSeqViewSpec();
    const a = doc.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(30) });
    doc.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' }).block({ sequence: a }).track().letters();
    const t = await lit(doc.getDocument(), { width: 400, height: 600 });
    const l = t.layout(t.section);
    expect(l.lines.count).toBeGreaterThan(4);
    const [k0, k1] = [7, l.lines.start[3] + 4];
    const rects = () =>
      passRects(t.overlay(), Pass.overlay).filter(
        (r) => r[4] === SkinToken.selection || r[4] === SkinToken.selectionBorder,
      );
    t.ctx.selection.set(seqLoci(t.set, 'a', k0, k1 + 1), 'api');
    const committed = rects();
    expect(new Set(committed.map((r) => r[1])).size).toBe(4);
    t.ctx.selection.set(null, 'api');
    const pd = t.ctx.selection.pending.current;
    Object.assign(pd, { active: true, section: 0, block: 0, kind: 1, startShown: k1, endShown: k0 });
    (t.ctx.selection.pending as unknown as { bump(): void }).bump();
    expect(rects()).toEqual(committed);
    // A drag that removes: its outline alone, no wash.
    pd.remove = true;
    (t.ctx.selection.pending as unknown as { bump(): void }).bump();
    const removing = rects();
    expect([removing.some((r) => r[4] === SkinToken.selection), removing.length > 0]).toEqual([false, true]);
    pd.remove = false;
    pd.active = false;
    (t.ctx.selection.pending as unknown as { bump(): void }).bump();
    expect(rects()).toEqual([]);
  });
});
