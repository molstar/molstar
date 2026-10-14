/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard map and the cursor (IN4; seqspec:viewer#behaviors, seqspec:viewer#state,
// seqspec:interaction#accessibility): one tab stop per section, exact modifiers, the macOS key codes, a cursor kept in
// the block coordinate system and clamped on an axis or layout change, brought into view, reaching headers and the rows
// of a lazy run; through filtered, reversed and cyclic pieces; and `cursorRect`. The old stack's behavior for each
// finding was compared until SW1 (its parity layer, deleted).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { OrderedSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case3, case13 } from '../../core/_test/format/cases.js';
import { EntryKind, HoverKind, KeyBinding, type KeyCursor, Modifier, type Selection } from '../src/index.js';
import { cursorRect } from '../src/keyboard.js';
import type { Layouts } from '../src/layout.js';
import { rangedDoc, setPiece } from './axis-fixtures.js';
import { type Driven, driven } from './input-fixtures.js';
import { lazyLoaders } from './layout-fixtures.js';

const S = 'family.main';
const cursor = (t: Driven) => t.ctx.ui.keyCursor.value;
const family = async () => driven(case3(), { size: { [S]: [1000, 150] } });
const columns = (s: Selection) => {
  const out: number[] = [];
  if (s.columns) OrderedSet.forEach(s.columns.columns, (i: number) => out.push(i), null);
  return out;
};
/** `cursorRect` of the current cursor: x, y, w, h, or null. */
function rectOf(t: Driven, c: KeyCursor | null = cursor(t)) {
  const out = new Float64Array(6);
  return c && cursorRect(t.ctx, c, out) ? Array.from(out.subarray(0, 4)) : null;
}
/** The keys of a track's repeat: `n` presses, the cursor's block indices after each. */
const walk = (t: Driven, section: string, key: string, n: number) =>
  Array.from({ length: n }, () => {
    t.key(section, key);
    const c = cursor(t);
    return c?.header ? 'header' : (c?.index ?? null);
  });

describe('exact bindings', () => {
  it('KeyBinding: Ctrl, Meta and Alt exact; Shift only for a named key; `code` for the physical key', () => {
    const k = (key: string, modifiers = 0, code = '') => ({ section: 0, key, code, modifiers, typing: false, time: 0 });
    const m = (b: object, key: string, modifiers = 0, code = '') => KeyBinding.matches(b, k(key, modifiers, code));
    expect(m({ key: '+' }, '+')).toBe(true);
    expect(m({ key: '+' }, '+', Modifier.shift)).toBe(true); // Shift is in the value of a printable key
    expect(m({ key: '+' }, '+', Modifier.ctrl)).toBe(false);
    expect(m({ key: '+' }, '+', Modifier.meta)).toBe(false);
    expect(m({ key: '-' }, '-', Modifier.alt)).toBe(false);
    expect(m({ key: 'ArrowLeft' }, 'ArrowLeft', Modifier.alt)).toBe(false);
    expect(m({ key: 'Enter' }, 'Enter', Modifier.shift)).toBe(false); // a named key: Shift is a different binding
    expect(m({ key: 'Enter', shift: true }, 'Enter', Modifier.shift)).toBe(true);
    expect(m({ key: 'Enter', shift: true }, 'Enter')).toBe(false);
    expect(m({ code: 'Equal', alt: true }, '≠', Modifier.alt, 'Equal')).toBe(true);
    expect(m({ code: 'Equal', alt: true }, '=', Modifier.alt, 'Equal')).toBe(true);
    expect(m({ code: 'Equal', alt: true }, '=', 0, 'Equal')).toBe(false);
    expect(m({ code: 'Equal', alt: true }, '≠', Modifier.alt, 'Minus')).toBe(false);
    expect(m({}, 'x')).toBe(false); // a trigger names a key or a code
  });

  it("Ctrl and Cmd with + and -, and Alt with an arrow, are the page's: not handled, nothing moves", async () => {
    const t = await family();
    t.focus(S);
    const [cell, at, scale] = [t.vp(S).live.current.cell, cursor(t), t.vp(S).live.current.rowScale];
    expect(t.key(S, '+', { ctrl: true })).toBe(false);
    expect(t.key(S, '-', { ctrl: true })).toBe(false);
    expect(t.key(S, '-', { meta: true })).toBe(false);
    expect(t.key(S, '=', { meta: true })).toBe(false);
    expect(t.key(S, 'ArrowLeft', { alt: true })).toBe(false);
    expect(t.key(S, 'ArrowRight', { alt: true })).toBe(false);
    expect(t.key(S, 'Enter', { ctrl: true })).toBe(false);
    expect([t.vp(S).live.current.cell, cursor(t), t.vp(S).live.current.rowScale]).toEqual([cell, at, scale]);
  });

  it('Shift+`+` zooms in and Shift+`_` out; plain `=` and `-` too; Shift+Enter is not Enter', async () => {
    const t = await family();
    t.focus(S);
    const s = t.vp(S).live.current;
    const c0 = s.cell;
    expect(t.key(S, '+', { shift: true })).toBe(true);
    const c1 = s.cell;
    expect(c1).toBeGreaterThan(c0);
    expect(t.key(S, '_', { shift: true })).toBe(true);
    expect(s.cell).toBeLessThan(c1);
    t.key(S, '=');
    expect(s.cell).toBeGreaterThan(c0 * 0.99);
    t.key(S, '-');
    expect(s.cell).toBeLessThan(c1);
  });

  it("Enter sets the selection and Shift+Enter adds to it, both as the select behavior's clicks", async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowUp'); // from the first row up onto the pinned columns
    expect(cursor(t)?.track).toBe('family.flvr-arginine');
    expect(t.ctx.hover.live.current.kind).toBe(HoverKind.column);
    t.key(S, 'Enter');
    expect(columns(t.ctx.selection.current.value)).toEqual([0]);
    t.key(S, 'ArrowRight');
    t.key(S, 'Enter', { shift: true });
    expect(columns(t.ctx.selection.current.value)).toEqual([0, 1]);
    t.key(S, 'ArrowRight');
    t.key(S, 'Enter');
    expect(columns(t.ctx.selection.current.value)).toEqual([2]);
    t.key(S, 'Escape');
    expect(t.ctx.selection.current.value.columns).toBeNull();
    // Space is Enter on data.
    t.key(S, ' ');
    expect(columns(t.ctx.selection.current.value)).toEqual([2]);
  });

  it('Alt+= and Alt+- are row zoom by `code`, whatever character macOS gives them', async () => {
    const t = await family();
    t.focus(S);
    const s = t.vp(S).live.current;
    const r0 = s.rowScale;
    expect(t.key(S, '≠', { alt: true, code: 'Equal' })).toBe(true);
    const r1 = s.rowScale;
    expect(r1).toBeGreaterThan(r0);
    expect(t.key(S, '–', { alt: true, code: 'Minus' })).toBe(true);
    expect(s.rowScale).toBeLessThan(r1);
    expect(t.key(S, '=', { alt: true, code: 'Equal' })).toBe(true); // Windows and Linux keep the character
    expect(s.rowScale).toBeGreaterThan(r0 * 0.99);
    expect(t.key(S, '≠', { alt: true, ctrl: true, code: 'Equal' })).toBe(false);
  });

  it("Tab and Shift+Tab are never handled; a key typed in an editable element is the page's", async () => {
    const t = await family();
    t.focus(S);
    const at = cursor(t);
    expect([t.key(S, 'Tab'), t.key(S, 'Tab', { shift: true }), t.key(S, 'ArrowRight', { typing: true })]).toEqual([
      false,
      false,
      false,
    ]);
    expect(cursor(t)).toBe(at);
    expect(t.key(S, 'q')).toBe(false);
  });

  it('Shift+F10 and ContextMenu are handled and leave the cursor where it is; F10 alone is not', async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowRight');
    const at = cursor(t);
    expect(t.key(S, 'F10', { shift: true })).toBe(true);
    expect(t.key(S, 'ContextMenu')).toBe(true);
    expect(t.key(S, 'F10')).toBe(false);
    expect(cursor(t)).toEqual(at);
    // The UI reads the same binding to open the cursor's track menu.
    const bindings = t.ctx.behaviors.get<{ bindings: Record<string, never> }>('keyboard')?.forSection(0)?.bindings;
    expect(
      KeyBinding.match(bindings?.menu ?? [], {
        section: 0,
        key: 'F10',
        code: 'F10',
        modifiers: Modifier.shift,
        typing: false,
        time: 0,
      }),
    ).toBe(true);
  });
});

describe('the cursor on the headers and rows', () => {
  it('focus shows it on the first position in view of the first scrolled entry; blur clears it and the hover; focus again brings it back', async () => {
    const t = await family();
    expect(cursor(t)).toBeNull();
    t.focus(S);
    expect(cursor(t)).toEqual({
      section: S,
      track: 'family.members',
      header: false,
      block: 0,
      index: 0,
      row: 0,
      displayRow: 0,
    });
    expect(t.ctx.hover.live.current.owner).toBe('keyboard');
    t.key(S, 'ArrowRight');
    t.key(S, 'ArrowRight');
    t.key(S, 'ArrowDown');
    const at = cursor(t);
    t.focus(S, false);
    expect([cursor(t), t.ctx.hover.live.current.owner, t.ctx.labels.current.value]).toEqual([null, '', null]);
    t.focus(S);
    expect(cursor(t)).toEqual(at);
    expect(t.ctx.hover.live.current.owner).toBe('keyboard');
  });

  it("Up and Down step the alignment's rows (display order) and the pinned tracks, and stop at the ends", async () => {
    const t = await family();
    t.focus(S);
    const row = () => [cursor(t)?.track, cursor(t)?.displayRow, cursor(t)?.row];
    t.key(S, 'ArrowDown');
    t.key(S, 'ArrowDown');
    t.key(S, 'ArrowDown');
    // Tree order ((SRC, GRB2), (PTN11, P85A)): data rows 0, 2, 1, 3.
    expect(row()).toEqual(['family.members', 3, 3]);
    expect(t.key(S, 'ArrowDown')).toBe(true);
    expect(row()).toEqual(['family.members', 3, 3]);
    t.key(S, 'ArrowUp');
    expect(row()).toEqual(['family.members', 2, 1]);
    for (let i = 0; i < 3; i++) t.key(S, 'ArrowUp');
    expect(cursor(t)?.track).toBe('family.flvr-arginine');
    const above: (string | null | undefined)[] = [];
    for (let i = 0; i < 6; i++) {
      t.key(S, 'ArrowUp');
      above.push(cursor(t)?.track);
    }
    expect(above).toEqual([
      'family.clade-b',
      'family.clade-a',
      'family.conservation',
      'family.logo',
      'family.logo',
      'family.logo',
    ]);
  });

  it('Left off the first position, and Home, go to the header; Right comes back; a header keeps its place under Up and Down', async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowLeft');
    expect(cursor(t)).toMatchObject({ header: true, block: -1, index: -1, track: 'family.members', displayRow: 0 });
    expect(t.ctx.hover.live.current.kind).toBe(HoverKind.header);
    expect(rectOf(t)).toEqual([0, t.layout(S).scrolledTop, t.layout(S).headerWidth, 18]);
    t.key(S, 'ArrowDown');
    expect(cursor(t)).toMatchObject({ header: true, track: 'family.members', displayRow: 1 });
    t.key(S, 'ArrowRight');
    expect(cursor(t)).toMatchObject({ header: false, index: 0, displayRow: 1 });
    t.key(S, 'ArrowRight');
    t.key(S, 'Home');
    expect(cursor(t)?.header).toBe(true);
    t.key(S, 'End');
    const last = (t.vp(S).axis?.shown[0] ?? 0) - 1;
    expect(cursor(t)).toMatchObject({ header: false, index: last, displayRow: 1 });
    t.key(S, 'ArrowRight'); // the last position: stays
    expect(cursor(t)?.index).toBe(last);
  });

  it('Enter on a header activates it; Shift+Enter too; a stale header never beats Enter on data', async () => {
    const t = await family();
    const activated: string[] = [];
    t.ctx.events.activate.subscribe((e) => activated.push(`${e.section} ${e.track} ${e.origin}`));
    t.focus(S);
    t.key(S, 'ArrowRight');
    t.key(S, 'Enter');
    expect(activated).toEqual([]);
    expect(t.ctx.selection.current.value.aligned).not.toBeNull();
    t.key(S, 'Home');
    t.key(S, 'Enter');
    t.key(S, 'Enter', { shift: true });
    expect(activated).toEqual([`${S} family.members keyboard`, `${S} family.members keyboard`]);
    t.focus(S, false);
    t.focus(S);
    t.key(S, 'ArrowRight');
    t.key(S, 'Enter');
    expect(activated.length).toBe(2);
  });

  it("a collapsed clade's wedge holds no element: focus and Down skip it on data, Down stops on it on a header, Space expands it", async () => {
    const t = await family();
    t.ctx.ui.toggleClade('sh2', 1);
    const rows = t.layout(S).rows;
    // Display row 0 is the wedge of (SRC, GRB2); rows 1 and 2 are PTN11 and P85A.
    expect([rows.count, Array.from(rows.node ?? [])]).toEqual([3, [1, -1, -1]]);
    t.focus(S);
    expect(cursor(t)).toMatchObject({ track: 'family.members', displayRow: 1, row: 1 });
    t.key(S, 'ArrowUp');
    expect(cursor(t)?.track).toBe('family.flvr-arginine'); // not onto the wedge
    t.key(S, 'Home');
    t.key(S, 'ArrowDown');
    expect(cursor(t)).toMatchObject({ header: true, track: 'family.members', displayRow: 0, row: -1 });
    t.key(S, ' ');
    expect(t.layout(S).rows.count).toBe(4);
    expect(Array.from(t.ctx.ui.clades.value.get('sh2') ?? [])).toEqual([]);
  });
});

describe('keep in view', () => {
  /** One long sequence at a readable zoom. */
  async function long() {
    const b = createSeqViewSpec();
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: b.data.sequence({ ref: 's', value: 'ACDEFGHIKL'.repeat(600) }) })
      .track({ ref: 'v.l' })
      .letters();
    const t = await driven(b.getDocument(), { size: { 'v.main': [600, 120] } });
    t.vp('v.main').zoomTo(0, 2990, 3030);
    return t;
  }

  it('a mouse pan to position 5,000, then an arrow: the cursor is back in view, centered; one more step moves a cell', async () => {
    const t = await long();
    const vp = t.vp('v.main');
    const l = t.layout('v.main');
    t.focus('v.main');
    t.key('v.main', 'ArrowRight');
    expect(cursor(t)?.index).toBe(2991);
    expect(rectOf(t)).not.toBeNull();
    vp.scrollTo(5000 * vp.live.current.cell, 0, 'pointer');
    expect(rectOf(t)).toBeNull();
    t.key('v.main', 'ArrowRight');
    const r = rectOf(t);
    expect(cursor(t)?.index).toBe(2992);
    expect(r && r[0] >= l.dataLeft && r[0] + r[2] <= l.dataRight).toBe(true);
    expect(Math.abs((r?.[0] ?? 0) + (r?.[2] ?? 0) / 2 - (l.dataLeft + l.dataRight) / 2)).toBeLessThan(2);
    const x = r?.[0] ?? 0;
    t.key('v.main', 'ArrowRight');
    expect(Math.abs((rectOf(t)?.[0] ?? 0) - x - vp.live.current.cell)).toBeLessThan(1); // snapped to device px
    // Past the right edge it scrolls by what is needed, not by four cells.
    const s0 = vp.live.current.scrollX;
    for (let i = 0; i < 200; i++) t.key('v.main', 'ArrowRight');
    expect(rectOf(t)).not.toBeNull();
    expect(vp.live.current.scrollX).toBeGreaterThan(s0);
    // End goes to the last position, far away, and Home to the header.
    t.key('v.main', 'End');
    expect([cursor(t)?.index, rectOf(t) !== null]).toEqual([5999, true]);
  });

  it('moving down onto a row scrolled out of view scrolls it in, and the hover is its cell, never a stale one', async () => {
    const t = await family();
    const vp = t.vp(S);
    t.focus(S);
    t.key(S, 'ArrowDown');
    expect(vp.live.current.scrollY).toBe(0);
    t.key(S, 'ArrowDown');
    t.key(S, 'ArrowDown');
    const l = t.layout(S);
    expect(vp.live.current.scrollY).toBeGreaterThan(0);
    const r = rectOf(t);
    expect(r && r[1] >= l.scrolledTop && r[1] + r[3] <= l.height + 1e-9).toBe(true);
    const h = t.ctx.hover.live.current;
    expect([h.owner, h.kind, h.lane, h.region]).toEqual(['keyboard', HoverKind.aligned, 3, 1]);
    // Scrolled away with the mouse, then Up: back in view.
    vp.scrollTo(vp.live.current.scrollX, 0, 'pointer');
    expect(rectOf(t)).toBeNull();
    t.key(S, 'ArrowUp');
    expect(rectOf(t)).not.toBeNull();
  });

  it('a wrapped section: End scrolls the last line in, and Left comes back up', async () => {
    const b = createSeqViewSpec();
    b.view({ ref: 'v', label: 'V', layout: 'text' })
      .block({ sequence: b.data.sequence({ ref: 's', value: 'ACDEFGHIKL'.repeat(40) }) })
      .track({ ref: 'v.l' })
      .letters();
    const t = await driven(b.getDocument(), { size: { 'v.main': [400, 120] } });
    t.focus('v.main');
    expect(rectOf(t)).not.toBeNull();
    t.key('v.main', 'End');
    expect(cursor(t)?.index).toBe(399);
    expect(t.vp('v.main').live.current.scrollY).toBeGreaterThan(0);
    expect(rectOf(t)).not.toBeNull();
    for (let i = 0; i < 300; i++) t.key('v.main', 'ArrowLeft');
    expect(cursor(t)?.index).toBe(99);
    expect(rectOf(t)).not.toBeNull();
  });

  it("a wrapped section with group-label bands: the rect is the track's band of the cursor's line; Home and Right stay on it", async () => {
    const b = createSeqViewSpec();
    b.view({ ref: 'v', label: 'V', layout: 'text' })
      .block({ sequence: b.data.sequence({ ref: 's', value: 'ACDEFGHIKL'.repeat(40) }) })
      .track({ ref: 'v.l' })
      .letters();
    const t = await driven(b.getDocument(), { size: { 'v.main': [400, 600] } });
    const l = t.layout('v.main');
    const top = (line: number) => {
      for (let e = 0; e < l.entryCount; e++)
        if (l.entryKind[e] === EntryKind.track && l.entryRow[e] === line) return l.scrolledTop + l.entryTop[e];
      return Number.NaN;
    };
    expect(Array.from(l.entryKind.subarray(0, l.entryCount)).filter((k) => k === EntryKind.labels).length).toBe(
      l.lines.count,
    );
    t.focus('v.main');
    // Across the first line break: the rect's band is the line's track entry, below its labels band, and the hit at its
    // center is the cursor's position.
    for (let i = 0; i < l.lines.start[1] + 2; i++) {
      const c = cursor(t) as KeyCursor;
      expect([rectOf(t)?.[1], t.ctx.hover.live.current.index]).toEqual([top(c.displayRow), c.index]);
      t.key('v.main', 'ArrowRight');
    }
    const line = cursor(t)?.displayRow ?? -1;
    expect(line).toBe(1);
    t.key('v.main', 'Home');
    expect(cursor(t)).toMatchObject({ header: true, displayRow: 1 });
    expect(rectOf(t)?.[1]).toBe(top(1));
    t.key('v.main', 'ArrowLeft');
    expect(cursor(t)).toMatchObject({ header: true, displayRow: 1 });
    t.key('v.main', 'ArrowRight');
    expect(cursor(t)).toMatchObject({ header: false, displayRow: 1, index: l.lines.start[1] });
  });
});

describe('a stacked section', () => {
  it("Down passes from a block's last track to the next block's first, keeping the position clamped; Right stays within its block", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'a', value: 'MKVLATGHEKMKVLATGHEK' });
    b.data.sequence({ ref: 'b', value: 'MKVLA' });
    const v = b.view({ ref: 'v', label: 'V', arrangement: 'stacked' });
    v.track({ ref: 'v.ruler' }).ruler();
    v.track({ ref: 'v.letters' }).letters();
    const t = await driven(b.getDocument(), { size: { 'v.main': [800, 400] } });
    t.focus('v.main');
    expect(cursor(t)).toMatchObject({ track: 'v.ruler', block: 0, index: 0 });
    t.key('v.main', 'End');
    expect(cursor(t)).toMatchObject({ block: 0, index: 19 });
    t.key('v.main', 'ArrowDown');
    expect(cursor(t)).toMatchObject({ track: 'v.letters', block: 0, index: 19 });
    t.key('v.main', 'ArrowDown');
    expect(cursor(t)).toMatchObject({ track: 'v.ruler', block: 1, index: 4 });
    expect(t.ctx.hover.live.current).toMatchObject({ owner: 'keyboard', block: 1, index: 4 });
    t.key('v.main', 'ArrowRight');
    expect(cursor(t)).toMatchObject({ block: 1, index: 4 });
    t.key('v.main', 'ArrowLeft');
    t.key('v.main', 'ArrowUp'); // past the next block's title
    expect(cursor(t)).toMatchObject({ track: 'v.letters', block: 0, index: 3 });
    t.key('v.main', 'ArrowUp');
    expect(cursor(t)).toMatchObject({ track: 'v.ruler', block: 0, index: 3 });
    expect(rectOf(t)?.every(Number.isFinite)).toBe(true);
  });
});

describe("a lazy run's rows", () => {
  it('Down and Up step the rows of a million-row alignment, scrolling them in; the row is its data row', async () => {
    const [uri, n] = ['https://example.org/huge', 1_000_000];
    const b = createSeqViewSpec();
    const aln = b.data.alignment({
      ref: 'huge',
      row_count: n,
      column_count: 200,
      source: { uri, protocol: 'seqstar-tiles/1' },
    });
    b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ ref: 'v.rows' }).rows();
    const t = await driven(b.getDocument(), { loaders: lazyLoaders(uri, n, 200) });
    t.focus('v.main');
    expect(cursor(t)).toMatchObject({ track: 'v.rows', row: 0, displayRow: 0, index: 0 });
    for (let i = 0; i < 40; i++) t.key('v.main', 'ArrowDown');
    expect(cursor(t)).toMatchObject({ row: 40, displayRow: 40 });
    const r = rectOf(t);
    const l = t.layout('v.main');
    expect(r && r[1] >= l.scrolledTop && r[1] + r[3] <= l.height).toBe(true);
    expect(t.vp('v.main').live.current.scrollY).toBeGreaterThan(0);
    expect(t.ctx.hover.live.current).toMatchObject({ owner: 'keyboard', kind: HoverKind.aligned });
    t.key('v.main', 'ArrowUp');
    expect(cursor(t)?.displayRow).toBe(39);
    // Row zoom keeps the cursor in view.
    t.key('v.main', '≠', { alt: true, code: 'Equal' });
    expect(rectOf(t)).not.toBeNull();
  });
});

describe('through filtered, reversed and cyclic pieces', () => {
  /** The block indices Right visits from the first position, then Left back, over a piece. */
  async function visits(view: string, piece?: [number, number, number, number]) {
    const t = await driven(rangedDoc(), { size: { [`${view}.main`]: [600, 120] } });
    t.ctx.views.setCurrent(view);
    await t.ctx.graph.update();
    const sec = `${view}.main`;
    const m = t.ctx.graph.assembled.value?.view;
    if (piece && m) {
      setPiece(m.resolved, 0, piece[0], piece[1], piece[2], piece[3]);
      (t.layouts as Layouts).update(m);
    }
    t.focus(sec);
    const n = t.vp(sec).axis?.shown[0] ?? 0;
    const right = [cursor(t)?.index, ...walk(t, sec, 'ArrowRight', n - 1)];
    const edge = walk(t, sec, 'ArrowRight', 1);
    const left = walk(t, sec, 'ArrowLeft', n);
    return { t, sec, n, right, edge, left };
  }

  it("a filter's break is a step: the filtered positions in order, across the gap, in block indices", async () => {
    const r = await visits('f');
    // site [7, 10) and [12, 14) over the unit [5, 15): five positions.
    expect([r.n, r.right, r.edge]).toEqual([5, [7, 8, 9, 12, 13], [13]]);
    expect(r.left).toEqual([12, 9, 8, 7, 'header']);
    const rect = rectOf(r.t);
    expect(rect).toEqual([0, rect?.[1], r.t.layout(r.sec).headerWidth, rect?.[3]]);
  });

  it('a reversed piece steps from its first shown position, which is its highest index', async () => {
    const r = await visits('w', [0, 20, -1, 20]);
    expect(r.right).toEqual(Array.from({ length: 20 }, (_, i) => 19 - i));
    expect(r.left.slice(0, 3)).toEqual([1, 2, 3]);
  });

  it('a piece across a cyclic origin steps over the origin', async () => {
    const r = await visits('w', [15, 15, 1, 20]);
    expect(r.right).toEqual([
      ...Array.from({ length: 5 }, (_, i) => 15 + i),
      ...Array.from({ length: 15 }, (_, i) => i),
    ]);
    // Pressing Right 19 times crossed the origin; the cursor's rect is finite and in view.
    expect(rectOf(r.t)).not.toBeNull();
  });
});

describe('clamped on an axis or layout change', () => {
  const patch = (t: Driven, ref: string, params: Record<string, unknown>) =>
    t.ctx.applyPatch({ format: 'seq_view_spec_patch', ops: [{ op: 'set', ref, params }] });

  it("a patch that shortens the unit's range moves the cursor to the nearest shown position, below; then above", async () => {
    const t = await driven(rangedDoc(), { size: { 'v.main': [600, 120] } });
    t.focus('v.main');
    for (let i = 0; i < 9; i++) t.key('v.main', 'ArrowRight');
    expect(cursor(t)).toMatchObject({ block: 0, index: 14 });
    await patch(t, 'u', { range: [5, 10] });
    expect(cursor(t)).toMatchObject({ block: 0, index: 9 });
    expect(t.ctx.hover.live.current).toMatchObject({ owner: 'keyboard', index: 9 });
    expect(rectOf(t)?.every(Number.isFinite)).toBe(true);
    await patch(t, 'u', { range: [12, 15] });
    expect(cursor(t)).toMatchObject({ block: 0, index: 12 });
    t.key('v.main', 'ArrowRight');
    expect(cursor(t)?.index).toBe(13);
    expect(rectOf(t)?.every(Number.isFinite)).toBe(true);
  });

  it("a relayout leaves the pointer's hover with the pointer; the cursor's comes back once nobody holds it (VS-01)", async () => {
    const t = await family();
    t.focus(S);
    expect(t.ctx.hover.live.current.owner).toBe('keyboard');
    t.pointer('move', S, t.xOf(S, 0, 8), t.layout(S).scrolledTop + 18 + 9);
    expect(t.ctx.hover.live.current).toMatchObject({ owner: 'pointer', index: 8 });
    t.ctx.viewports.get(S)?.setSize(1000, 200, 1);
    t.advance(50);
    expect(t.ctx.hover.live.current).toMatchObject({ owner: 'pointer', index: 8 });
    t.ctx.setRenderProps({ rowHeight: 20 } as never);
    t.advance(50);
    expect(t.ctx.hover.live.current.owner).toBe('pointer');
    expect(cursor(t)?.index).toBe(0);
    // The pointer leaves: nobody holds the hover, and the next relayout shows the cursor's again.
    t.ctx.input.leave(t.ctx.layout.sectionIndex(S));
    expect(t.ctx.hover.live.current.owner).toBe('');
    t.ctx.viewports.get(S)?.setSize(1000, 150, 1);
    t.advance(50);
    expect(t.ctx.hover.live.current).toMatchObject({ owner: 'keyboard', index: 0 });
  });

  it("a section with fewer blocks: the cursor's block is clamped, its position the nearest, nothing NaN", async () => {
    const b = createSeqViewSpec();
    const v = b.view({ ref: 'v', label: 'V' });
    for (const i of [0, 1, 2]) v.block({ sequence: b.data.sequence({ ref: `c${i}`, value: 'MKVLAAGIVG' }) });
    v.track({ ref: 'v.l' }).letters();
    const t = await driven(b.getDocument(), { size: { 'v.main': [800, 120] } });
    t.focus('v.main');
    t.key('v.main', 'End');
    expect(cursor(t)).toMatchObject({ block: 2, index: 9 });
    const m = t.ctx.graph.assembled.value?.view;
    if (!m) throw new Error('no view');
    (m.resolved.sections.blockCount as Int32Array)[0] = 1;
    (t.layouts as Layouts).update(m);
    t.layout('v.main');
    expect(cursor(t)).toMatchObject({ block: 0, index: 9 });
    expect(rectOf(t)?.every(Number.isFinite)).toBe(true);
    expect(t.key('v.main', 'ArrowLeft')).toBe(true);
    expect(cursor(t)).toMatchObject({ block: 0, index: 8 });
  });

  it("a collapse that removes the cursor's row, and a hidden track, move it to the nearest entry still there", async () => {
    const t = await family();
    t.focus(S);
    for (let i = 0; i < 3; i++) t.key(S, 'ArrowDown');
    expect(cursor(t)).toMatchObject({ displayRow: 3, row: 3 });
    // Collapsing the clade (SRC, GRB2) folds display rows 0 and 1 into a wedge: the cursor's row 3 becomes row 2.
    t.ctx.ui.toggleClade('sh2', 1);
    expect(t.layout(S).rows.count).toBe(3);
    const c = cursor(t);
    expect([c?.track, c?.displayRow, c?.row]).toEqual(['family.members', 2, 3]);
    expect(rectOf(t)).not.toBeNull();
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'family.members', params: { hidden: true } }],
    });
    const h = cursor(t);
    expect(h?.track).not.toBe('family.members');
    expect(h?.track).not.toBeNull();
    expect(rectOf(t)).not.toBeNull();
  });

  it('a collapse keeps the cursor on its data row, wherever its display row went', async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowDown');
    t.key(S, 'ArrowDown');
    expect(cursor(t)).toMatchObject({ displayRow: 2, row: 1 }); // PTN11
    t.ctx.ui.toggleClade('sh2', 1);
    t.layout(S);
    expect(cursor(t)).toMatchObject({ displayRow: 1, row: 1 });
    const h = t.ctx.hover.live.current;
    expect([h.owner, h.aligned.row]).toEqual(['keyboard', 1]);
  });

  it("Space on a leaf's header collapses its clade and stays on the wedge; Space again expands it; Right on a wedge stays", async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowDown'); // GRB2, data row 2, in the clade (SRC, GRB2): node 1
    t.key(S, 'Home');
    t.key(S, ' ');
    expect(Array.from(t.ctx.ui.clades.value.get('sh2') ?? [])).toEqual([1]);
    t.layout(S);
    expect(cursor(t)).toMatchObject({ header: true, track: 'family.members', displayRow: 0, row: -1 });
    t.key(S, 'ArrowRight');
    expect(cursor(t)?.header).toBe(true);
    t.key(S, ' ');
    t.layout(S);
    expect(Array.from(t.ctx.ui.clades.value.get('sh2') ?? [])).toEqual([]);
    expect(cursor(t)).toMatchObject({ header: true, displayRow: 0 });
  });

  it("collapsing a split track's rows from its header keeps the cursor on the track", async () => {
    const F = 'features.main';
    const NV = 'features.natural-variants';
    const t = await driven(case13(), { size: { [F]: [1000, 900] } });
    t.focus(F);
    for (let i = 0; i < 40 && !(cursor(t)?.track === NV && (cursor(t)?.displayRow ?? 0) > 2); i++)
      t.key(F, 'ArrowDown');
    expect(cursor(t)?.track).toBe(NV);
    t.key(F, 'Home');
    t.key(F, ' ');
    t.layout(F);
    expect(cursor(t)).toMatchObject({ header: true, track: NV });
    expect(rectOf(t)).not.toBeNull();
  });

  it('a cursor in a section the view switch removed goes, with its hover', async () => {
    const t = await driven(rangedDoc(), { size: { 'v.main': [600, 120] } });
    t.focus('v.main');
    t.key('v.main', 'ArrowRight');
    expect(t.ctx.hover.live.current.owner).toBe('keyboard');
    t.ctx.views.setCurrent('w');
    await t.ctx.graph.update();
    expect([cursor(t), t.ctx.hover.live.current.owner]).toEqual([null, '']);
  });

  it('a cursor naming a block, a row or a track that is not there is refused by cursorRect, never NaN', async () => {
    const t = await family();
    t.focus(S);
    const c = cursor(t) as KeyCursor;
    expect(rectOf(t, { ...c, block: 4 })).toBeNull();
    expect(rectOf(t, { ...c, index: 9999 })).toBeNull();
    expect(rectOf(t, { ...c, displayRow: 99, row: 99 })).toBeNull();
    // The data row is the cursor's: a stale display row is drawn where that row is now.
    expect(rectOf(t, { ...c, displayRow: 99 })).toEqual(rectOf(t, c));
    expect(rectOf(t, { ...c, track: 'nope' })).toBeNull();
    expect(rectOf(t, { ...c, section: 'nope' })).toBeNull();
    expect(rectOf(t, c)?.every(Number.isFinite)).toBe(true);
  });
});

describe("cursorRect (the ring and the tooltip's anchor)", () => {
  it("at an entry, a row and a header: the cell of the position as drawn in the entry's band", async () => {
    const t = await family();
    const l = t.layout(S);
    const v = t.ctx.layout.axis(S);
    if (!v) throw new Error('no view');
    t.focus(S);
    t.key(S, 'ArrowRight');
    t.key(S, 'ArrowRight');
    // The cell as drawn: less its one device px gap on the right (DPR 1).
    const cell = t.vp(S).live.current.cell - 1;
    // A row of the alignment: scrolled, 18 px high at the layout's top plus the row's.
    expect(rectOf(t)).toEqual([v.x(0, 2 as never), l.scrolledTop, cell, 18]);
    // A pinned entry: its own band, less the lane's padding above and below.
    t.key(S, 'ArrowUp');
    const pad = l.entryHeight[4] > 10 ? 2 : 0;
    expect(rectOf(t)).toEqual([v.x(0, 2 as never), l.entryTop[4] + pad, cell, l.entryHeight[4] - 2 * pad]);
    // A header: the header column of the entry.
    t.key(S, 'Home');
    expect(rectOf(t)).toEqual([0, l.entryTop[4], l.headerWidth, l.entryHeight[4]]);
  });
});

describe("the behavior's lifecycle", () => {
  it("removing it mid-cursor clears the cursor and the keyboard's hover; the keys are no longer handled", async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowRight');
    expect([cursor(t) !== null, t.ctx.hover.live.current.owner]).toEqual([true, 'keyboard']);
    t.ctx.behaviors.get('keyboard')?.remove();
    expect([cursor(t), t.ctx.hover.live.current.owner]).toEqual([null, '']);
    expect(t.key(S, 'ArrowRight')).toBe(false);
  });

  it("params: a binding changed per section; the pointer's hover is not the cursor's (a leave keeps the cursor's)", async () => {
    const t = await family();
    const handle = t.ctx.behaviors.get<{ bindings: Record<string, readonly object[]> }>('keyboard');
    handle?.setParams({ bindings: { ...handle.params.value.bindings, right: [{ key: 'l' }] } });
    t.focus(S);
    expect(t.key(S, 'ArrowRight')).toBe(false);
    expect(t.key(S, 'l')).toBe(true);
    expect(cursor(t)?.index).toBe(1);
    t.ctx.input.leave(t.ctx.layout.sectionIndex(S));
    expect(t.ctx.hover.live.current.owner).toBe('keyboard');
  });

  it("the keyboard's hover has a label: the live region announces the cursor's position, not the pointer's", async () => {
    const t = await family();
    t.focus(S);
    t.key(S, 'ArrowRight');
    expect(t.ctx.labels.current.value?.title).toBeTruthy();
    t.pointer('move', S, 500, t.layout(S).scrolledTop + 27);
    expect(t.ctx.hover.live.current.owner).toBe('pointer');
    // The cursor remembers its place while the pointer hovers.
    expect(cursor(t)?.index).toBe(1);
  });
});
