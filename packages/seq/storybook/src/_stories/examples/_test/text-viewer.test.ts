/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text viewer, headless on a fake clock, composed as the text stories compose it: the interaction matrix's cells
// for a text section as binding tests (the letter size about the pointer's line, `0` back to the preset's, Space+drag
// and Shift+Alt+arrows scrolling lines, a drag selecting in reading order, ↑↓ a line, Shift+↑↓ extending by a line, a
// double-click revealing the element's first line with the size kept, a drag past the box's edge scrolling it by the
// clock), a linked selection revealed, and the letters copied as shown.

import { PlacePoint, type TextPlacement } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { P53 } from '@molstar/seqx-examples/p53/data';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const System = '/data/p53';

/** The text view at 360 px: its lines overflow their box, so it scrolls. */
async function open(id = 'text') {
  viewer = await openExample(P53Views.find((v) => v.id === id)!, { gene: '' });
  viewer.input.mac = false;
  viewer.resize(360);
  viewer.draw();
  return viewer;
}

const textOf = (v: SeqViewer) => v.frame!.layout.sections.find((s) => s.section.layout === 'text')!;
const placement = (v: SeqViewer) => textOf(v).placement as TextPlacement;
const address = (v: SeqViewer) => textOf(v).section.address;
const upl = (v: SeqViewer) => placement(v).geometry!.unitsPerLine;

/** Where a display unit is drawn, on the page (the middle of its cell in its line's first band). */
function at(v: SeqViewer, u: number): { x: number; y: number } {
  const pt = new PlacePoint();
  expect(placement(v).place(u, pt)).toBe(true);
  const a = placement(v).strips.axes[pt.strip];
  return { x: a.xOfUnit(u + 0.5), y: pt.y };
}

async function flush(v: SeqViewer, ms = 50) {
  await vi.advanceTimersByTimeAsync(ms);
  v.draw();
}

const selected = (v: SeqViewer) => [...(v.selection.value.value.shown.parts[0]?.ranges ?? [])];

describe('a text section: the interaction matrix', () => {
  it('Cmd/Ctrl+wheel steps the letter size keeping the pointer’s line; 0 the preset’s; - and = step it', async () => {
    const v = await open();
    const g = placement(v).geometry!;
    const u = 2 * g.unitsPerLine + 3;
    const p = at(v, u);
    expect(v.input.wheel(p.x, p.y, { deltaX: 0, deltaY: -100, ctrlKey: true })).toBe(true);
    await flush(v);
    const q = placement(v);
    expect(v.view.get(address(v))!.scale).toBeGreaterThan(g.cell);
    expect(q.geometry!.cell).toBeGreaterThan(g.cell);
    // the line holding the unit stays two lines below the top
    expect(Math.floor(u / q.geometry!.unitsPerLine) - Math.floor(q.anchor / q.geometry!.unitsPerLine)).toBe(2);
    expect(v.input.key({ key: '0' })).toBe(true);
    await flush(v);
    expect(v.view.get(address(v))!.scale).toBeUndefined();
    expect(placement(v).geometry!.cell).toBe(g.cell);
    v.input.key({ key: '=' });
    await flush(v, 400);
    expect(placement(v).geometry!.cell).toBeGreaterThan(g.cell);
    // never below the letters' threshold
    for (let i = 0; i < 12; i++) v.input.key({ key: '-' });
    await flush(v, 400);
    expect(placement(v).geometry!.cell).toBeGreaterThanOrEqual(10);
  });

  it('Space+drag scrolls lines; Shift+Alt+↑↓ scroll a line', async () => {
    const v = await open();
    const p = at(v, 3);
    v.input.setSpace(true);
    expect(v.input.down(p.x, p.y)).toBe(true);
    v.input.move(p.x, p.y - 90);
    v.input.up();
    v.input.setSpace(false);
    await flush(v);
    const scrolled = placement(v).anchor;
    expect(scrolled).toBeGreaterThan(upl(v));
    expect(v.input.key({ key: 'ArrowUp', shiftKey: true, altKey: true })).toBe(true);
    await flush(v);
    expect(placement(v).anchor).toBeCloseTo(scrolled - upl(v), 6);
    v.input.key({ key: 'ArrowDown', shiftKey: true, altKey: true });
    await flush(v);
    expect(placement(v).anchor).toBeCloseTo(scrolled, 6);
  });

  it('a drag selects in reading order across lines, to a line’s last unit past its end', async () => {
    const v = await open();
    const n = upl(v);
    const a = at(v, 5);
    const b = at(v, 2 * n + 3);
    expect(v.input.down(a.x, a.y)).toBe(true);
    v.input.move(b.x, b.y);
    expect(selected(v)).toEqual([5, 2 * n + 4]);
    // past the second line's end (its number gutter, the space after it): the line's last unit
    const s = placement(v).strips;
    const line = s.axes[1];
    v.input.move(line.x0 + line.width + 6, b.y - (b.y - a.y) / 2);
    expect(selected(v)).toEqual([5, 2 * n]);
    v.input.up();
  });

  it('↑↓ move the cursor a line, Shift+↑↓ extend by a line; ↑ on the first line has no section above', async () => {
    const v = await open();
    const n = upl(v);
    const p = at(v, 5);
    v.input.click(p.x, p.y);
    expect(v.selection.cursor.value?.position).toBe(5);
    expect(v.input.key({ key: 'ArrowDown' })).toBe(true);
    expect(v.selection.cursor.value?.position).toBe(5 + n);
    v.input.key({ key: 'ArrowUp' });
    expect(v.selection.cursor.value?.position).toBe(5);
    v.input.key({ key: 'ArrowUp' });
    expect(v.selection.cursor.value?.position).toBe(5);
    v.input.key({ key: 'ArrowDown', shiftKey: true });
    expect(selected(v)).toEqual([5, 5 + n + 1]);
  });

  it('a double-click on a letter reveals its line at the top, the letter size kept; past the box the cursor scrolls in', async () => {
    const v = await open();
    const n = upl(v);
    const before = placement(v).geometry!.cell;
    const box = placement(v).box;
    // the cursor stepped below the box: its line scrolls in, the least move
    v.input.click(at(v, 5).x, at(v, 5).y);
    for (let i = 0; i < 12; i++) v.input.key({ key: 'ArrowDown' });
    await flush(v, 400);
    const c = v.selection.cursor.value!.position;
    const shown = at(v, c);
    expect(shown.y).toBeGreaterThan(box.y);
    expect(shown.y).toBeLessThan(box.y + box.height);
    // a double-click focuses the letter: its line at the top, the size kept
    v.input.doubleClick(shown.x, shown.y);
    await flush(v, 400);
    expect(placement(v).anchor).toBe(Math.floor(c / n) * n);
    expect(placement(v).geometry!.cell).toBe(before);
  });

  it('a drag-select held past the box’s bottom scrolls it, by the clock, and extends; it stops at the button’s release', async () => {
    const v = await open();
    const box = placement(v).box;
    const a = at(v, 5);
    v.input.down(a.x, a.y);
    // 40 px past the bottom: 320 px a second
    v.input.move(a.x, box.y + box.height + 40);
    await flush(v, 250);
    const pitch = placement(v).geometry!.lineHeight + placement(v).geometry!.lineGap;
    const first = placement(v).anchor / upl(v);
    expect(first).toBeGreaterThan(0);
    await flush(v, 250);
    const second = placement(v).anchor / upl(v);
    // the scroll is the time's: about the same distance each quarter second (80 px)
    expect((second - first) * pitch).toBeGreaterThan(40);
    expect((second - first) * pitch).toBeLessThan(120);
    expect(selected(v)[1]).toBeGreaterThan(Math.floor(second) * upl(v));
    v.input.up();
    await flush(v, 250);
    expect(placement(v).anchor / upl(v)).toBeCloseTo(second, 6);
  });
});

describe('a text section and the rest of the view', () => {
  it('a selection made elsewhere scrolls the text to its first line (reveal.linked); one made in it does not', async () => {
    const v = await open('text-under-tracks');
    expect(placement(v).anchor).toBe(0);
    await v.execute('selection.add', { range: { system: System, start: 350, end: 360 } });
    await flush(v);
    const w = placement(v).window;
    expect(w.start).toBeLessThanOrEqual(350);
    expect(w.end).toBeGreaterThan(350);
    expect(v.view.history.value.back).toBe(0);
    // with the cursor in the text, a selection is made in it: the text stays where it is
    const top = placement(v).anchor;
    const p = at(v, Math.ceil(top / upl(v)) * upl(v) + 2);
    v.input.click(p.x, p.y);
    await v.execute('selection.add', { range: { system: System, start: 0, end: 5 } });
    await flush(v);
    expect(placement(v).anchor).toBe(top);
  });

  it('copies the selection’s letters as shown (selection.copyLetters), each run a line', async () => {
    const v = await open();
    await v.execute('selection.add', { range: { system: System, start: 10, end: 25 } });
    await v.execute('selection.add', { range: { system: System, start: 40, end: 44 } });
    expect(await v.execute('selection.copy-letters', {})).toBe(`${P53.slice(10, 25)}\n${P53.slice(40, 44)}`);
    expect(await v.execute('selection.copy-letters', { stored: true })).toBe(P53.slice(10, 44));
  });
});
