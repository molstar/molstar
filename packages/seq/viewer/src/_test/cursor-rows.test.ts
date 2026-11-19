/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard cursor's rows (grids-and-rows, The keyboard cursor; _spec/README.md "Behaviors", Cursor), headless on a
// fake clock: ↑ ↓ step display rows (the pinned band to the scrolled rows, a summary row one stop), Page Up and Page
// Down a page of rows with the band, Cmd/Ctrl+↑ ↓ the first and the last row, ← from the first column onto the header
// and → back, Shift+arrows the block, Space a cell or (on the header) the row's clade, Enter opens a summary row; a row
// never loaded is stepped over by its index; the cursor is drawn on its row, said by its row's label and carded there;
// on a grid ↑ ↓ step its y.

import type { Node } from '@molstar/core/repr/tree';
import { cellsSource, type FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { partOf } from '@molstar/seqx-core/query';
import { type RowsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { Announcer } from '../announcer.js';
import { HoverCardProvider } from '../behaviors/hover-card.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/v/s';
const Track = `${Section}/rows`;
const Columns = '/data/aln';
const Rows = '/data/aln/(rows)';
const Labels = ['HBA', 'HBB', 'MYG', 'GLB', 'CYG', 'NGB'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLSDG', 'MKAVLG', 'MEKVPG', 'MERPEP'];
const Newick = '((HBA,HBB)hb,(MYG,(GLB,(CYG,NGB)cn)gc)rest)root;';

/** Six globins in tree order (60 px: three rows a window), or `n` rows `R1…` with two pinned (200 px). */
function treeSpec(collapsed: string[] = []): Node {
  const b = createSeqSpec({ title: 'rows' });
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v
      .section('s', { of: aln })
      .track('rows', { label: '', height: 60 })
      .layer('r', { type: 'rows', data: aln, tree, order: 'tree', collapsed }),
  );
  return b.build();
}

function pinnedSpec(n = 40): Node {
  const labels = Array.from({ length: n }, (_, i) => `R${i + 1}`);
  const b = createSeqSpec({ title: 'pinned' });
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: labels, letters: labels.map(() => 'MVLSPA') },
  });
  b.views.view('v', (v) =>
    v
      .section('s', { of: aln })
      .track('rows', { label: '', height: 200 })
      .layer('r', { type: 'rows', data: aln, pinned: ['R1', 'R2'] } as never),
  );
  return b.build();
}

/** A lazy alignment of `n` rows (6 columns): its rows read by index, only those in view. */
const Lazy: FormatRuntime = {
  name: 'lazy-rows',
  formats: ['lazy-rows'],
  async open(source) {
    const n = (source.value as { rows: number }).rows;
    return cellsSource(n, 6, (start, end) => ({
      codes: new Uint8Array((end - start) * 6).fill(65),
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
    v.section('s', { of: aln }).track('rows', { label: '', height: 200 }).layer('r', { type: 'rows', data: aln }),
  );
  return b.build();
}

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames?.restore();
  frames = undefined;
  vi.useRealTimers();
});

async function make(spec: Node, lazy = false) {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  frames = installFrames();
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    ...(lazy ? { registry: [...(base.registry ?? []), { kind: 'format-runtime', runtime: Lazy }] } : {}),
    config: [[viewerItem('view.jumpMs'), 0]],
  } as never);
  viewer = v;
  await v.load(spec);
  v.resize(400);
  await v.settle();
  frames.flush();
  v.input.mac = false;
  return v;
}

const lay = (v: SeqViewer) => {
  for (const s of v.frame!.layout.sections)
    for (const band of s.lays) for (const l of band ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
};
const key = (v: SeqViewer, k: string, e: { shiftKey?: boolean; ctrlKey?: boolean } = {}) => {
  const ran = v.input.key({ key: k, ...e });
  frames!.flush();
  return ran;
};
const row = (v: SeqViewer) => v.selection.cursor.value?.yPosition;
/** The cursor pass's rects: x, y, w, h. */
const cursorRects = (v: SeqViewer) =>
  recordDrawLists(v.frame!.lists, ['cursor'])
    .split('\n')
    .filter((r) => r.startsWith('rect '))
    .map((r) => r.slice(5).split(' ').slice(0, 4).map(Number));

describe('the keyboard cursor on an alignment’s rows', () => {
  it('↑ ↓ step display rows from the first shown; the cursor is drawn on its row’s cell; ← onto the header, → back', async () => {
    const v = await make(treeSpec());
    // the first key places it: the first position and the first row shown
    expect(key(v, 'ArrowDown')).toBe(true);
    expect(v.selection.cursor.value).toMatchObject({
      section: Section,
      system: Columns,
      position: 0,
      ySystem: Rows,
      yPosition: 0,
      keyboard: true,
    });
    key(v, 'ArrowDown');
    key(v, 'ArrowDown');
    // tree order: HBA, HBB, MYG…
    expect(row(v)).toBe(2);
    key(v, 'ArrowUp');
    expect(row(v)).toBe(1);
    const l = lay(v);
    const [r] = cursorRects(v);
    expect(r[1]).toBeGreaterThanOrEqual(l.rowY(1));
    expect(r[1] + r[3]).toBeLessThanOrEqual(l.rowY(1) + l.rowScale);
    expect(r[0]).toBeGreaterThanOrEqual(l.box.x);
    // ← at the first column: the header (the labels' column, left of the data); → back to the cell
    key(v, 'ArrowLeft');
    expect(v.selection.cursor.value?.part).toBe('header');
    const [h] = cursorRects(v);
    expect(h[0] + h[2]).toBeLessThanOrEqual(l.box.x);
    key(v, 'ArrowLeft');
    expect(v.selection.cursor.value?.part).toBe('header');
    key(v, 'ArrowRight');
    expect(v.selection.cursor.value).toMatchObject({ position: 0, yPosition: 1 });
    expect(v.selection.cursor.value?.part ?? 'cell').toBe('cell');
    key(v, 'ArrowRight');
    expect(v.selection.cursor.value?.position).toBe(1);
  });

  it('a summary row is one stop (its first shown row); a step past it reveals the band’s next row', async () => {
    const v = await make(treeSpec(['hb']));
    expect(lay(v).display.count).toBe(5);
    key(v, 'ArrowDown');
    // the summary row of hb stands for HBA and HBB: the cursor on HBA
    expect(row(v)).toBe(0);
    key(v, 'ArrowDown');
    expect(row(v)).toBe(2);
    key(v, 'ArrowUp');
    expect(row(v)).toBe(0);
    key(v, 'ArrowDown');
    key(v, 'ArrowDown');
    key(v, 'ArrowDown');
    // display row 3 (CYG) is below the window of three rows: the band moves by one
    expect(row(v)).toBe(4);
    expect(v.view.band(Track)?.bandTop).toBe(1);
  });

  it('the pinned band steps to the scrolled rows; Page Down moves a page with the band and shows the row', async () => {
    const v = await make(pinnedSpec());
    const l = lay(v);
    expect(l.display.pinned).toBe(2);
    key(v, 'ArrowDown');
    expect(row(v)).toBe(0);
    key(v, 'ArrowDown');
    key(v, 'ArrowDown');
    // R1 and R2 in the band, then the scrolled rows from R1 in its place (the pinned rows are kept in place too)
    expect(v.selection.cursor.value?.yPosition).toBe(0);
    // the cursor on R1's place in the order is drawn there and in the band
    expect(cursorRects(v).length).toBe(2);
    // a page: 8 scrolled rows in the window less one
    expect(key(v, 'PageDown')).toBe(true);
    expect(row(v)).toBe(7);
    expect(v.view.band(Track)?.bandTop).toBe(7);
    const after = lay(v);
    const i = 2 + 7;
    expect(after.rowY(i)).toBeGreaterThanOrEqual(after.box.y + after.pinnedHeight);
    expect(after.rowY(i) + after.rowScale).toBeLessThanOrEqual(after.box.y + after.box.height);
    expect(key(v, 'PageUp')).toBe(true);
    expect(row(v)).toBe(0);
    expect(v.view.band(Track)?.bandTop).toBe(0);
  });

  it('Page Down counts a summary row one row; at the end it stays on the last', async () => {
    const v = await make(treeSpec(['gc']));
    // HBA, HBB, MYG, then GLB +2: four display rows, three a window, a page of two
    expect(lay(v).display.count).toBe(4);
    key(v, 'ArrowDown');
    key(v, 'PageDown');
    expect(row(v)).toBe(2);
    key(v, 'PageDown');
    // the summary row of gc, by its first row
    expect(row(v)).toBe(3);
    expect(v.view.band(Track)?.bandTop).toBe(1);
    key(v, 'PageDown');
    expect(row(v)).toBe(3);
  });

  it('Cmd/Ctrl+↑ ↓ go to the first and the last row, the band with them', async () => {
    const v = await make(pinnedSpec());
    key(v, 'ArrowDown');
    expect(key(v, 'ArrowDown', { ctrlKey: true })).toBe(true);
    expect(row(v)).toBe(39);
    const l = lay(v);
    expect(v.view.band(Track)?.bandTop).toBe(l.maxTop(l.rowScale, l.display.count));
    expect(key(v, 'ArrowUp', { ctrlKey: true })).toBe(true);
    // the first display row: R1 in the pinned band (the band keeps its scroll: the pinned band is always shown)
    expect(row(v)).toBe(0);
  });

  it('with no cursor shown Page Down is the rows’ scroll, Cmd/Ctrl+↓ the page’s', async () => {
    const v = await make(pinnedSpec());
    expect(key(v, 'ArrowDown', { ctrlKey: true })).toBe(false);
    expect(v.selection.cursor.value).toBeUndefined();
    expect(key(v, 'PageDown')).toBe(true);
    expect(v.view.band(Track)?.bandTop).toBe(7);
    expect(v.selection.cursor.value).toBeUndefined();
  });

  it('Space (on its release) toggles the cell in the selection; on the header the row’s clade folds and unfolds', async () => {
    const v = await make(treeSpec());
    key(v, 'ArrowDown');
    key(v, 'ArrowRight');
    const space = () => {
      v.input.setSpace(true);
      const ran = v.input.spaceUp();
      frames!.flush();
      return ran;
    };
    expect(space()).toBe(true);
    const rects = v.selection.value.value.location.rects;
    expect(rects.length).toBe(1);
    const sides = [rects[0].a, rects[0].b];
    expect([...sides.find((s) => s.system === Columns)!.ranges]).toEqual([1, 2]);
    expect([...sides.find((s) => s.system === Rows)!.ranges]).toEqual([0, 1]);
    space();
    expect(v.selection.value.value.location.rects.length).toBe(0);
    // a Space held through a drag pans: no toggle
    v.input.setSpace(true);
    v.input.down(lay(v).box.x + 5, lay(v).rowY(0) + 3, {});
    v.input.move(lay(v).box.x + 25, lay(v).rowY(0) + 3);
    v.input.up();
    expect(v.input.spaceUp()).toBe(false);
    // the header: HBA's parent clade (hb) folds, the cursor on its summary row; again unfolds
    key(v, 'ArrowLeft');
    key(v, 'ArrowLeft');
    expect(v.selection.cursor.value?.part).toBe('header');
    space();
    expect(v.view.band(Track)?.folds).toEqual(['hb']);
    expect(lay(v).display.displayRow(row(v)!)).toBe(0);
    space();
    expect(v.view.band(Track)?.folds).toEqual([]);
  });

  it('Shift+arrows extend the block from where it began; Enter opens a summary row', async () => {
    const v = await make(treeSpec(['cn']));
    key(v, 'ArrowDown');
    key(v, 'ArrowDown', { shiftKey: true });
    key(v, 'ArrowRight', { shiftKey: true });
    key(v, 'ArrowDown', { shiftKey: true });
    const r = v.selection.value.value.location.rects[0];
    const sides = [r.a, r.b];
    expect([...sides.find((s) => s.system === Columns)!.ranges]).toEqual([0, 2]);
    expect([...sides.find((s) => s.system === Rows)!.ranges]).toEqual([0, 3]);
    // Enter on the summary row of cn (CYG, NGB) opens it
    key(v, 'ArrowDown', { ctrlKey: true });
    expect(lay(v).display.rowAt(lay(v).display.displayRow(row(v)!))).toBeLessThan(0);
    expect(key(v, 'Enter')).toBe(true);
    expect(v.view.band(Track)?.folds).toEqual([]);
    expect(partOf(v.selection.value.value.location, Rows)).toBeUndefined();
  });

  it('on the header Shift+↓ extends a run of rows and Enter focuses the row’s letters', async () => {
    const b = createSeqSpec({ title: 'header' });
    const L = 'MVLSPA'.repeat(5);
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: { label: ['A', 'B', 'C', 'D'], letters: [`${'-'.repeat(20)}MVLSPAMVLS`, L, L, L] },
    });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }).track('rows', { label: '', height: 60 }).layer('r', { type: 'rows', data: aln }),
    );
    const v = await make(b.build());
    key(v, 'ArrowDown');
    key(v, 'ArrowLeft');
    expect(v.selection.cursor.value?.part).toBe('header');
    expect(v.runCommand('cursor.focus', {})).toBe('Focused');
    expect(v.view.get(Section)).toMatchObject({ start: 20, end: 30 });
    key(v, 'ArrowDown', { shiftKey: true });
    key(v, 'ArrowDown', { shiftKey: true });
    expect([...partOf(v.selection.value.value.location, Rows)!.ranges]).toEqual([0, 3]);
    expect(v.selection.cursor.value).toMatchObject({ part: 'header', yPosition: 2 });
  });

  it('the announcement and the card at the cursor name its row', async () => {
    const v = await make(treeSpec(['hb']));
    const a = new Announcer(v);
    key(v, 'ArrowDown');
    expect(a.cursorText()).toContain('row HBA +1');
    key(v, 'ArrowDown');
    expect(a.cursorText()).toContain('row MYG');
    const card = v.behavior(HoverCardProvider)!.shown.value!;
    expect(card.card.title).toBe('MYG');
    const l = lay(v);
    expect(card.place.anchor.y).toBe(l.rowY(1));
    expect(card.place.anchor.height).toBe(l.rowScale);
    key(v, 'ArrowLeft');
    expect(a.cursorText()).toBe('row MYG');
    a.dispose();
  });
});

describe('the keyboard cursor on rows never loaded', () => {
  it('steps by index over rows never read: the last of 100,000 rows, said pending until its tile comes', async () => {
    const v = await make(lazySpec(100_000), true);
    expect(lay(v).total).toBe(100_000);
    key(v, 'ArrowDown');
    expect(row(v)).toBe(0);
    key(v, 'ArrowDown', { ctrlKey: true });
    expect(row(v)).toBe(99_999);
    const a = new Announcer(v);
    key(v, 'ArrowUp');
    expect(row(v)).toBe(99_998);
    // the band at its end shows it: its tile asked for and read
    await v.settle();
    frames!.flush();
    expect(a.cursorText()).toContain('row L99998');
    key(v, 'PageUp');
    expect(row(v)).toBe(99_998 - (Math.floor(200 / lay(v).rowScale) - 1));
    a.dispose();
  });
});

describe('the keyboard cursor on a grid', () => {
  function gridSpec() {
    const b = createSeqSpec({ title: 'grid y' });
    const x = b.data.sequence('x', { alphabet: 'rna', letters: 'ACGU'.repeat(15) });
    const y = b.data.sequence('y', { alphabet: 'rna', letters: 'GUCA'.repeat(10) });
    const g = b.data.grid('g', { x, y, values: Float32Array.from({ length: 60 * 40 }, (_, i) => i) });
    b.views.view('v', (v) => {
      v.section('plot', { layout: 'grid', of: x, of_y: y }).track('plot').layer('grid', { type: 'grid', data: g });
    });
    return b.build();
  }

  it('↑ ↓ step its y, Cmd/Ctrl+↑ ↓ its ends; the cell is drawn and said with both positions', async () => {
    const v = await make(gridSpec());
    key(v, 'ArrowDown');
    expect(v.selection.cursor.value).toMatchObject({
      system: '/data/x',
      position: 0,
      ySystem: '/data/y',
      yPosition: 0,
    });
    key(v, 'ArrowDown');
    key(v, 'ArrowDown');
    expect(row(v)).toBe(2);
    key(v, 'ArrowUp');
    expect(row(v)).toBe(1);
    expect(cursorRects(v).length).toBe(1);
    key(v, 'ArrowDown', { ctrlKey: true });
    expect(row(v)).toBe(39);
    key(v, 'ArrowUp', { ctrlKey: true });
    expect(row(v)).toBe(0);
    const a = new Announcer(v);
    expect(a.cursorText()).toContain('× y');
    a.dispose();
  });
});
