/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { Op } from '@molstar/seqx-core/render';
import { SquareNavigator } from '@molstar/seqx-core/render/navigator-square';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import type { NavigatorModel } from '../navigator.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The square over a ring is a window (_spec/README.md "Composition", Navigators, The square over a set; Decided 16 of
// the design of the genome ring): one click on any cell shows that pair at once, the window is the blocks shown (a run one box, a
// scattered set its cells), its box dragged, resized and keyed, each gesture one patch and one undo step.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Ring = '/views/main/ring';
/** Eight chromosomes, equal lengths: each cell of the 200 px square's lattice is the same size. */
const Names = Array.from({ length: 8 }, (_, i) => `chr${i + 1}`);

function spec(include?: readonly string[]): Node {
  const b = createSeqSpec({ title: 'the square over a ring' });
  const g = b.data.assembly(
    'g',
    assembly(
      'Eight',
      Names.map((name) => ({ name, length: 1_000_000 })),
    ),
  );
  b.views.view('main', (v) => {
    genome(v, 'ring', g, { layout: 'rings', ...(include ? { include } : {}) } as never);
    v.navigator('square', { of: [Ring, Ring], axes: 'both', system: g, switch: 'include', size: [200, 200] } as never);
  });
  return b.build();
}

async function make(include?: readonly string[]) {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, history: uiHistory() });
  viewers.push(viewer);
  await viewer.load(spec(include));
  viewer.resize(600, 600);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  const n = viewer.navigatorModel('square')!;
  n.resize(200, 200);
  n.draw();
  return { v: viewer, n };
}

async function settle(v: SeqViewer, n: NavigatorModel) {
  await v.settle();
  frames?.flush();
  v.draw();
  n.draw();
}

const square = (n: NavigatorModel) => (n as unknown as { renderer: SquareNavigator }).renderer;
const include = (v: SeqViewer) => {
  const views = v.spec.root.value.children!.find((c) => c.kind === 'views')!;
  return views.children![0].children!.find((c) => c.id === 'ring')!.params?.include as string[] | undefined;
};
/** A cell's centre in the navigator's px (block indices along x and y). */
const cell = (n: NavigatorModel, i: number, j: number) => {
  const e = edge(n);
  return [(e(i, true) + e(i + 1, true)) / 2, (e(j, false) + e(j + 1, false)) / 2] as const;
};
/** A lattice edge's px (block `k`'s start; `k` the count: the end), read from the square's lattice. */
const edge = (n: NavigatorModel) => {
  const s = square(n) as unknown as {
    at(u: number, horizontal: boolean): number;
    lattice: { blocks: { offset?: number }[] };
    total: number;
  };
  const blocks = s.lattice.blocks;
  return (k: number, h: boolean) => s.at(k < blocks.length ? (blocks[k].offset ?? 0) : s.total, h);
};
/** The window pass's stroked rects (the outline's), as `[x0, y0, x1, y1]` rounded. */
const outline = (n: NavigatorModel) => {
  const pass = n.frame!.lists.passes.window;
  const r = pass.rects;
  const out: number[][] = [];
  let k = 0;
  for (let o = 0; o < pass.opCount; o++) {
    if (pass.ops[o] !== Op.rect) continue;
    const i = k++;
    if (r.strokeWidth[i] > 0) out.push([r.x[i], r.y[i], r.x[i] + r.w[i], r.y[i] + r.h[i]].map(Math.round));
  }
  return out;
};
/** The outline a set of blocks should have: one box for a run, its cells for a scattered set. */
const cellsOf = (n: NavigatorModel, blocks: readonly number[]) => {
  const e = edge(n);
  const sorted = [...blocks].sort((a, b) => a - b);
  const run = sorted.every((b, k) => k === 0 || b === sorted[k - 1] + 1);
  const box = (i0: number, i1: number, j0: number, j1: number) =>
    [e(i0, true), e(j0, false), e(i1, true), e(j1, false)].map(Math.round);
  if (run) return [box(sorted[0], sorted.at(-1)! + 1, sorted[0], sorted.at(-1)! + 1)];
  return sorted.flatMap((i) => sorted.map((j) => box(i, i + 1, j, j + 1)));
};
const click = (n: NavigatorModel, x: number, y: number, detail = 1) => {
  n.down(x, y);
  n.up();
  return n.click(x, y, detail);
};
const undo = (v: SeqViewer) => (v.spec.history as UiHistory).stack.undo();

describe('The square over a ring: one click shows a pair', () => {
  it('a click on any cell, inside the window or out, shows that pair at once and outlines exactly it', async () => {
    const { v, n } = await make();
    // the whole genome shown: the window is every cell, one box
    expect(outline(n)).toEqual(cellsOf(n, [0, 1, 2, 3, 4, 5, 6, 7]));
    // inside the window (every cell is): chr2 × chr6, a pair off the diagonal
    expect(click(n, ...cell(n, 1, 5))).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr2', 'chr6']);
    expect(outline(n)).toEqual(cellsOf(n, [1, 5]));
    // outside the window now: chr3 × chr4, a run on the diagonal
    expect(click(n, ...cell(n, 3, 2))).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr3', 'chr4']);
    expect(outline(n)).toEqual(cellsOf(n, [2, 3]));
    // a quick second click (the platform's count 2) on another cell is a click: never a select first
    expect(click(n, ...cell(n, 6, 6), 2)).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr7']);
    expect(outline(n)).toEqual(cellsOf(n, [6]));
    // and its double-click adds no zoom
    expect(n.doubleClick(...cell(n, 6, 6))).toBe(false);
    // each a step: undone back to the whole
    await undo(v);
    await undo(v);
    await undo(v);
    await settle(v, n);
    expect(include(v)).toBeUndefined();
  });
});

const drag = (n: NavigatorModel, from: readonly [number, number], to: readonly [number, number]) => {
  expect(n.down(...from)).toBe(true);
  const steps = 4;
  for (let s = 1; s <= steps; s++)
    n.move(from[0] + ((to[0] - from[0]) * s) / steps, from[1] + ((to[1] - from[1]) * s) / steps);
  const cursor = n.cursor;
  n.up();
  expect(n.click(...to)).toBe(false);
  return cursor;
};
/** The ring's shown chromosomes, as the resolved view has them. */
const shown = (v: SeqViewer) =>
  v.view.resolved
    .value!.sections.find((s) => s.address === Ring)!
    .blocks!.map((b) => v.spec.model.value.systems.get(b.axis.system)!.name);
const toggleClick = (n: NavigatorModel, at: readonly [number, number], mac = false) => {
  n.down(...at);
  n.up();
  return n.click(...at, 1, mac ? { metaKey: true } : { ctrlKey: true });
};

describe('The square over a ring: its box a window', () => {
  it('dragging the box moves the run along the diagonal, snapped to whole blocks; one patch, one undo step', async () => {
    const { v, n } = await make(['chr2', 'chr3']);
    const e = edge(n);
    const mid = [e(2, true), e(2, false)] as const;
    // the cursor over the body says a move
    n.move(...mid);
    expect(n.cursor).toBe('move');
    // three blocks down the diagonal
    const to = [mid[0] + 3 * (e(1, true) - e(0, true)), mid[1] + 3 * (e(1, false) - e(0, false))] as const;
    expect(drag(n, mid, to)).toBe('move');
    await settle(v, n);
    expect(include(v)).toEqual(['chr5', 'chr6']);
    expect(shown(v)).toEqual(['chr5', 'chr6']);
    expect(outline(n)[0]).toEqual(cellsOf(n, [4, 5])[0]);
    // past the end it stops at the last block
    const far = [e(5, true), e(5, false)] as const;
    drag(n, far, [far[0] + 500, far[1] + 500]);
    await settle(v, n);
    expect(include(v)).toEqual(['chr7', 'chr8']);
    await undo(v);
    await settle(v, n);
    expect(include(v)).toEqual(['chr5', 'chr6']);
    await undo(v);
    await settle(v, n);
    expect(include(v)).toEqual(['chr2', 'chr3']);
  });

  it('an edge or a corner grows or shrinks the run; the cursor says which; Escape cancels the preview', async () => {
    const { v, n } = await make(['chr3', 'chr4']);
    const e = edge(n);
    const right = [e(4, true), (e(2, false) + e(4, false)) / 2] as const;
    n.move(...right);
    expect(n.cursor).toBe('ew-resize');
    // the right edge two blocks on: chr3 to chr6
    expect(drag(n, right, [e(6, true) - 1, right[1]])).toBe('ew-resize');
    await settle(v, n);
    expect(include(v)).toEqual(['chr3', 'chr4', 'chr5', 'chr6']);
    // the top edge down one: chr4 to chr6
    const top = [(e(2, true) + e(6, true)) / 2, e(2, false)] as const;
    n.move(...top);
    expect(n.cursor).toBe('ns-resize');
    drag(n, top, [top[0], e(3, false) + 1]);
    await settle(v, n);
    expect(include(v)).toEqual(['chr4', 'chr5', 'chr6']);
    // the top-left corner along the diagonal: from chr1
    const corner = [e(3, true), e(3, false)] as const;
    n.move(...corner);
    expect(n.cursor).toBe('nwse-resize');
    drag(n, corner, [e(0, true), e(0, false)]);
    await settle(v, n);
    expect(include(v)).toEqual(['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6']);
    // a resize never leaves fewer than one block
    const end = [e(6, true), e(6, false)] as const;
    drag(n, end, [e(0, true), e(0, false)]);
    await settle(v, n);
    expect(include(v)).toEqual(['chr1']);
    // Escape during a drag: the preview gone, nothing written
    const box = cell(n, 0, 0);
    expect(n.down(...box)).toBe(true);
    n.move(box[0] + 60, box[1] + 60);
    n.draw();
    expect(outline(n)[0]).not.toEqual(cellsOf(n, [0])[0]);
    expect(n.cancel()).toBe(true);
    n.draw();
    expect(outline(n)[0]).toEqual(cellsOf(n, [0])[0]);
    expect(include(v)).toEqual(['chr1']);
    // one step a gesture
    for (const want of [
      ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6'],
      ['chr4', 'chr5', 'chr6'],
    ]) {
      await undo(v);
      await settle(v, n);
      expect(include(v)).toEqual(want);
    }
  });

  it('a drag to every block writes no include; Home shows all', async () => {
    const { v, n } = await make(['chr1', 'chr2']);
    const e = edge(n);
    drag(n, [e(2, true), e(2, false)], [e(8, true) + 4, e(8, false) + 4]);
    await settle(v, n);
    expect(include(v)).toBeUndefined();
    expect(shown(v)).toEqual(Names);
    await undo(v);
    await settle(v, n);
    expect(include(v)).toEqual(['chr1', 'chr2']);
    expect(n.key('Home')).toBe(true);
    await settle(v, n);
    expect(include(v)).toBeUndefined();
  });

  it('the keyboard: arrows move the run, Shift+arrows resize it, Enter shows the keyboard’s cell', async () => {
    const { v, n } = await make(['chr2', 'chr3']);
    n.setFocus(true);
    n.draw();
    expect(n.key('ArrowRight')).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr3', 'chr4']);
    expect(n.key('ArrowDown', { shiftKey: true })).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr3', 'chr4', 'chr5']);
    expect(n.key('ArrowLeft', { shiftKey: true })).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr3', 'chr4']);
    // the keyboard's cell: on the run's first diagonal cell, moved by Ctrl+arrows, said; Enter shows its pair
    expect(n.key('ArrowRight', { ctrlKey: true })).toBe(true);
    expect(n.key('ArrowRight', { ctrlKey: true })).toBe(true);
    expect(n.readout.value).toBe('chr5 × chr3');
    expect(n.key('Enter')).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr3', 'chr5']);
    // at the lattice's start the run stops
    await v.execute('section.include', { address: Ring, include: ['chr1'] });
    await settle(v, n);
    expect(n.key('ArrowLeft')).toBe(false);
  });
});

describe('The square over a ring: any set', () => {
  it('Cmd or Ctrl+click on a cell adds or removes its two chromosomes; a scattered set is outlined cell by cell', async () => {
    const { v, n } = await make(['chr4']);
    // chr4 and chr8 (the author's chr4 and chr9 on a smaller assembly), added by Cmd+click on their cell
    expect(toggleClick(n, cell(n, 3, 7), true)).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr4', 'chr8']);
    expect(shown(v)).toEqual(['chr4', 'chr8']);
    // the outline: each one's diagonal cell and the two cross cells; no grips (no filled squares), no box
    expect(outline(n)).toEqual(cellsOf(n, [3, 7]));
    const rects = n.frame!.lists.passes.window.rects;
    let fills = 0;
    for (let i = 0; i < rects.n; i++) if (rects.strokeWidth[i] === 0) fills++;
    expect(fills).toBe(0);
    // over a scattered set a press drags nothing: the cursor a click's
    n.move(...cell(n, 3, 7));
    expect(n.cursor).toBe('pointer');
    // Ctrl+click on chr1 × chr4: chr1 added; again: both removed
    expect(toggleClick(n, cell(n, 0, 3))).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr1', 'chr4', 'chr8']);
    expect(toggleClick(n, cell(n, 0, 3))).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr8']);
    // removing the last one leaves nothing: refused
    expect(toggleClick(n, cell(n, 7, 7))).toBe(false);
    expect(include(v)).toEqual(['chr8']);
  });

  it('a label: a click shows only it, Cmd or Ctrl+click adds or removes it; each one undo step', async () => {
    const { v, n } = await make();
    const e = edge(n);
    // chr5's name above its column, chr2's left of its row
    const top = [(e(4, true) + e(5, true)) / 2, e(0, false) - 4] as const;
    const left = [e(0, true) - 4, (e(1, false) + e(2, false)) / 2] as const;
    n.move(...top);
    expect(n.readout.value).toBe('chr5');
    expect(n.down(...top)).toBe(true);
    n.up();
    expect(n.click(...top)).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr5']);
    expect(toggleClick(n, left, true)).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr2', 'chr5']);
    expect(outline(n)).toEqual(cellsOf(n, [1, 4]));
    expect(toggleClick(n, top)).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr2']);
    for (const want of [['chr2', 'chr5'], ['chr5'], undefined]) {
      await undo(v);
      await settle(v, n);
      expect(include(v)).toEqual(want);
    }
  });

  it('the keyboard adds and removes: Ctrl+Enter on the keyboard’s cell; plain arrows move the cell over a scattered set', async () => {
    const { v, n } = await make(['chr1', 'chr3']);
    n.setFocus(true);
    n.draw();
    // scattered: arrows move the keyboard's cell (from the first shown block), Shift+arrows nothing
    expect(n.key('ArrowDown')).toBe(true);
    expect(n.key('ArrowDown')).toBe(true);
    expect(n.key('ArrowDown')).toBe(true);
    expect(n.readout.value).toBe('chr1 × chr4');
    expect(n.key('ArrowDown', { shiftKey: true })).toBe(false);
    expect(n.key('Enter', { ctrlKey: true })).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr1', 'chr3', 'chr4']);
    expect(n.key('Home')).toBe(true);
    await settle(v, n);
    expect(include(v)).toBeUndefined();
  });

  it('the outline is only the cells the ring shows: no stray cell on chr1 × chr1, focused or not', async () => {
    const { v, n } = await make();
    // the part focused while the whole genome is shown (its keyboard's cell would be chr1 × chr1)
    n.setFocus(true);
    n.draw();
    expect(click(n, ...cell(n, 4, 6))).toBe(true);
    await settle(v, n);
    expect(include(v)).toEqual(['chr5', 'chr7']);
    // a pointer's click: the outline is the set's cells, nothing more
    expect(outline(n)).toEqual(cellsOf(n, [4, 6]));
    n.setFocus(false);
    n.draw();
    expect(outline(n)).toEqual(cellsOf(n, [4, 6]));
    // a key drives the cell: drawn inset and rounded, on a cell of the set (never chr1 × chr1)
    n.setFocus(true);
    expect(n.key('ArrowRight', { ctrlKey: true })).toBe(true);
    n.draw();
    const e = edge(n);
    const extra = outline(n).filter((r) => !cellsOf(n, [4, 6]).some((c) => c.join() === r.join()));
    expect(extra).toHaveLength(1);
    expect(extra[0][0]).toBeGreaterThan(Math.round(e(4, true)));
    expect(extra[0][1]).toBeGreaterThan(Math.round(e(4, false)));
    expect(n.frame!.lists.passes.window.rects.radius[n.frame!.lists.passes.window.rects.n - 1]).toBeGreaterThan(0);
  });

  it('EC-110: a ring’s view zoomed (its window idle unless followed): the square’s window stays its set', async () => {
    const { v, n } = await make(['chr2', 'chr3']);
    expect(outline(n)).toEqual(cellsOf(n, [1, 2]));
    v.view.jump(() => v.view.zoom(Ring, 4), true, 'zoom');
    await settle(v, n);
    expect(v.view.get(Ring)!.end - v.view.get(Ring)!.start).toBeLessThan(2_000_000);
    // the ring draws its set whole: the outline is the set, the run's grips kept
    expect(outline(n)).toEqual(cellsOf(n, [1, 2]));
    n.move(...cell(n, 1, 2));
    expect(n.cursor).toBe('move');
  });
});
