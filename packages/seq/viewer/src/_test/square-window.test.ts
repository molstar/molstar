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
