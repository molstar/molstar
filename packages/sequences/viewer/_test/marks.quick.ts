/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of marks.heap.ts (`test:perf`): the marks' and the projector's per-move updates (marks-cases.ts, the
// guard's cases), each with its control.

import { beforeAll, describe, expect, it } from 'vitest';
import { MarkSet } from '../src/marks.js';
import { quickCase } from './heap.js';
import { marksCases, sorting } from './marks-cases.js';

let c: Awaited<ReturnType<typeof marksCases>>;

beforeAll(async () => {
  c = await marksCases();
}, 60_000);

describe('marks and the projector, quick', () => {
  it('a hover closed and projected into every section', ({ task }) => {
    quickCase(task, 'marks: hover', c.hover());
    expect(c.sink()).not.toBe(0);
  });

  it('loci converted and projected', ({ task }) => {
    quickCase(task, 'marks: loci', c.loci());
  });

  it('Coordinates through a registered translator', ({ task }) => {
    quickCase(task, 'marks: coordinates', c.coordinates());
  });

  it('a MarkSet sorted and merged at close', ({ task }) => {
    quickCase(task, 'marks: sort', sorting(new MarkSet()));
  });

  it("row rectangles and a ligand's atoms", ({ task }) => {
    quickCase(task, 'marks: row rects', c.rowRects());
  });
});
