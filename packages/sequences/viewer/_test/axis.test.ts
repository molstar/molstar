/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section axis (VC3; seqspec:model#pieces-and-layout-maps): every index names its coordinate system, `k` the
// layout's shown index and `i` the block coordinate system's. A ranged unit at every index, through a filter, case 10's
// abl1-kd; reversed pieces and a cyclic origin through the run converters; breaks per filter run; the block gap and
// equal axes.

import { describe, expect, it } from 'vitest';
import { case10 } from '../../core/_test/format/cases.js';
import { Axis } from '../src/axis.js';
import {
  axisOf,
  blocksDoc,
  k,
  type Resolved,
  rangedDoc,
  readBack,
  resolved,
  setPiece,
  toRuns,
  toShown,
} from './axis-fixtures.js';

const range = (a: number, b: number) => Array.from({ length: b - a }, (_, i) => a + i);
const shownOver = (domain: number, f: (i: number) => number) => Array.from({ length: domain }, (_, i) => f(i));

describe("a ranged unit's block coordinate system (layout k ↔ block i), at every index", () => {
  it('a unit over [5, 15): shown k is block index 5 + k, every block index outside the piece is -1', async () => {
    const ax = axisOf(await resolved(rangedDoc(), 'v'));
    expect([ax.base[0], ax.first[0], ax.dir[0], ax.shown[0]]).toEqual([5, 5, 1, 10]);
    const r = readBack(ax, 0, 20);
    expect(r.index).toEqual(range(5, 15));
    expect(r.shown).toEqual(shownOver(20, (i) => (i >= 5 && i < 15 ? i - 5 : -1)));
    expect(r.breaks).toEqual([]);
    // Block-system runs to shown runs, clipped to the piece, and back.
    expect(toShown(ax, 0, 0, 20)).toEqual([[0, 10]]);
    expect(toShown(ax, 0, 7, 10)).toEqual([[2, 5]]);
    expect(toShown(ax, 0, 0, 5)).toEqual([]);
    expect(toRuns(ax, 0, 0, 10)).toEqual([[5, 15]]);
    expect(toRuns(ax, 0, 3, 4)).toEqual([[8, 9]]);
  });

  it("through a filter: the shown positions are the filter's within the piece, a break between its runs", async () => {
    const ax = axisOf(await resolved(rangedDoc(), 'f'));
    expect(ax.shown[0]).toBe(5);
    const r = readBack(ax, 0, 20);
    expect(r.index).toEqual([7, 8, 9, 12, 13]);
    expect(r.shown).toEqual(shownOver(20, (i) => [7, 8, 9, 12, 13].indexOf(i)));
    expect(r.breaks).toEqual([3]);
    expect(toShown(ax, 0, 8, 13)).toEqual([[1, 4]]);
    expect(toShown(ax, 0, 10, 12)).toEqual([]);
    expect(toRuns(ax, 0, 0, 5)).toEqual([
      [7, 10],
      [12, 14],
    ]);
    expect(toRuns(ax, 0, 2, 4)).toEqual([
      [9, 10],
      [12, 13],
    ]);
  });

  it("case 10's abl1-kd [241, 493) at every index; T315 (block index 314) is shown at 73", async () => {
    const ax = axisOf(await resolved(case10(), 'text'));
    const r = readBack(ax, 0, 600);
    expect(r.index).toEqual(range(241, 493));
    expect(r.shown).toEqual(shownOver(600, (i) => (i >= 241 && i < 493 ? i - 241 : -1)));
    expect([ax.shownOf(0, 314), ax.indexOf(0, k(73))]).toEqual([73, 314]);
  });
});

describe("pieces the format gains in AX1 and CY1, through the axis (layoutOf's maps)", () => {
  async function whole(edit: (r: Resolved) => void) {
    const r = await resolved(rangedDoc(), 'w');
    edit(r);
    return axisOf(r);
  }

  it('a reversed piece (dir -1) over [3, 8): shown k is block index 7 - k; runs map with their order flipped', async () => {
    const ax = await whole((r) => setPiece(r.view, 0, 3, 8, -1, 20));
    expect([ax.first[0], ax.dir[0], ax.shown[0]]).toEqual([7, -1, 5]);
    const r = readBack(ax, 0, 20);
    expect(r.index).toEqual([7, 6, 5, 4, 3]);
    expect(r.shown).toEqual(shownOver(20, (i) => (i >= 3 && i < 8 ? 7 - i : -1)));
    expect(toShown(ax, 0, 4, 6)).toEqual([[2, 4]]);
    expect(toShown(ax, 0, 0, 20)).toEqual([[0, 5]]);
    expect(toRuns(ax, 0, 1, 3)).toEqual([[5, 7]]);
  });

  it("a reversed piece through a filter: the filter's positions in the piece's order, breaks where they skip", async () => {
    const r = await resolved(rangedDoc(), 'f');
    setPiece(r.view, 0, 5, 15, -1, 20);
    const ax = axisOf(r);
    const b = readBack(ax, 0, 20);
    expect(b.index).toEqual([13, 12, 9, 8, 7]);
    expect(b.shown).toEqual(shownOver(20, (i) => [13, 12, 9, 8, 7].indexOf(i)));
    expect(b.breaks).toEqual([2]);
    expect(toShown(ax, 0, 8, 13)).toEqual([[1, 4]]);
    // In shown order: the higher run first.
    expect(toRuns(ax, 0, 0, 5)).toEqual([
      [12, 14],
      [7, 10],
    ]);
  });

  it('across a cyclic origin: [15, 20) then [0, 5); two shown runs for a run on both sides, one where they meet', async () => {
    const ax = await whole((r) => setPiece(r.view, 0, 15, 5, 1, 20));
    expect([ax.first[0], ax.shown[0], ax.wrapAt(0)]).toEqual([15, 10, 5]);
    const r = readBack(ax, 0, 20);
    expect(r.index).toEqual([15, 16, 17, 18, 19, 0, 1, 2, 3, 4]);
    expect(r.shown).toEqual(shownOver(20, (i) => (i >= 15 ? i - 15 : i < 5 ? i + 5 : -1)));
    expect(r.breaks).toEqual([]);
    expect(toShown(ax, 0, 3, 17)).toEqual([
      [0, 2],
      [8, 10],
    ]);
    expect(toShown(ax, 0, 0, 20)).toEqual([[0, 10]]);
    expect(toRuns(ax, 0, 3, 7)).toEqual([
      [18, 20],
      [0, 2],
    ]);
  });

  it('a cyclic origin reversed, and a whole cyclic system from an origin', async () => {
    const rev = await whole((r) => setPiece(r.view, 0, 15, 5, -1, 20));
    expect(readBack(rev, 0, 20).index).toEqual([4, 3, 2, 1, 0, 19, 18, 17, 16, 15]);
    expect(toShown(rev, 0, 3, 17)).toEqual([
      [0, 2],
      [8, 10],
    ]);
    const all = await whole((r) => setPiece(r.view, 0, 4, 4, 1, 20));
    expect(readBack(all, 0, 20).index).toEqual([...range(4, 20), ...range(0, 4)]);
    expect(toShown(all, 0, 0, 20)).toEqual([[0, 20]]);
  });

  it('a cyclic origin through a filter: no break where both sides reach the origin, one where they skip', async () => {
    const r = await resolved(rangedDoc(), 'f');
    // The filter [7, 10) ∪ [12, 14) on a piece from 12 around to 10 (the unit's [5, 15) as a domain of 15).
    setPiece(r.view, 0, 12, 10, 1, 15);
    const ax = axisOf(r);
    const b = readBack(ax, 0, 15);
    expect(b.index).toEqual([12, 13, 7, 8, 9]);
    expect(b.breaks).toEqual([2]);
    setPiece(r.view, 0, 9, 8, 1, 15);
    const meet = readBack(axisOf(r), 0, 15);
    expect(meet.index).toEqual([9, 12, 13, 7]);
    expect(meet.breaks).toEqual([1, 3]);
  });
});

describe('blocks side by side', () => {
  it("24 blocks share 7 gaps' room, at least 4 px each, never more than blockGap", async () => {
    const r = await resolved(blocksDoc(24), 'v');
    expect(new Axis(r.view, r.set, 0, 28, 3).gap).toBe(9);
    expect(new Axis(r.view, r.set, 0, 2, 3).gap).toBe(2);
    const few = await resolved(blocksDoc(3), 'v');
    const ax = axisOf(few);
    expect([ax.gap, ax.cells, Array.from(ax.offset)]).toEqual([28, 30, [0, 10, 20]]);
  });

  it('equal axes compare blocks, pieces and shown positions', async () => {
    const r = await resolved(rangedDoc(), 'w');
    const a = axisOf(r);
    expect(a.equals(axisOf(r))).toBe(true);
    setPiece(r.view, 0, 0, 20, -1, 20);
    expect(a.equals(axisOf(r))).toBe(false);
    const f = await resolved(rangedDoc(), 'f');
    const v = await resolved(rangedDoc(), 'v');
    expect(axisOf(f).equals(axisOf(v))).toBe(false);
  });
});
