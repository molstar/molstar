/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Interval, SortedArray } from '@molstar/core/data/int';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { CopyMarkerBlock } from '@molstar/graphics/gl/culled-data';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { createCopyMarkers, DefaultCopyBlockBudget, getCopyBlockBudget } from '../target/copy-markers.js';
import { StorePageSize } from '../target/page-store.js';

/** two templates: the first on store pages 2 and 3 (two copies), the second on page 5 (three copies) */
function create() {
  const m = createCopyMarkers();
  m.layout([
    { pages: Int32Array.of(2, 3), copyCount: 2 },
    { pages: Int32Array.of(5), copyCount: 3 },
  ]);
  return m;
}

/** a marker as the vertex stage reads it (`sCopyMarkers`: the words, then the blocks' bytes) */
function read(m: ReturnType<typeof create>, c: number, index: number) {
  const word = m.words[2 * c];
  if (!(word & CopyMarkerBlock)) return word & 0xff;
  const b = (m.words[2 * c + 1] | 0) + index - m.words.byteLength;
  return m.bytes[b];
}

describe('copy markers', () => {
  it('marks a copy whole by its word, no block', () => {
    const m = create();
    expect(m.copyCount).toBe(5);
    expect(m.markCopy(3, MarkerAction.Highlight)).toBe(true);
    expect(m.markCopy(3, MarkerAction.Highlight)).toBe(false);
    expect(m.markedCount).toBe(1);
    expect(m.words[6]).toBe(1);
    expect(m.byteLength).toBe(5 * 8);
    expect(read(m, 3, 5 * StorePageSize + 7)).toBe(1);
    expect(read(m, 4, 5 * StorePageSize + 7)).toBe(0);
    expect(m.markCopy(3, MarkerAction.RemoveHighlight)).toBe(true);
    expect(m.markedCount).toBe(0);
  });

  it("marks a copy's particles in a block by its template's store indices, let go once clear", () => {
    const m = create();
    const a = 2 * StorePageSize + 1;
    const b = 3 * StorePageSize + 4;
    expect(m.markIndices(1, SortedArray.ofSortedArray([a, b]), MarkerAction.Select)).toBe(true);
    expect(m.markedCount).toBe(1);
    // the block the template's two pages
    expect(m.byteLength).toBe(5 * 8 + 2 * StorePageSize);
    expect(read(m, 1, a)).toBe(2);
    expect(read(m, 1, b)).toBe(2);
    expect(read(m, 1, a + 1)).toBe(0);
    // the other copy of the template not
    expect(read(m, 0, a)).toBe(0);
    expect(m.get(1, a)).toBe(2);
    // a copy whole over its block
    expect(m.markCopy(1, MarkerAction.Highlight)).toBe(true);
    expect(read(m, 1, a)).toBe(3);
    expect(read(m, 1, a + 1)).toBe(1);
    expect(m.markCopy(1, MarkerAction.RemoveHighlight)).toBe(true);
    expect(read(m, 1, a + 1)).toBe(0);
    // cleared: its block let go (its word its marker again), reused by the template's other copy
    expect(m.markIndices(1, Interval.ofBounds(a, b + 1), MarkerAction.Deselect)).toBe(true);
    expect(m.markedCount).toBe(0);
    expect(m.words[2] & CopyMarkerBlock).toBe(0);
    expect(m.markIndices(0, SortedArray.ofSortedArray([b]), MarkerAction.Highlight)).toBe(true);
    expect(m.byteLength).toBe(5 * 8 + 2 * StorePageSize);
    expect(read(m, 0, b)).toBe(1);
    expect(read(m, 0, a)).toBe(0);
  });

  it('a copy marked whole keeps its marker in a block made for a part', () => {
    const m = create();
    m.markCopy(2, MarkerAction.Select);
    const i = 5 * StorePageSize + 3;
    expect(m.markIndices(2, SortedArray.ofSortedArray([i]), MarkerAction.Highlight)).toBe(true);
    expect(read(m, 2, i)).toBe(3);
    expect(read(m, 2, i + 1)).toBe(2);
    expect(m.markedCount).toBe(1);
  });

  it("takes a template's marks (by store index) into a copy's own", () => {
    const m = create();
    const values = new Uint8Array(8 * StorePageSize);
    const i = 3 * StorePageSize + 2;
    values[i] = 2;
    // the first template's copies: a block each with the mark; the second's (page 5) none to take
    expect(m.include(0, values)).toBe(true);
    expect(m.include(3, values)).toBe(false);
    expect(m.isMarked(0)).toBe(true);
    expect(m.isMarked(3)).toBe(false);
    expect(read(m, 0, i)).toBe(2);
    expect(read(m, 0, i - 1)).toBe(0);
    // or-ed with its own
    m.markIndices(1, SortedArray.ofSortedArray([i]), MarkerAction.Highlight);
    m.include(1, values);
    expect(read(m, 1, i)).toBe(3);
    expect(m.markedCount).toBe(2);
  });

  it('marks every copy, and lays out again cleared', () => {
    const m = create();
    m.markIndices(4, SortedArray.ofSortedArray([5 * StorePageSize]), MarkerAction.Highlight);
    expect(m.markAll(MarkerAction.Select)).toBe(true);
    expect(m.markedCount).toBe(5);
    for (let c = 0; c < 4; ++c) expect(m.words[2 * c]).toBe(2);
    expect(read(m, 4, 5 * StorePageSize)).toBe(3);
    expect(m.markAll(MarkerAction.Clear)).toBe(true);
    expect(m.markedCount).toBe(0);
    m.markCopy(0, MarkerAction.Select);
    m.layout([{ pages: Int32Array.of(0), copyCount: 1 }]);
    expect(m.copyCount).toBe(1);
    expect(m.markedCount).toBe(0);
    expect(m.words[0]).toBe(0);
  });

  it('takes away from a copy marked in part by its marked particles, as from its block whole', () => {
    // a block's marked particles: highlighted and selected, one marked twice (cleared and marked again between)
    const p = 2 * StorePageSize;
    const run = (positions: number[], clearAll: boolean) => {
      const m = create();
      m.markIndices(0, SortedArray.ofSortedArray(positions), MarkerAction.Highlight);
      m.markIndices(0, SortedArray.ofSortedArray([positions[0]]), MarkerAction.RemoveHighlight);
      m.markIndices(0, SortedArray.ofSortedArray([positions[0]]), MarkerAction.Select);
      m.markIndices(0, SortedArray.ofSortedArray([positions[0]]), MarkerAction.Highlight);
      // more than a sixteenth of the block marked: its block whole
      if (clearAll) m.markIndices(0, Interval.ofBounds(p, p + StorePageSize), MarkerAction.Highlight);
      const removed = m.markCopy(0, MarkerAction.RemoveHighlight);
      const left = [...positions, p + 9].map((i) => m.get(0, i));
      const marked = m.markedCount;
      const deselected = m.markCopy(0, MarkerAction.Deselect);
      return { removed, left, marked, deselected, after: m.markedCount };
    };
    const sparse = run([p + 3, p + 17, 3 * StorePageSize + 5], false);
    expect(sparse).toEqual({ removed: true, left: [2, 0, 0, 0], marked: 1, deselected: true, after: 0 });
    const dense = run([p + 3, p + 17, 3 * StorePageSize + 5], true);
    expect(dense).toEqual(sparse);
    // a whole copy highlighted after a part: every particle, then taken away by them all
    const m = create();
    m.markIndices(0, SortedArray.ofSortedArray([p + 3]), MarkerAction.Select);
    expect(m.markCopy(0, MarkerAction.Highlight)).toBe(true);
    expect(m.get(0, p + 9)).toBe(1);
    expect(m.get(0, p + 3)).toBe(3);
    expect(m.markCopy(0, MarkerAction.RemoveHighlight)).toBe(true);
    expect(m.get(0, p + 9)).toBe(0);
    expect(m.get(0, p + 3)).toBe(2);
  });

  it("tells where a copy is marked: whole, its block's positions, or not known once filled", () => {
    const m = create();
    const p = 2 * StorePageSize;
    m.markCopy(0, MarkerAction.Highlight);
    expect(m.getMarkedPositions(0)).toBeUndefined();
    m.markIndices(1, SortedArray.ofSortedArray([p + 3, p + 17]), MarkerAction.Select);
    expect(m.getMarkedPositions(1)).toEqual([3, 17]);
    // past a sixteenth of its block: not known
    m.markIndices(1, Interval.ofBounds(p, p + StorePageSize), MarkerAction.Highlight);
    expect(m.getMarkedPositions(1)).toBeNull();
  });

  it('hands out the bytes changed: all once laid out, then what an action changed', () => {
    const m = create();
    expect(m.takeChanged()).toBeUndefined();
    expect(m.takeChanged()).toEqual([]);
    const p = 2 * StorePageSize;
    // copy 1's block after copy 0's (two pages each): its particles 3 and 2 pages apart, two ranges
    m.markIndices(0, SortedArray.ofSortedArray([p + 3]), MarkerAction.Highlight);
    m.markIndices(1, SortedArray.ofSortedArray([p + 5, p + 9]), MarkerAction.Highlight);
    const b1 = 2 * StorePageSize;
    expect(m.takeChanged()).toEqual([3, 4, b1 + 5, b1 + 10]);
    // a hover's clear: the marked bytes alone
    m.markAll(MarkerAction.RemoveHighlight);
    expect(m.takeChanged()).toEqual([3, 4, b1 + 5, b1 + 10]);
    expect(m.markedCount).toBe(0);
    // a block made of a copy marked whole: its every byte
    m.markCopy(2, MarkerAction.Select);
    m.markIndices(2, SortedArray.ofSortedArray([5 * StorePageSize]), MarkerAction.Highlight);
    const b2 = m.byteLength - 5 * 8 - StorePageSize;
    expect(m.takeChanged()).toEqual([b2, b2 + StorePageSize]);
  });

  it('marks a copy in part whole past its budget: no block, a warning once a layout', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    let budget = 2 * StorePageSize;
    const m = createCopyMarkers(() => budget);
    const layout = () =>
      m.layout([
        { pages: Int32Array.of(2, 3), copyCount: 2 },
        { pages: Int32Array.of(5), copyCount: 3 },
      ]);
    layout();
    const a = 2 * StorePageSize + 1;
    // copy 0's block the budget's whole
    expect(m.markIndices(0, SortedArray.ofSortedArray([a]), MarkerAction.Select)).toBe(true);
    expect(m.byteLength).toBe(5 * 8 + 2 * StorePageSize);
    expect(m.coarseCount).toBe(0);
    // copy 1 and copy 2 (the other template) none: marked whole, their words
    expect(m.markIndices(1, SortedArray.ofSortedArray([a]), MarkerAction.Select)).toBe(true);
    expect(m.markIndices(2, SortedArray.ofSortedArray([5 * StorePageSize + 3]), MarkerAction.Highlight)).toBe(true);
    expect(m.byteLength).toBe(5 * 8 + 2 * StorePageSize);
    expect(m.words[2]).toBe(2);
    expect(m.get(1, a + 1)).toBe(2);
    expect(m.getMarkedPositions(1)).toBeUndefined();
    expect(m.get(2, 5 * StorePageSize + 9)).toBe(1);
    expect(m.markedCount).toBe(3);
    expect(m.coarseCount).toBe(2);
    expect(warn).toHaveBeenCalledTimes(1);
    // indices outside its template's pages: nothing
    expect(m.markIndices(3, SortedArray.ofSortedArray([2 * StorePageSize]), MarkerAction.Select)).toBe(false);
    expect(m.isMarked(3)).toBe(false);
    // a removal on a copy marked whole clears it whole, not counted
    expect(m.markIndices(1, SortedArray.ofSortedArray([a]), MarkerAction.Deselect)).toBe(true);
    expect(m.isMarked(1)).toBe(false);
    expect(m.coarseCount).toBe(2);
    // copy 0's block let go: copy 1 (its template's) takes it, in part again
    m.markIndices(0, SortedArray.ofSortedArray([a]), MarkerAction.Deselect);
    expect(m.markIndices(1, SortedArray.ofSortedArray([a]), MarkerAction.Select)).toBe(true);
    expect(m.get(1, a)).toBe(2);
    expect(m.get(1, a + 1)).toBe(0);
    // a template's marks or-ed into a copy with no room: whole, every value's
    const values = new Uint8Array(6 * StorePageSize);
    values[5 * StorePageSize + 4] = 1;
    values[5 * StorePageSize + 6] = 2;
    expect(m.include(3, values)).toBe(true);
    expect(m.words[6]).toBe(3);
    expect(m.coarseCount).toBe(3);
    expect(warn).toHaveBeenCalledTimes(1);
    // laid out again: the budget read again, the count cleared, warned again once past it
    budget = 4 * StorePageSize;
    layout();
    expect(m.coarseCount).toBe(0);
    m.markIndices(0, SortedArray.ofSortedArray([a]), MarkerAction.Select);
    m.markIndices(1, SortedArray.ofSortedArray([a]), MarkerAction.Select);
    expect(m.coarseCount).toBe(0);
    m.markIndices(2, SortedArray.ofSortedArray([5 * StorePageSize]), MarkerAction.Select);
    expect(m.coarseCount).toBe(1);
    expect(warn).toHaveBeenCalledTimes(2);
  });

  it("takes half a device's copy marker binding as its budget", () => {
    expect(getCopyBlockBudget(undefined)).toBe(DefaultCopyBlockBudget);
    const webgpu = { compute: { maxStorageBufferSize: 2 ** 31 - 4 } } as unknown as GraphicsContext;
    expect(getCopyBlockBudget(webgpu)).toBe(2 ** 30 - 4);
    const webgl = { resources: { limits: { maxTextureSize: 16384 } } } as unknown as GraphicsContext;
    expect(getCopyBlockBudget(webgl)).toBe(128 * 2 ** 20);
  });
});
