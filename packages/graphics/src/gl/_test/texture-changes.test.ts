/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util/value-cell';
import { TextureChanges } from '../renderable/texture-changes.js';
import { createDataArray, type TextureImage } from '../renderable/util.js';

function layer(count: number) {
  return ValueCell.create(createDataArray(count, 1, Uint8Array));
}

describe('texture changes', () => {
  it('hand out the ranges of the logged updates between two versions', () => {
    const cell = layer(10_000);
    const from = cell.ref.version;
    TextureChanges.update(cell, [100, 110]);
    TextureChanges.update(cell, [5000, 5001]);
    const n = cell.ref.value.count;
    expect(TextureChanges.get(cell.ref.value, from, cell.ref.version, n)).toEqual([100, 110, 5000, 5001]);
    // since the first update only
    expect(TextureChanges.get(cell.ref.value, from + 1, cell.ref.version, n)).toEqual([5000, 5001]);
  });

  it('upload all after an update not logged, or a whole change', () => {
    const cell = layer(10_000);
    const from = cell.ref.version;
    TextureChanges.update(cell, [100, 110]);
    ValueCell.update(cell, cell.ref.value);
    const n = 10_000;
    expect(TextureChanges.get(cell.ref.value, from, cell.ref.version, n)).toBeUndefined();
    const whole = cell.ref.version;
    TextureChanges.update(cell, undefined);
    expect(TextureChanges.get(cell.ref.value, whole, cell.ref.version, n)).toBeUndefined();
    // a new image of the cell (resized) has no log
    const next = cell.ref.version;
    ValueCell.update(cell, createDataArray(20_000, 1, Uint8Array));
    TextureChanges.update(cell, [0, 1]);
    expect(TextureChanges.get(cell.ref.value, next, cell.ref.version, n)).toBeUndefined();
  });

  it('upload the changes made without an update with the next update', () => {
    const cell = layer(10_000);
    const n = 10_000;
    const from = cell.ref.version;
    TextureChanges.defer(cell, [100, 110]);
    expect(cell.ref.version).toBe(from);
    TextureChanges.update(cell, [5000, 5001]);
    expect(TextureChanges.get(cell.ref.value, from, cell.ref.version, n)).toEqual([100, 110, 5000, 5001]);
    // pending once only
    const next = cell.ref.version;
    TextureChanges.update(cell, [7000, 7001]);
    expect(TextureChanges.get(cell.ref.value, next, cell.ref.version, n)).toEqual([7000, 7001]);
    // all of it changed meanwhile: the next update uploads all
    const all = cell.ref.version;
    TextureChanges.defer(cell, undefined);
    TextureChanges.defer(cell, [1, 2]);
    TextureChanges.update(cell, [3, 4]);
    expect(TextureChanges.get(cell.ref.value, all, cell.ref.version, n)).toBeUndefined();
  });

  it('merge near ranges, and upload all where scattered or filling much', () => {
    const n = 100_000;
    expect(TextureChanges.mergeRanges([10, 20, 30, 40, 1000, 1010], n)).toEqual([10, 40, 1000, 1010]);
    // every other group of a range: one range
    const alternate: number[] = [];
    for (let i = 0; i < 2000; i += 2) alternate.push(i, i + 1);
    expect(TextureChanges.mergeRanges(alternate, n)).toEqual([0, 1999]);
    // many far apart: all
    const scattered: number[] = [];
    for (let i = 0; i < 20; ++i) scattered.push(i * 1000, i * 1000 + 1);
    expect(TextureChanges.mergeRanges(scattered, n)).toBeUndefined();
    // more than half: all
    expect(TextureChanges.mergeRanges([0, 60_000], n)).toBeUndefined();
  });

  it('clamp to the element count by default', () => {
    const cell = layer(1000);
    const from = cell.ref.version;
    TextureChanges.update(cell, [10, 20, 990, 2000]);
    expect(TextureChanges.get(cell.ref.value, from, cell.ref.version)).toEqual([10, 20, 990, 1000]);
  });
});

describe('texture changes, the rows of an image written in place', () => {
  function image(width: number, height: number) {
    const result: TextureImage<Uint8Array> = { array: new Uint8Array(width * height), width, height };
    return result;
  }

  it('advance its version and hand out the rows written since a version, merged', () => {
    const atlas = image(64, 100);
    TextureChanges.writeRows(atlas, 0, 10);
    expect(atlas.version).toBe(1);
    TextureChanges.writeRows(atlas, 0, 10);
    TextureChanges.writeRows(atlas, 10, 20);
    TextureChanges.writeRows(atlas, 30, 35);
    expect(atlas.version).toBe(4);
    expect(TextureChanges.getImageRows(atlas, 1, 4)).toEqual([0, 20, 30, 35]);
    expect(TextureChanges.getImageRows(atlas, 3, 4)).toEqual([30, 35]);
    // clamped to the image
    TextureChanges.writeRows(atlas, 95, 120);
    expect(TextureChanges.getImageRows(atlas, 4, 5)).toEqual([95, 100]);
  });

  it('hand out none where a version is not logged or much of the image is written', () => {
    const atlas = image(64, 100);
    TextureChanges.writeRows(atlas, 0, 10);
    // a write without rows (the version advanced by hand)
    atlas.version! += 1;
    TextureChanges.writeRows(atlas, 10, 20);
    expect(TextureChanges.getImageRows(atlas, 0, 3)).toBeUndefined();
    expect(TextureChanges.getImageRows(atlas, 2, 3)).toEqual([10, 20]);
    TextureChanges.writeRows(atlas, 20, 80);
    expect(TextureChanges.getImageRows(atlas, 3, 4)).toBeUndefined();
    // nothing in between, or an image never written
    expect(TextureChanges.getImageRows(atlas, 4, 4)).toBeUndefined();
    expect(TextureChanges.getImageRows(image(4, 4), 0, 1)).toBeUndefined();
  });
});
