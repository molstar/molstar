/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util/value-cell';
import { createCullView, LodView } from '../renderable.js';

describe('cull view', () => {
  it('skips a cull of the same view and data only: view version, cell versions, no occlusion test', () => {
    const lod = { ...LodView(), version: 3 };
    const cell = ValueCell.create(1);
    const cells = [cell, undefined];

    const view = createCullView();
    expect(view.isSame(lod, null, cells)).toBe(false);
    view.set(lod, cells);
    expect(view.isSame(lod, null, cells)).toBe(true);

    // an occlusion test reads levels that change between frames
    expect(view.isSame(lod, () => false, cells)).toBe(false);
    // another view (the renderer's version moved), or none known
    expect(view.isSame({ ...lod, version: 4 }, null, cells)).toBe(false);
    view.set(LodView(), cells);
    expect(view.isSame(LodView(), null, cells)).toBe(false);
    view.set(lod, cells);
    // the data: a cell updated, even to the same value
    ValueCell.update(cell, 1);
    expect(view.isSame(lod, null, cells)).toBe(false);
    view.set(lod, cells);
    expect(view.isSame(lod, null, cells)).toBe(true);
    // the lists changed otherwise (`uncull`, `cullSimple`)
    view.invalidate();
    expect(view.isSame(lod, null, cells)).toBe(false);
  });
});
