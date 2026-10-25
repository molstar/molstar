/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Interval } from '@molstar/core/data/int/interval';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { Loci } from '@molstar/model/model/loci';
import { createMarkers } from '@molstar/graphics/geo/geometry/marker-data';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { TextureChanges } from '@molstar/graphics/gl/renderable/texture-changes';
import { Visual } from '../visual.js';

const Count = 1000;

/** The values `Visual.mark` reads: markers per group of one instance */
function renderObject() {
  const values = {
    ...createMarkers(Count, 'groupInstance'),
    uGroupCount: ValueCell.create(Count),
    instanceCount: ValueCell.create(1),
    uInstanceGranularity: ValueCell.create(false),
  };
  return { values } as unknown as GraphicsRenderObject;
}

/** A loci of the groups `[start, end)` */
const range = (start: number, end: number) => ({ kind: 'test-range', start, end }) as unknown as Loci;
const lociApply: Visual.LociApply = (loci, apply) => {
  const { start, end } = loci as unknown as { start: number; end: number };
  return apply(Interval.ofBounds(start, end));
};

describe('marking uploads', () => {
  it('a mark cleared while the uniform stands for the layer is uploaded with the next one', () => {
    const ro = renderObject();
    const { tMarker } = ro.values as unknown as { tMarker: ValueCell<{ array: Uint8Array }> };
    const image = tMarker.ref.value as Parameters<typeof TextureChanges.get>[0];

    // a residue highlighted: the layer is read (status -1), its range uploaded
    Visual.mark(ro, range(10, 20), MarkerAction.Highlight, lociApply);
    const uploaded = tMarker.ref.version;
    // cleared: all clear, the uniform stands for the layer (status 0) and the cell is not updated
    Visual.mark(ro, range(10, 20), MarkerAction.RemoveHighlight, lociApply);
    expect(tMarker.ref.version).toBe(uploaded);
    // another residue highlighted: what is uploaded since the first upload clears the first residue too
    Visual.mark(ro, range(500, 510), MarkerAction.Highlight, lociApply);
    const ranges = TextureChanges.get(image, uploaded, tMarker.ref.version, Count);
    expect(ranges).toBeDefined();
    const covers = (i: number) => ranges!.some((_, k) => k % 2 === 0 && ranges![k] <= i && i < ranges![k + 1]);
    expect(covers(10) && covers(19)).toBe(true);
    expect(covers(500) && covers(509)).toBe(true);
    expect(tMarker.ref.value.array[10]).toBe(0);
  });
});
