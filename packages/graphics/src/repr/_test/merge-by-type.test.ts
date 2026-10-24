/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { createRenderObject, type GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { type MergedEntry, mergeRenderObjectsByType } from '../merged.js';

function spheres(n: number) {
  const centers = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) centers[i * 3] = i;
  const s = Spheres.create(centers, new Float32Array(n), n);
  const props = PD.getDefaultValues(Spheres.Params);
  const values = Spheres.Utils.createValuesSimple(s, props, ColorNames.orange, 1);
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

describe('mergeRenderObjectsByType', () => {
  it('skips the merge check for the same members, unchanged since it passed (perf-audit-2 07/P17)', () => {
    const entries = [spheres(3), spheres(2), spheres(4)].map((renderObject, visual) => ({ visual, renderObject }));
    const out: GraphicsRenderObject[] = [];
    let merged: Map<string, MergedEntry<number>> = mergeRenderObjectsByType(entries, true, new Map(), out);
    expect(out.length).toBe(1);
    const object = out[0];

    // the check sorts each member's define keys (getMergeKey)
    const sort = vi.spyOn(Array.prototype, 'sort');
    merged = mergeRenderObjectsByType(entries, true, merged, out);
    expect(sort).not.toHaveBeenCalled();
    sort.mockRestore();
    expect(out).toEqual([object]);

    // a member whose defines changed is checked again: here it no longer merges
    ValueCell.update(entries[1].renderObject.values.dXrayShaded, 'on');
    merged = mergeRenderObjectsByType(entries, true, merged, out);
    expect(out.length).toBe(3);
    expect(merged.size).toBe(0);
  });
});
