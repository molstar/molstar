/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Vec3 } from '@molstar/core/math/linear-algebra';
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

/**
 * Runs `merge` and returns whether it checked the members again: the check (`canMergeRenderObjects`) calls
 * `Array.prototype.includes` (`isMergeableType`), the skip does not. `getMergeKey` caches its sorted keys, so `sort`
 * tells nothing.
 */
function checks(merge: () => void) {
  const includes = vi.spyOn(Array.prototype, 'includes');
  try {
    merge();
    return includes.mock.calls.length > 0;
  } finally {
    includes.mockRestore();
  }
}

describe('mergeRenderObjectsByType', () => {
  it('skips the merge check for the same members, unchanged since it passed (perf-audit-2 07/P17)', () => {
    const entries = [spheres(3), spheres(2), spheres(4)].map((renderObject, visual) => ({ visual, renderObject }));
    const out: GraphicsRenderObject[] = [];
    let merged: Map<string, MergedEntry<number>> = mergeRenderObjectsByType(entries, true, new Map(), out);
    expect(out.length).toBe(1);
    const object = out[0];

    expect(
      checks(() => {
        merged = mergeRenderObjectsByType(entries, true, merged, out);
      }),
    ).toBe(false);
    expect(out).toEqual([object]);

    // a member whose defines changed is checked again: here it no longer merges
    ValueCell.update(entries[1].renderObject.values.dXrayShaded, 'on');
    merged = mergeRenderObjectsByType(entries, true, merged, out);
    expect(out.length).toBe(3);
    expect(merged.size).toBe(0);
  });

  it('does not check again after a visual update that changes nothing it reads (perf-audit-3 07/P17)', () => {
    const entries = [spheres(3), spheres(2), spheres(4)].map((renderObject, visual) => ({ visual, renderObject }));
    const out: GraphicsRenderObject[] = [];
    let merged: Map<string, MergedEntry<number>> = mergeRenderObjectsByType(entries, true, new Map(), out);
    const object = out[0];
    const props = PD.getDefaultValues(Spheres.Params);
    // what every update of a visual runs (a recolor, a new alpha): `uLod`, `lodLevels` set again
    const update = (p = props) => {
      for (const { renderObject } of entries) {
        Spheres.Utils.updateValues(renderObject.values as any, p);
        ValueCell.update(renderObject.values.lodLevels, renderObject.values.lodLevels.ref.value);
      }
    };

    for (let i = 0; i < 3; i++) {
      update();
      expect(
        checks(() => {
          merged = mergeRenderObjectsByType(entries, true, merged, out);
        }),
      ).toBe(false);
    }
    update({ ...props, alpha: 0.5 });
    expect(
      checks(() => {
        merged = mergeRenderObjectsByType(entries, true, merged, out);
      }),
    ).toBe(false);
    expect(out).toEqual([object]);

    // a new LOD is read by the check: checked again, and still mergeable
    update({ ...props, lod: Vec3.create(1, 100, 0) });
    expect(
      checks(() => {
        merged = mergeRenderObjectsByType(entries, true, merged, out);
      }),
    ).toBe(true);
    expect(out).toEqual([object]);
  });

  it('checks again when any `*Type` cell the merge check reads is set (perf-audit-3 07/P17)', () => {
    const entries = [spheres(3), spheres(2)].map((renderObject, visual) => ({ visual, renderObject }));
    const values = entries[0].renderObject.values as unknown as Record<string, ValueCell<any>>;
    // every `*Type` cell of a member other than the defines (watched as all `d*` are): the ones `canMergeValues`
    // (gl/renderable/merged.ts) reads, and those it does not (`markerType` only shapes the merged layout, in `sync`).
    // A new `*Type` cell fails here until it is put in one of the lists (and in `MergeTypeKeys` if read).
    const read = [
      'colorType',
      'sizeType',
      'overpaintType',
      'transparencyType',
      'emissiveType',
      'substanceType',
      'clippingType',
      'wiggleType',
    ];
    const notRead = ['markerType', 'uColorIndexType', 'uSizeIndexType'];
    expect(
      Object.keys(values)
        .filter((k) => k.endsWith('Type') && !/^d[A-Z]/.test(k))
        .sort(),
    ).toEqual([...read, ...notRead].sort());

    const out: GraphicsRenderObject[] = [];
    let merged: Map<string, MergedEntry<number>> = mergeRenderObjectsByType(entries, true, new Map(), out);
    for (const k of [...read, ...notRead]) {
      ValueCell.update(values[k], values[k].ref.value);
      expect(
        checks(() => {
          merged = mergeRenderObjectsByType(entries, true, merged, out);
        }),
        k,
      ).toBe(read.includes(k));
    }
  });
});
