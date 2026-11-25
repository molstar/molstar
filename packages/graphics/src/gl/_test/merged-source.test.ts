/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util/value-cell';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import { Color } from '@molstar/core/util/color';
import { createObjectColor } from '@molstar/graphics/geo/geometry/color-data';
import { createDataArray } from '../renderable/util.js';
import {
  createMergedValues,
  fillMergeRecord,
  getSourceRecordFillKeys,
  type MergeableValues,
  type MergedInstanceSource,
} from '../renderable/merged.js';
import { SegmentedArray } from '../segmented-array.js';

/** `member`'s record filled with `values` as a source's member is (`getSourceRecordFillKeys`) */
function fill(member: MergeableValues, values: MergeableValues) {
  fillMergeRecord(member, values, getSourceRecordFillKeys('spheres', member, new Set(['tColor', 'tMarker'])));
}

/** `n` members of spheres, their colors and markers an instance each (a page store's members) */
function createMembers(n: number) {
  const p = PD.getDefaultValues(Spheres.Params);
  const members: MergeableValues[] = [];
  for (let i = 0; i < n; ++i) {
    const values = Spheres.Utils.createValuesSimple(
      Spheres.createEmpty(),
      p,
      ColorNames.orange,
      1,
    ) as unknown as MergeableValues;
    ValueCell.update(values.colorType, 'instance');
    ValueCell.update(values.markerType, 'instance');
    members.push(values);
  }
  return members;
}

/** the store's counts as the members' instance counts */
function setCounts(members: readonly MergeableValues[], counts: number[]) {
  counts.forEach((c, i) => ValueCell.updateIfChanged(members[i].instanceCount, c));
}

function createSource(capacity: number): MergedInstanceSource {
  return {
    layers: new Map([
      ['tColor', ValueCell.create(createDataArray(capacity, 4, Uint8Array))],
      ['tMarker', ValueCell.create(createDataArray(capacity, 1, Uint8Array))],
    ]),
    // an empty store's: the merged values alone, no culler
    paged: {
      source: {
        format: undefined,
        records: new Uint32Array(0),
        pageCapacity: 0,
        runs: new Uint32Array(0),
        runCapacity: 0,
        extras: new Uint32Array(0),
        extraCapacity: 0,
        values: [],
        pageCount: 0,
        pool: { centers: new Float64Array(0), count: 0, counts: new Int32Array(0) },
      },
      takeChanges: () => ({ written: [], released: [], runs: [], extras: [], grown: false }),
      reaches: [],
    },
  };
}

/**
 * the merged values' arrays of the members' geometry and layers merged from them: their vertex attributes' parts, and
 * the layers' (spheres': their positions, a texture) as segmented arrays (not the source's)
 */
function geometryArrays(values: Record<string, ValueCell<unknown>>) {
  const out: unknown[] = [];
  for (const [k, c] of Object.entries(values)) {
    const v = c.ref.value as { array?: unknown } | undefined;
    if (k === 'aTransform') continue;
    if (SegmentedArray.is(v)) out.push(v);
    else if (v && SegmentedArray.is(v.array)) out.push(v.array);
  }
  return out;
}

describe('merged values of an instance source', () => {
  it("lays out the segments by the members' counts, its layers the source's cells, no instance data of the members", () => {
    // three targets and the gaps' member last
    const members = createMembers(4);
    setCounts(members, [3, 0, 5, 2]);
    const source = createSource(64);
    const merged = createMergedValues('spheres', members[0], members, true, false, source);
    const { values, segments } = merged;
    expect(values.tColor).toBe(source.layers.get('tColor'));
    expect(values.tMarker).toBe(source.layers.get('tMarker'));
    expect(segments.instanceBases).toEqual([0, 3, 3, 8]);
    expect(segments.instanceCounts).toEqual([3, 0, 5, 2]);
    expect(values.instanceCount.ref.value).toBe(10);
    // a page store's: no array of its instances' ids (its culled poses' words hold them)
    expect((values.aInstance.ref.value as Uint32Array).length).toBe(0);
    expect((values.aTransform.ref.value as SegmentedArray).length).toBe(0);
    expect(values.instancePoses.ref.value).toBeUndefined();
  });

  it('takes an instance count changed alone as the segments again, the geometry and layers kept', () => {
    const members = createMembers(4);
    setCounts(members, [3, 0, 5, 2]);
    const source = createSource(64);
    const merged = createMergedValues('spheres', members[0], members, true, false, source);
    const values = merged.values as unknown as Record<string, ValueCell<unknown>>;
    const geometry = geometryArrays(values);
    expect(geometry.length).toBeGreaterThan(0);
    const cull = merged.cullVersion.ref.value;
    const color = values.tColor.ref.value;
    // a tile added: two members more instances, a new one's first
    setCounts(members, [4, 6, 5, 2]);
    expect(merged.sync()).toBe(true);
    expect(merged.segments.instanceBases).toEqual([0, 4, 10, 15]);
    expect(merged.values.instanceCount.ref.value).toBe(17);
    expect(merged.cullVersion.ref.value).toBe(cull + 1);
    // the geometry's arrays the same objects, the layer the source's as it is
    expect(geometryArrays(values)).toEqual(geometry);
    geometryArrays(values).forEach((a, i) => expect(a).toBe(geometry[i]));
    expect(values.tColor.ref.value).toBe(color);
    // nothing changed: no layout
    expect(merged.sync()).toBe(false);
  });

  it('reads the members its source says it changed alone (takeMembers), their changes counted', () => {
    const members = createMembers(4);
    setCounts(members, [3, 0, 5, 2]);
    let changed: number[] | undefined;
    const source: MergedInstanceSource = {
      ...createSource(64),
      takeMembers: () => {
        const taken = changed;
        changed = [];
        return taken;
      },
    };
    const merged = createMergedValues('spheres', members[0], members, true, false, source);
    // the first sync reads all of them
    expect(merged.sync()).toBe(false);
    const seen = Array.from(merged.memberChanges);
    // two members' counts changed, the source tells of one: the other's cells not read
    setCounts(members, [4, 1, 5, 2]);
    changed = [1];
    expect(merged.sync()).toBe(true);
    expect(merged.segments.instanceCounts).toEqual([3, 1, 5, 2]);
    expect(merged.segments.instanceBases).toEqual([0, 3, 4, 9]);
    expect(Array.from(merged.memberChanges)).toEqual([seen[0], seen[1] + 1, seen[2], seen[3]]);
    // told of the other: read
    changed = [0];
    expect(merged.sync()).toBe(true);
    expect(merged.segments.instanceCounts).toEqual([4, 1, 5, 2]);
    expect(merged.values.instanceCount.ref.value).toBe(12);
    expect(Array.from(merged.memberChanges)).toEqual([seen[0] + 1, seen[1] + 1, seen[2], seen[3]]);
    // none: nothing read, no layout
    setCounts(members, [4, 1, 5, 7]);
    expect(merged.sync()).toBe(false);
    expect(merged.segments.instanceCounts).toEqual([4, 1, 5, 2]);
    // undefined: all read
    changed = undefined;
    expect(merged.sync()).toBe(true);
    expect(merged.segments.instanceCounts).toEqual([4, 1, 5, 7]);
  });

  it("reads the counts and bounds its source's copies hold of the members it counted (counts), not their cells", () => {
    const members = createMembers(3);
    setCounts(members, [3, 0, 2]);
    const instanceCounts = new Int32Array([3, 0, 2]);
    const spheres = new Float64Array(12);
    members.forEach((m, i) =>
      spheres.set([...m.boundingSphere.ref.value.center, m.boundingSphere.ref.value.radius], i * 4),
    );
    let counted: number[] = [];
    const source: MergedInstanceSource = {
      ...createSource(64),
      takeMembers: () => [],
      counts: {
        instanceCounts,
        spheres,
        take: () => {
          const taken = counted;
          counted = [];
          return taken;
        },
      },
    };
    const merged = createMergedValues('spheres', members[0], members, true, false, source);
    expect(merged.sync()).toBe(false);
    const seen = Array.from(merged.memberChanges);
    const cull = merged.cullVersion.ref.value;
    // the copies of member 1 changed (its cells left as they were: not read)
    instanceCounts[1] = 4;
    spheres.set([100, 0, 0, 5], 4);
    counted = [1];
    expect(merged.sync()).toBe(true);
    expect(merged.segments.instanceCounts).toEqual([3, 4, 2]);
    expect(merged.values.instanceCount.ref.value).toBe(9);
    const bounds = merged.values.boundingSphere.ref.value;
    expect(bounds.center[0] + bounds.radius).toBeGreaterThanOrEqual(105 - 1e-6);
    expect(merged.cullVersion.ref.value).toBe(cull + 1);
    // their cells not read as changed: what the culler reads of them (their rows) not changed
    expect(Array.from(merged.memberChanges)).toEqual(seen);
    // its bounds alone: no layout, the bounds and the culling's version
    spheres.set([200, 0, 0, 5], 4);
    counted = [1];
    expect(merged.sync()).toBe(false);
    expect(
      merged.values.boundingSphere.ref.value.center[0] + merged.values.boundingSphere.ref.value.radius,
    ).toBeGreaterThanOrEqual(205 - 1e-6);
    expect(merged.cullVersion.ref.value).toBe(cull + 2);
    // none counted: nothing read
    spheres.set([300, 0, 0, 5], 4);
    expect(merged.sync()).toBe(false);
    expect(merged.cullVersion.ref.value).toBe(cull + 2);
  });

  it("bounds a source's object by its members' spheres, one sphere its own (not its box's, √3 as large)", () => {
    const members = createMembers(3);
    setCounts(members, [3, 0, 2]);
    // the members' bounds none but the last's (a page store's object: the object's sphere its last member's)
    const spheres = new Float64Array(12).fill(Number.NaN);
    spheres.set([10, 20, 30, 50], 8);
    let counted: number[] = [];
    const source: MergedInstanceSource = {
      ...createSource(64),
      takeMembers: () => [],
      counts: {
        instanceCounts: new Int32Array([3, 0, 2]),
        spheres,
        take: () => {
          const taken = counted;
          counted = [];
          return taken;
        },
      },
    };
    const merged = createMergedValues('spheres', members[0], members, true, false, source);
    merged.sync();
    // the copies of all three read
    counted = [0, 1, 2];
    merged.sync();
    const bounds = merged.values.boundingSphere.ref.value;
    expect(Array.from(bounds.center)).toEqual([10, 20, 30]);
    expect(bounds.radius).toBeCloseTo(50, 6);
    // a count changed alone (its sphere NaN as before: NaN as itself): the bounds not unioned again
    const version = merged.values.boundingSphere.ref.version;
    source.counts!.instanceCounts[0] = 4;
    counted = [0];
    merged.sync();
    expect(merged.segments.instanceCounts[0]).toBe(4);
    expect(merged.values.boundingSphere.ref.version).toBe(version);
  });

  it("grows a source's arrays where a member's new geometry outgrows their room: extended longer, the parts kept", () => {
    const members = createMembers(4);
    setCounts(members, [1, 1, 1, 0]);
    const merged = createMergedValues('spheres', members[0], members, true, false, createSource(64));
    const values = merged.values as unknown as Record<string, ValueCell<unknown>>;
    const before = geometryArrays(values) as SegmentedArray[];
    // a member of another geometry, more vertices than the arrays have room for (at least 16k)
    const count = 20_000;
    const builder = SpheresBuilder.create(count);
    for (let i = 0; i < count; ++i) builder.add(i, 0, 0, i);
    const p = PD.getDefaultValues(Spheres.Params);
    const big = Spheres.Utils.createValuesSimple(
      builder.getSpheres(),
      p,
      ColorNames.orange,
      1,
    ) as unknown as MergeableValues;
    ValueCell.update(big.colorType, 'instance');
    ValueCell.update(big.markerType, 'instance');
    fill(members[1], big);
    merged.replaceGeometry(1);
    merged.sync();
    const after = geometryArrays(values) as SegmentedArray[];
    expect(after.length).toBe(before.length);
    // the layer laid out by allocation (the spheres' positions): an extension of the array before (an uploader writes
    // the added part alone), its parts the same objects where they were, longer (room again); the others (the first
    // member's) as they were
    const positions = (values.tPositionGroup.ref.value as { array: SegmentedArray }).array;
    const i = after.indexOf(positions);
    expect(i).toBeGreaterThanOrEqual(0);
    expect(positions.from).toBe(before[i]);
    before[i].parts.forEach((part, k) => expect(positions.parts[k]).toBe(part));
    expect(positions.parts.length).toBe(before[i].parts.length + 1);
    expect(positions.length).toBeGreaterThan(before[i].length);
    expect(positions.parts.at(-1)!.count).toBe(count * 4);
    after.forEach((a, k) => {
      if (k !== i) expect(a).toBe(before[k]);
    });
  });

  it("lays out a member's new geometry after the others', theirs kept; the draws' max of the live ones", () => {
    const members = createMembers(4);
    setCounts(members, [1, 1, 1, 0]);
    const merged = createMergedValues('spheres', members[0], members, true, false, createSource(64));
    const p = PD.getDefaultValues(Spheres.Params);
    const spheres = (count: number) => {
      const builder = SpheresBuilder.create(count);
      for (let i = 0; i < count; ++i) builder.add(i, 0, 0, i);
      const v = Spheres.Utils.createValuesSimple(
        builder.getSpheres(),
        p,
        ColorNames.orange,
        1,
      ) as unknown as MergeableValues;
      ValueCell.update(v.colorType, 'instance');
      ValueCell.update(v.markerType, 'instance');
      return v;
    };
    const fillWith = (i: number, count: number) => {
      const v = spheres(count);
      fill(members[i], v);
      merged.replaceGeometry(i);
      merged.sync();
      return v;
    };
    const { segments, values } = merged;
    const a = fillWith(1, 100);
    const end = segments.vertexBases[1] + a.uVertexCount.ref.value;
    const drawEnd = segments.drawBases[1] + a.drawCount.ref.value;
    expect(values.drawCountMax.ref.value).toBe(a.drawCount.ref.value);
    const kept = Array.from(segments.vertexBases);
    const b = fillWith(2, 50);
    // appended: after the first's, the others where they were
    expect(segments.vertexBases[2]).toBe(end);
    expect(segments.drawBases[2]).toBe(drawEnd);
    expect(segments.elementOffsets[2]).toBe(drawEnd * 4);
    for (const i of [0, 1, 3]) expect(segments.vertexBases[i]).toBe(kept[i]);
    expect(values.drawCountMax.ref.value).toBe(a.drawCount.ref.value);
    // the first's smaller again: its unused allocation more than the used, laid out whole; the max the live draws'
    fillWith(1, 10);
    expect(values.drawCountMax.ref.value).toBe(b.drawCount.ref.value);
    expect(segments.vertexBases[2]).toBe(segments.vertexBases[1] + members[1].uVertexCount.ref.value);
  });

  it("keeps a member's color by member (colorType 'object') at its index: filled again, recolored, after the layout", () => {
    // a member its color (as a page store's object colors its members: `setMemberColor`)
    const setColor = (values: MergeableValues, color: Color) => {
      const c = createObjectColor(color) as unknown as Record<string, ValueCell<unknown>>;
      const cells = values as unknown as Record<string, ValueCell<unknown>>;
      for (const k in c) if (cells[k]) ValueCell.update(cells[k], c[k].ref.value);
    };
    const rgb = (color: Color) => [(color >> 16) & 255, (color >> 8) & 255, color & 255];
    const members = createMembers(4);
    const colors = [Color(0x112233), Color(0x445566), Color(0x778899), Color(0xaabbcc)];
    members.forEach((m, i) => setColor(m, colors[i]));
    setCounts(members, [1, 1, 1, 0]);
    const merged = createMergedValues('spheres', members[0], members, true, false, createSource(64));
    // member i's texel at i: what the shader reads (`ColorTypeIndex` 'object', the segment's index)
    const texel = (i: number) => {
      const array = (merged.values.tColor.ref.value as { array: SegmentedArray }).array;
      return Array.from(SegmentedArray.read(array, i * 4, i * 4 + 3));
    };
    for (let i = 0; i < 4; ++i) expect(texel(i)).toEqual(rgb(colors[i]));

    // a member filled with another target's geometry and colored (allocated again, at the arrays' end): its texel its
    // index's still
    const builder = SpheresBuilder.create(8);
    for (let i = 0; i < 8; ++i) builder.add(i, 0, 0, i);
    const other = Spheres.Utils.createValuesSimple(
      builder.getSpheres(),
      PD.getDefaultValues(Spheres.Params),
      ColorNames.orange,
      1,
    ) as unknown as MergeableValues;
    fill(members[1], other);
    setColor(members[1], Color(0x0000ff));
    merged.replaceGeometry(1);
    merged.sync();
    expect(texel(1)).toEqual(rgb(Color(0x0000ff)));
    // a member recolored alone
    setColor(members[2], Color(0x00ff00));
    merged.sync();
    expect(texel(2)).toEqual(rgb(Color(0x00ff00)));
    for (const i of [0, 3]) expect(texel(i)).toEqual(rgb(colors[i]));
  });

  it('without a source an instance count changed is a new layout, the geometry merged again', () => {
    const members = createMembers(2);
    setCounts(members, [3, 2]);
    const merged = createMergedValues('spheres', members[0], members, true, false);
    const values = merged.values as unknown as Record<string, ValueCell<unknown>>;
    const geometry = geometryArrays(values);
    expect(geometry.length).toBeGreaterThan(0);
    setCounts(members, [4, 2]);
    expect(merged.sync()).toBe(true);
    expect(geometryArrays(values).some((a, i) => a !== geometry[i])).toBe(true);
  });
});
