/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util/value-cell';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { Points } from '@molstar/graphics/geo/geometry/points/points';
import { createRenderObject } from '../render-object.js';
import { RenderableState, type GraphicsRenderable } from '../renderable.js';
import { createDrawTable, DrawFact, DrawRow, DrawState, getDrawFacts } from '../draw-table.js';

/** A render object standing in for its renderable: the table reads only the values and the state */
function meshObject(props: Partial<PD.Values<Mesh.Params>> = {}) {
  const p = { ...PD.getDefaultValues(Mesh.Params), ...props };
  const values = Mesh.Utils.createValuesSimple(Mesh.createEmpty(), p, ColorNames.orange, 1);
  ValueCell.update(values.drawCount, 3);
  return createRenderObject('mesh', values, Mesh.Utils.createRenderableState(p), -1) as unknown as GraphicsRenderable;
}

function pointsObject() {
  const p = { ...PD.getDefaultValues(Points.Params), pointStyle: 'fuzzy' as const };
  const values = Points.Utils.createValuesSimple(Points.createEmpty(), p, ColorNames.orange, 1);
  return createRenderObject(
    'points',
    values,
    Points.Utils.createRenderableState(p),
    -1,
  ) as unknown as GraphicsRenderable;
}

function row(table: ReturnType<typeof createDrawTable>, r: GraphicsRenderable) {
  const o = table.getSlot(r) * DrawRow.Stride;
  return {
    sphere: Array.from(table.floats.subarray(o, o + 4)),
    lod: [table.floats[o + DrawRow.LodMin], table.floats[o + DrawRow.LodMax]],
    facts: table.bits[o + DrawRow.Facts],
    state: table.bits[o + DrawRow.State],
  };
}

describe('draw table', () => {
  it('writes the rows of added objects, in the group order', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const b = pointsObject();
    ValueCell.update(a.values.boundingSphere, Sphere3D.create(Vec3.create(1, 2, 3), 4));
    ValueCell.update(a.values.uLod, Vec4.create(10, 20, 0, 0));
    for (const r of [a, b]) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();

    expect(table.count).toBe(2);
    expect(Array.from(table.order.subarray(0, 2))).toEqual([table.getSlot(a), table.getSlot(b)]);
    expect(row(table, a).sphere).toEqual([1, 2, 3, 4]);
    expect(row(table, a).lod).toEqual([10, 20]);
    expect(row(table, a).facts).toBe(getDrawFacts(a.values));
    expect(row(table, a).facts & DrawFact.DrawCount).toBeTruthy();
    expect(row(table, a).facts & DrawFact.AlphaOne).toBeTruthy();
    expect(row(table, b).facts & DrawFact.FuzzyPoints).toBeTruthy();
    expect(row(table, b).facts & DrawFact.DrawCount).toBe(0);
    // sides: a mesh's given (`uDoubleSided`), points have none
    const sides = DrawFact.SingleSided | DrawFact.DoubleSided;
    expect(row(table, a).facts & sides).toBe(
      a.values.uDoubleSided.ref.value ? DrawFact.DoubleSided : DrawFact.SingleSided,
    );
    expect(row(table, b).facts & sides).toBe(0);
    expect(row(table, a).state).toBe(DrawState.Visible | DrawState.Pickable | DrawState.Opaque | DrawState.WriteDepth);
  });

  it('follows the values the facts read, and only those', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    renderables.push(a);
    table.add(a);
    table.refresh();
    const slot = table.getSlot(a);
    const version = table.version;
    const rowVersion = table.rowVersions[slot];

    ValueCell.update(a.values.uAlpha, 0.5);
    ValueCell.update(a.values.markerAverage, 1);
    ValueCell.update(a.values.dTransparentBackfaces, 'opaque');
    expect(table.rowVersions[slot]).toBe(rowVersion);
    table.refresh();
    expect(table.rowVersions[slot]).toBe(rowVersion + 1);
    expect(table.version).toBeGreaterThan(version);
    const { facts } = row(table, a);
    expect(facts & DrawFact.AlphaOne).toBe(0);
    expect(facts & DrawFact.MarkerAny).toBeTruthy();
    expect(facts & DrawFact.MarkerFull).toBeTruthy();
    expect(facts & DrawFact.BackfacesOpaque).toBeTruthy();
    expect(facts).toBe(getDrawFacts(a.values));

    // a cell the facts don't read
    ValueCell.update(a.values.uMetalness, 0.3);
    table.refresh();
    expect(table.rowVersions[slot]).toBe(rowVersion + 1);
  });

  it('advances the values versions with any value, at once; volatile cells on refresh', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const b = meshObject();
    const volatile = { ...ValueCell.create(0), volatile: true };
    (b.values as unknown as { [k: string]: unknown }).uVolatileTest = volatile;
    for (const r of [a, b]) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();
    const slotA = table.getSlot(a);
    const slotB = table.getSlot(b);
    const versionA = table.valuesVersions[slotA];
    const versionB = table.valuesVersions[slotB];
    const epoch = table.valuesEpoch;
    const rowVersion = table.rowVersions[slotA];

    // a cell the facts don't read: the values version, not the row
    ValueCell.update(a.values.uMetalness, 0.3);
    expect(table.valuesVersions[slotA]).toBe(versionA + 1);
    expect(table.valuesVersions[slotB]).toBe(versionB);
    expect(table.valuesEpoch).toBe(epoch + 1);
    table.refresh();
    expect(table.rowVersions[slotA]).toBe(rowVersion);

    // changed without notifying: compared on refresh
    volatile.ref = { ...volatile.ref, version: volatile.ref.version + 1 };
    expect(table.valuesVersions[slotB]).toBe(versionB);
    table.refresh();
    expect(table.valuesVersions[slotB]).toBe(versionB + 1);
    table.refresh();
    expect(table.valuesVersions[slotB]).toBe(versionB + 1);

    // removed: no longer followed
    renderables.splice(0, 1);
    table.remove(a);
    const removed = table.valuesEpoch;
    ValueCell.update(a.values.uMetalness, 0.6);
    expect(table.valuesEpoch).toBe(removed);
  });

  it('follows the state through its setters', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    renderables.push(a);
    table.add(a);
    table.refresh();
    const slot = table.getSlot(a);
    const rowVersion = table.rowVersions[slot];

    // unchanged writes don't queue
    a.state.visible = true;
    table.refresh();
    expect(table.rowVersions[slot]).toBe(rowVersion);

    a.state.visible = false;
    Object.assign(a.state, { colorOnly: true });
    table.refresh();
    expect(table.rowVersions[slot]).toBe(rowVersion + 1);
    expect(row(table, a).state & DrawState.Visible).toBe(0);
    expect(row(table, a).state & DrawState.ColorOnly).toBeTruthy();
  });

  it('compares plain states when refreshing', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const state: RenderableState = {
      disposed: false,
      visible: true,
      pickable: true,
      colorOnly: false,
      opaque: true,
      writeDepth: true,
    };
    const plain = { values: a.values, state } as unknown as GraphicsRenderable;
    expect(RenderableState.isObservable(state)).toBe(false);
    renderables.push(plain);
    table.add(plain);
    table.refresh();
    state.pickable = false;
    table.refresh();
    expect(row(table, plain).state & DrawState.Pickable).toBe(0);
  });

  it('frees the slots of removed objects for reuse, and stops following them', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const b = meshObject();
    for (const r of [a, b]) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();
    const slotA = table.getSlot(a);

    renderables.splice(0, 1);
    table.remove(a);
    table.refresh();
    expect(table.getSlot(a)).toBe(-1);
    expect(table.objects[slotA]).toBeUndefined();
    expect(table.bits[slotA * DrawRow.Stride + DrawRow.State]).toBe(DrawState.Disposed);
    expect(Array.from(table.order.subarray(0, table.count))).toEqual([table.getSlot(b)]);

    const version = table.version;
    ValueCell.update(a.values.uAlpha, 0.2);
    a.state.visible = false;
    table.refresh();
    expect(table.version).toBe(version);

    const c = meshObject();
    renderables.push(c);
    table.add(c);
    table.refresh();
    expect(table.getSlot(c)).toBe(slotA);
    expect(table.capacity).toBeGreaterThanOrEqual(2);
  });

  it('follows a reordered group when told', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const b = meshObject();
    for (const r of [a, b]) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();
    renderables.reverse();
    table.invalidateOrder();
    const version = table.version;
    table.refresh();
    expect(Array.from(table.order.subarray(0, 2))).toEqual([table.getSlot(b), table.getSlot(a)]);
    expect(table.version).toBeGreaterThan(version);
  });

  it('writes what objects add to the scene averages with their rows', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const b = pointsObject();
    for (const r of [a, b]) {
      renderables.push(r);
      table.add(r);
    }
    ValueCell.update(a.values.uAlpha, 0.8);
    ValueCell.update(a.values.markerAverage, 0.25);
    ValueCell.update(a.values.emissiveAverage, 0.1);
    ValueCell.update(a.values.uEmissive, 0.2);
    ValueCell.update(a.values.transparencyAverage, 0.5);
    table.refresh();
    const c = table.contributions;
    const sa = table.getSlot(a);
    const sb = table.getSlot(b);
    expect(c.marker[sa]).toBe(0.25);
    expect(c.emissive[sa]).toBeCloseTo(0.3);
    expect(c.opacity[sa]).toBeCloseTo(0.4);
    expect(c.transparencyMin[sa]).toBeCloseTo(0.2);
    expect(c.opaque[sa]).toBe(0);
    // fuzzy points: half transparent
    expect(c.opacity[sb]).toBeCloseTo(0.5);
    expect(c.transparencyMin[sb]).toBe(0.5);
    expect(c.opaque[sb]).toBe(1);

    expect(row(table, a).facts & DrawFact.Animated).toBe(0);
    ValueCell.update(a.values.uWiggleAmplitude, 1);
    ValueCell.update(a.values.uWiggleSpeed, 1);
    ValueCell.update(a.values.uWiggleFrequency, 1);
    table.refresh();
    expect(row(table, a).facts & DrawFact.Animated).toBeTruthy();
  });

  it('lists the slots whose values changed, each once', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const objects = [meshObject(), meshObject(), meshObject()];
    for (const r of objects) {
      renderables.push(r);
      table.add(r);
    }
    expect(table.takeChanged().sort()).toEqual(objects.map((r) => table.getSlot(r)).sort());
    expect(table.takeChanged()).toEqual([]);
    ValueCell.update(objects[1].values.markerAverage, 1);
    ValueCell.update(objects[1].values.uMetalness, 0.3);
    expect(table.takeChanged()).toEqual([table.getSlot(objects[1])]);
  });

  it('advances the visible version with what is visible only', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const a = meshObject();
    const b = meshObject();
    for (const r of [a, b]) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();
    let version = table.visibleVersion;
    // values and other state: the same visibility
    ValueCell.update(a.values.uAlpha, 0.5);
    a.state.pickable = false;
    table.refresh();
    expect(table.visibleVersion).toBe(version);
    a.state.visible = false;
    table.refresh();
    expect(table.visibleVersion).toBeGreaterThan(version);
    version = table.visibleVersion;
    // removed: no longer visible
    renderables.splice(renderables.indexOf(b), 1);
    table.remove(b);
    table.refresh();
    expect(table.visibleVersion).toBeGreaterThan(version);
  });

  it('grows, keeping the rows', () => {
    const renderables: GraphicsRenderable[] = [];
    const table = createDrawTable(renderables);
    const objects = Array.from({ length: 100 }, () => meshObject());
    objects.forEach((r, i) => ValueCell.update(r.values.boundingSphere, Sphere3D.create(Vec3.create(i, 0, 0), 1)));
    for (const r of objects.slice(0, 10)) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();
    for (const r of objects.slice(10)) {
      renderables.push(r);
      table.add(r);
    }
    table.refresh();
    expect(table.capacity).toBeGreaterThanOrEqual(100);
    for (let i = 0; i < 100; ++i) expect(row(table, objects[i]).sphere[0]).toBe(i);
  });
});
