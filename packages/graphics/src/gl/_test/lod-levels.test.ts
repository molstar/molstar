/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { getLodPixelLevels, getLodPixelLevelsValue } from '@molstar/graphics/geo/geometry/lod-levels';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { ParticleTargetRepresentationParams } from '@molstar/graphics/repr/particles/representation/target/common';
import { createSegmentedMdbList, LodView, type LodLevelsValue } from '../renderable.js';

const Fov = Math.PI / 4;

function perspective() {
  const t = Math.tan(Fov / 2);
  return Mat4.perspective(Mat4(), (-t * 16) / 9, (t * 16) / 9, t, -t, 1, 1000);
}

/** levels by pixel size: 0 to 1, 1 to 4, 4 and up Å a pixel covers */
function pixelLevels(): LodLevelsValue {
  return Object.assign(
    [
      [0, 1, 600, 1],
      [1, 4, 60, 2],
      [4, Infinity, 6, 4],
    ] as LodLevelsValue,
    { pixelSize: 3 },
  );
}

describe('levels of detail by size on screen', () => {
  it('a view of the Å a CSS pixel covers from the projection', () => {
    const view = LodView.set(LodView(), perspective(), 720, false);
    expect(view.depthScale).toBeCloseTo((2 * Math.tan(Fov / 2)) / 720, 12);
    expect(view.pixelSize).toBe(0);

    // a view 300 Å high
    const ortho = LodView.set(LodView(), Mat4.ortho(Mat4(), -200, 200, 150, -150, 1, 1000), 600, true);
    expect(ortho.depthScale).toBe(0);
    expect(ortho.pixelSize).toBeCloseTo(0.5, 12);
  });

  it('ranges in distances: as they are, divided by the depth scale, one level orthographic', () => {
    const range: number[] = [];
    const distances = [[100, 2000, 600, 1]] as LodLevelsValue;
    expect(LodView.getRange(distances, 0, { depthScale: 0.01, pixelSize: 0 }, range)).toEqual([100, 2000]);

    const levels = pixelLevels();
    const view = { depthScale: 0.002, pixelSize: 0 };
    expect(LodView.getRange(levels, 0, view, range)).toEqual([0, 500]);
    expect(LodView.getRange(levels, 1, view, range)).toEqual([500, 2000]);
    // the open end of the last level
    expect(LodView.getRange(levels, 2, view, range)).toEqual([2000, 1e9]);

    // orthographic at 2 Å a pixel: all of level 1, none of the others
    const ortho = { depthScale: 0, pixelSize: 2 };
    expect(LodView.getRange(levels, 0, ortho, range)).toEqual([-2, -1]);
    expect(LodView.getRange(levels, 1, ortho, range)).toEqual([0, 1e9]);
    expect(LodView.getRange(levels, 2, ortho, range)).toEqual([-2, -1]);
  });

  it('level uniforms converted again where the view changes', () => {
    const mdb = createSegmentedMdbList();
    const cell = ValueCell.create(pixelLevels());
    const uLod = () => mdb.list.map((l) => Array.from(l.uniforms[0][1].ref.value as Vec4));

    mdb.prepare(cell, 1, false, { depthScale: 0.002, pixelSize: 0 });
    expect(uLod()).toEqual([
      [0, 500, 0, 1],
      [500, 2000, 0, 2],
      [2000, 1e9, 0, 4],
    ]);
    mdb.prepare(cell, 1, false, { depthScale: 0.004, pixelSize: 0 });
    expect(uLod()[1]).toEqual([250, 1000, 0, 2]);
  });

  it('a level of other uniforms (levels off or on) another object: a draw block is made once a level (WebGPU)', () => {
    const mdb = createSegmentedMdbList();
    const view = { depthScale: 0.002, pixelSize: 0 };
    const cell = ValueCell.create(pixelLevels());
    mdb.prepare(cell, 1, false, view);
    const on = mdb.list[0];
    const onUniforms = on.uniforms.slice();

    ValueCell.update(cell, Object.assign([] as unknown as LodLevelsValue, { pixelSize: 0 }));
    mdb.prepare(cell, 1, false, view);
    const off = mdb.list[0];
    expect(off).not.toBe(on);
    expect(off.uniforms).toEqual([]);
    // the level drawn with levels on as it was
    expect(on.uniforms).toEqual(onUniforms);

    ValueCell.update(cell, pixelLevels());
    mdb.prepare(cell, 1, false, view);
    expect(mdb.list[0]).not.toBe(off);
    expect(off.uniforms).toEqual([]);
    expect(Array.from(mdb.list[0].uniforms[0][1].ref.value as Vec4)).toEqual([0, 500, 0, 1]);
    // the same uniforms again: the same level
    const again = mdb.list[0];
    mdb.prepare(cell, 1, false, view);
    expect(mdb.list[0]).toBe(again);
  });

  it('sphere levels by pixel size from the strides and scale biases, none at 0', () => {
    const spheres = Spheres.createEmpty();
    const lodLevels = getLodPixelLevels(3);
    const props = { ...PD.getDefaultValues(Spheres.Params), lodPixelSize: 3 };
    const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);

    const value = values.lodLevels.ref.value as LodLevelsValue;
    expect(value.pixelSize).toBe(3);
    // level j ends where its spheres (radius 1.7 Å times stride^(1 / scaleBias)) are 3 px across
    const boundary = (j: number) => (2 * 1.7 * Math.pow(lodLevels[j].stride, 1 / lodLevels[j].scaleBias)) / 3;
    expect(value.map((l) => [l[0], l[1]])).toEqual([
      [0, boundary(0)],
      [boundary(0), boundary(1)],
      [boundary(1), boundary(2)],
      [boundary(2), Infinity],
    ]);

    Spheres.Utils.updateValues(values, { ...props, lodPixelSize: 0 });
    const none = values.lodLevels.ref.value as LodLevelsValue;
    expect(none.pixelSize).toBeUndefined();
    expect(none.length).toBe(0);
  });

  it("keeps coarse items' base through levels off: on again, the levels of the whole", () => {
    const spheres = Spheres.createEmpty();
    const props = { ...PD.getDefaultValues(Spheres.Params), lodPixelSize: 3 };
    const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);
    const count = values.uVertexCount.ref.value;
    ValueCell.update(values.lodLevels, getLodPixelLevelsValue(count, 3, 1, undefined, 10));
    Spheres.Utils.updateValues(values, { ...props, lodPixelSize: 0 });
    expect((values.lodLevels.ref.value as LodLevelsValue).length).toBe(0);
    Spheres.Utils.updateValues(values, props);
    const on = values.lodLevels.ref.value as LodLevelsValue & { base?: number };
    expect(on.base).toBe(10);
    expect(on.map((l) => l[4])).toEqual(getLodPixelLevelsValue(count, 3, 1, undefined, 10).map((l) => l[4]));
  });

  it('coarse items (every 10th of a whole): each level drawn as its stride of the whole, the finer ones as the first', () => {
    const whole = getLodPixelLevelsValue(1000, 3, 1);
    const coarse = getLodPixelLevelsValue(100, 3, 1, undefined, 10);
    expect(coarse.base).toBe(10);
    // the ranges as the whole's: merged with whole members
    expect(coarse.map((l) => [l[0], l[1]])).toEqual(whole.map((l) => [l[0], l[1]]));
    // per level: count (6 an item), scale, stride in the coarse items
    expect(coarse.map((l) => [l[2], l[4]])).toEqual([
      [600, 1],
      [600, 1],
      [60, 10],
      [6, 100],
    ]);
    // the finest level drawn as the stride 10's: its scale
    expect(coarse[0][3]).toBeCloseTo(whole[1][3], 12);
    for (let j = 1; j < 4; ++j) expect(coarse[j][3]).toBeCloseTo(whole[j][3], 12);
    // the base above every stride: the coarsest's scale bias, every item at each level
    const coarser = getLodPixelLevelsValue(1, 3, 1, undefined, 10_000);
    expect(coarser.map((l) => l[4])).toEqual([1, 1, 1, 1]);
    expect(coarser[0][3]).toBeCloseTo(Math.pow(10_000, 1 / 2.3), 12);
  });

  it('presets: pixel sizes of the chain of strides 1, 10, 100, 1000', () => {
    expect(Spheres.LodPixelPresets).toEqual({ performance: 6, balanced: 4, quality: 3, ultra: 1.5 });
    for (const pixelSize of Object.values(Spheres.LodPixelPresets)) {
      expect(getLodPixelLevels(pixelSize).map((l) => l.stride)).toEqual([1, 10, 100, 1000]);
    }
    expect(getLodPixelLevels(0)).toEqual([]);
    // particles by size on screen by default, spheres without levels
    const particles = PD.getDefaultValues(ParticleTargetRepresentationParams);
    expect(particles.lodPixelSize).toBe(3);
    expect(PD.getDefaultValues(Spheres.Params).lodPixelSize).toBe(0);
  });

  it('dots: the boundaries for their own size (lodRadius 0), atoms for the nominal radius', () => {
    const lodLevels = getLodPixelLevels(3);
    const create = (lodRadius: number, size: number) => {
      const props = { ...PD.getDefaultValues(Spheres.Params), lodPixelSize: 3, lodRadius };
      const values = Spheres.Utils.createValuesSimple(Spheres.createEmpty(), props, ColorNames.orange, size);
      return { values, props };
    };
    const boundary = (radius: number) => Spheres.getLodPixelBoundary(lodLevels[0], 3, radius);

    // dots of radius 4: their levels switch where they are 3 px across
    const dots = create(0, 4);
    expect((dots.values.lodLevels.ref.value as LodLevelsValue)[0][1]).toBeCloseTo(boundary(4), 12);
    // atoms: the nominal radius whatever their size
    const atoms = create(PD.getDefaultValues(Spheres.Params).lodRadius, 4);
    expect((atoms.values.lodLevels.ref.value as LodLevelsValue)[0][1]).toBeCloseTo(boundary(1.7), 12);

    // the dots resized: their boundaries follow on the next props update
    dots.values.uSize.ref.value = 2;
    Spheres.Utils.updateValues(dots.values, dots.props);
    expect((dots.values.lodLevels.ref.value as LodLevelsValue)[0][1]).toBeCloseTo(boundary(2), 12);
  });
});
