/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ValueCell } from '@molstar/core/util';
import type { DataArray } from '@molstar/graphics/gl/renderable/util';
import { BaseGeometry } from './base.js';
import { orderPositionGroup } from './center-group.js';

/**
 * Levels of detail of the geometries of one item per quad (spheres, points): items ordered into levels by strides
 * (`orderPositionGroup`), each level drawn in its range of distances (or sizes on screen), its items scaled up
 * to make up for the fewer of them
 */

/** Levels of detail (`lodLevels` param) */
export type LodLevels = {
  minDistance: number;
  maxDistance: number;
  overlap: number;
  stride: number;
  scaleBias: number;
}[];

export function areLodLevelsEqual(a: LodLevels, b: LodLevels) {
  if (a.length !== b.length) return false;
  for (let i = 0, il = a.length; i < il; ++i) {
    if (a[i].maxDistance !== b[i].maxDistance) return false;
    if (a[i].minDistance !== b[i].minDistance) return false;
    if (a[i].overlap !== b[i].overlap) return false;
    if (a[i].stride !== b[i].stride) return false;
    if (a[i].scaleBias !== b[i].scaleBias) return false;
  }
  return true;
}

/**
 * Per level its range (distances in Å, or with `pixelSize` above 0 the Å a CSS pixel covers), fade, vertex count and
 * sphere scale (`renderable.ts` `LodLevelsValue`), then the level's props
 */
export type LodLevelsValue = [
  min: number,
  max: number,
  overlap: number,
  count: number,
  scale: number,
  stride: number,
  scaleBias: number,
  minDistance: number,
  maxDistance: number,
][] & { readonly pixelSize?: number; readonly radius?: number };

/** A carbon's van der Waals radius: the nominal sphere radius levels by size on screen are chosen for (atoms) */
export const LodNominalRadius = 1.7;

/**
 * Where level `j`'s spheres get smaller than `pixelSize` CSS pixels across, in Å a pixel covers: its nominal spheres
 * (`LodNominalRadius` times the stride's scale) `pixelSize` pixels across (`screen-space-lod.md`). The preset's stride:
 * the stride adjusted for `sizeFactor` draws spheres the size of those of the preset's stride at size factor 1
 */
export function getLodPixelBoundary(level: LodLevels[0], pixelSize: number, radius = LodNominalRadius) {
  return (2 * radius * Math.pow(Math.max(1, level.stride), 1 / level.scaleBias)) / pixelSize;
}

/** The levels' value (`LodLevelsValue`) of the levels `prop` ordered by `getStrides` (their vertex `offsets`) */
export function getLodLevelsValue(
  prop: LodLevels,
  sizeFactor: number,
  offsets: number[],
  count: number,
  pixelSize: number,
  radius: number,
): LodLevelsValue {
  const pixel = pixelSize > 0;
  const value = prop.map((l, i): LodLevelsValue[0] => {
    const stride = getAdjustedStride(l, sizeFactor);
    return [
      pixel ? (i === 0 ? 0 : getLodPixelBoundary(prop[i - 1], pixelSize, radius)) : l.minDistance,
      pixel ? (i === prop.length - 1 ? Infinity : getLodPixelBoundary(l, pixelSize, radius)) : l.maxDistance,
      l.overlap,
      offsets[offsets.length - 1 - i],
      Math.pow(Math.min(count, stride), 1 / l.scaleBias),
      l.stride,
      l.scaleBias,
      l.minDistance,
      l.maxDistance,
    ];
  });
  return pixel ? Object.assign(value, { pixelSize, radius }) : value;
}

export function getLodLevels(lodLevelsValue: LodLevelsValue): LodLevels {
  return lodLevelsValue.map((l) => ({
    minDistance: l[7],
    maxDistance: l[8],
    overlap: l[2],
    stride: l[5],
    scaleBias: l[6],
  }));
}

function getAdjustedStride(lodLevel: LodLevels[0], sizeFactor: number) {
  return Math.max(1, Math.round(lodLevel.stride / Math.pow(sizeFactor, lodLevel.scaleBias)));
}

/** The strides of the levels, coarsest first (`orderPositionGroup`), adjusted for `sizeFactor` */
export function getStrides(lodLevels: LodLevels, sizeFactor: number) {
  return lodLevels.map((l) => getAdjustedStride(l, sizeFactor)).reverse();
}

export const LodLevelsPresets: { [key in 'performance' | 'balanced' | 'quality' | 'ultra']: LodLevels } = {
  performance: [
    { minDistance: 1, maxDistance: 300, overlap: 0, stride: 1, scaleBias: 1 },
    { minDistance: 300, maxDistance: 2000, overlap: 0, stride: 40, scaleBias: 3 },
    { minDistance: 2000, maxDistance: 6000, overlap: 0, stride: 150, scaleBias: 3 },
    { minDistance: 6000, maxDistance: 10000000, overlap: 0, stride: 300, scaleBias: 2.5 },
  ],
  balanced: [
    { minDistance: 1, maxDistance: 500, overlap: 0, stride: 1, scaleBias: 1 },
    { minDistance: 500, maxDistance: 2000, overlap: 0, stride: 15, scaleBias: 3 },
    { minDistance: 2000, maxDistance: 6000, overlap: 0, stride: 70, scaleBias: 2.7 },
    { minDistance: 6000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2.5 },
  ],
  quality: [
    { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1, scaleBias: 1 },
    { minDistance: 1000, maxDistance: 4000, overlap: 0, stride: 10, scaleBias: 3 },
    { minDistance: 4000, maxDistance: 10000, overlap: 0, stride: 50, scaleBias: 2.7 },
    { minDistance: 10000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2.3 },
  ],
  ultra: [
    { minDistance: 1, maxDistance: 5000, overlap: 0, stride: 1, scaleBias: 1 },
    { minDistance: 5000, maxDistance: 10000, overlap: 0, stride: 10, scaleBias: 3 },
    { minDistance: 10000, maxDistance: 30000, overlap: 0, stride: 50, scaleBias: 2.5 },
    { minDistance: 30000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2 },
  ],
};

/** The `lodLevels` param: levels by distance, none by default, with the presets of `LodLevelsPresets` */
export function getLodLevelsParam() {
  return PD.ObjectList(
    {
      minDistance: PD.Numeric(0),
      maxDistance: PD.Numeric(0),
      overlap: PD.Numeric(0),
      stride: PD.Numeric(0),
      scaleBias: PD.Numeric(3, { min: 0.1, max: 10, step: 0.1 }),
    },
    (o) => `${o.stride}`,
    {
      ...BaseGeometry.CullingLodCategory,
      defaultValue: [] as LodLevels,
      presets: Object.entries(LodLevelsPresets).map(([k, v]) => [v, k]),
    },
  );
}

/**
 * Orders the `count` items of `positionGroup` into the levels `lodLevels` (in place, `orderPositionGroup`), their
 * strides adjusted for `sizeFactor`; returns the levels' value (empty without levels). `pixelSize` and `radius`: levels
 * by size on screen (`getLodPixelBoundary`)
 */
export function orderLodLevels(
  positionGroup: ValueCell<DataArray<Float32Array>>,
  count: number,
  lodLevels: LodLevels,
  sizeFactor: number,
  pixelSize = 0,
  radius = 0,
): LodLevelsValue {
  const offsets = orderPositionGroup(positionGroup, count, getStrides(lodLevels, sizeFactor));
  return offsets ? getLodLevelsValue(lodLevels, sizeFactor, offsets, count, pixelSize, radius) : [];
}
