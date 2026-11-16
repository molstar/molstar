/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { BaseGeometry } from './base.js';

/**
 * Levels of detail of the geometries of one item per quad (spheres, points), by size on screen (`lodPixelSize`):
 * the levels of a chain of nested strides, level `j` every `stride`-th item (the items in their
 * order, the shaders read item `p * stride` at the level's quad `p`), each level drawn where its spheres are at least the pixel size across, its items scaled up to make up for the fewer of them
 */

/** Levels of detail: a level's stride and scale bias (its distances unused: the levels are by size on screen) */
export type LodLevels = {
  minDistance: number;
  maxDistance: number;
  stride: number;
  scaleBias: number;
}[];

export function areLodLevelsEqual(a: LodLevels, b: LodLevels) {
  if (a.length !== b.length) return false;
  for (let i = 0, il = a.length; i < il; ++i) {
    if (a[i].maxDistance !== b[i].maxDistance) return false;
    if (a[i].minDistance !== b[i].minDistance) return false;
    if (a[i].stride !== b[i].stride) return false;
    if (a[i].scaleBias !== b[i].scaleBias) return false;
  }
  return true;
}

/**
 * Per level its range (distances in Å, or with `pixelSize` above 0 the Å a CSS pixel covers), vertex count (6 an item: the
 * level draws its first `count / 6` quads, quad `p` item `p * stride`), sphere scale and stride (adjusted for the size
 * factor and nested, `getNestedStrides`); `renderable.ts` `LodLevelsValue`
 */
export type LodLevelsValue = [min: number, max: number, count: number, scale: number, stride: number][] & {
  readonly pixelSize?: number;
  readonly radius?: number;
  /** the items every `base`-th of a whole's (`getLodLevelsValue`), 1 without */
  readonly base?: number;
};

/** A carbon's van der Waals radius: the nominal sphere radius levels by size on screen are chosen for (atoms) */
export const LodNominalRadius = 1.7;

/**
 * Where level `j`'s spheres get smaller than `pixelSize` CSS pixels across, in Å a pixel covers: its nominal spheres
 * (`LodNominalRadius` times the stride's scale) `pixelSize` pixels across. The preset's stride:
 * the stride adjusted for `sizeFactor` draws spheres the size of those of the preset's stride at size factor 1
 */
export function getLodPixelBoundary(level: LodLevels[0], pixelSize: number, radius = LodNominalRadius) {
  return (2 * radius * Math.pow(Math.max(1, level.stride), 1 / level.scaleBias)) / pixelSize;
}

/**
 * The levels' value (`LodLevelsValue`) of the levels `prop` for `count` items: level `j` draws the items its nested
 * stride takes (`ceil(count / stride)`). `base` above 1: the items every `base`-th of a whole's (`count` times
 * `base`, a coarse read), each level drawn as its stride of the whole is, a level finer than `base` as the first one
 * that is not (its scale; its items all, its stride 1). Its levels' ranges as the whole's: a coarse member merges with
 * whole ones (`getMergeKey` compares the ranges), its own counts, scales and strides by member
 */
export function getLodLevelsValue(
  prop: LodLevels,
  sizeFactor: number,
  count: number,
  pixelSize: number,
  radius: number,
  base = 1,
): LodLevelsValue {
  const pixel = pixelSize > 0;
  const strides = getNestedStrides(prop, sizeFactor);
  base = Math.max(1, base);
  const whole = count * base;
  // the scale bias of a level finer than `base`: the first not finer's, else the coarsest's
  let coarse = strides.findIndex((s) => s >= base);
  if (coarse < 0) coarse = prop.length - 1;
  const value = prop.map((l, i): LodLevelsValue[0] => {
    const of = Math.max(strides[i], base);
    const stride = Math.max(1, Math.round(of / base));
    const scaleBias = strides[i] >= base ? l.scaleBias : prop[coarse].scaleBias;
    return [
      pixel ? (i === 0 ? 0 : getLodPixelBoundary(prop[i - 1], pixelSize, radius)) : l.minDistance,
      pixel ? (i === prop.length - 1 ? Infinity : getLodPixelBoundary(l, pixelSize, radius)) : l.maxDistance,
      Math.ceil(count / stride) * 6,
      Math.pow(Math.min(whole, of), 1 / scaleBias),
      stride,
    ];
  });
  const extra = base > 1 ? { base } : {};
  return pixel ? Object.assign(value, { pixelSize, radius }, extra) : Object.assign(value, extra);
}

function getAdjustedStride(lodLevel: LodLevels[0], sizeFactor: number) {
  return Math.max(1, Math.round(lodLevel.stride / Math.pow(sizeFactor, lodLevel.scaleBias)));
}

/**
 * The strides of the levels, finest first, adjusted for `sizeFactor` and nested: each a multiple of the one before (the
 * finer), the adjusted stride rounded to one (at least the finer's). Nested levels are subsets of each other, so the
 * items drawn at a level are every stride-th: a coarser level's items a subset of a finer's
 */
export function getNestedStrides(lodLevels: LodLevels, sizeFactor: number) {
  const strides = lodLevels.map((l) => getAdjustedStride(l, sizeFactor));
  for (let k = 1; k < strides.length; ++k) {
    const finer = strides[k - 1];
    strides[k] = Math.max(finer, Math.round(strides[k] / finer) * finer);
  }
  return strides;
}

/**
 * Strides growing by 10 (scale biases 1, 3, 2.7, 2.3): levels by size on screen that reach a whole mesoscale scene with
 * four levels. Nested (each a multiple of the finer), as the size factor's adjustment keeps
 * them (`getNestedStrides`)
 */
const LodPixelChain = [
  [1, 1],
  [10, 3],
  [100, 2.7],
  [1000, 2.3],
] as const;

/** The levels by size on screen at `pixelSize` (`LodPixelChain`), none at 0 */
export function getLodPixelLevels(pixelSize: number): LodLevels {
  if (!(pixelSize > 0)) return [];
  return LodPixelChain.map(([stride, scaleBias]) => ({
    minDistance: 0,
    maxDistance: 0,
    stride,
    scaleBias,
  }));
}

/** The number of levels by size on screen at `pixelSize` (`getLodPixelLevels`), 0 at 0 */
export function getLodPixelLevelCount(pixelSize: number) {
  return pixelSize > 0 ? LodPixelChain.length : 0;
}

/** The pixel sizes of the graphics modes (`lodPixelSize`): the smallest sphere drawn, in CSS pixels */
export const LodPixelPresets = { performance: 6, balanced: 4, quality: 3, ultra: 1.5 } as const;
export type LodPixelPreset = keyof typeof LodPixelPresets;

/** The `lodPixelSize` param: levels of detail by the size of their spheres on screen, none at 0 */
export function getLodPixelSizeParam(defaultValue = 0) {
  return PD.Numeric(
    defaultValue,
    { min: 0, max: 10, step: 0.5 },
    {
      ...BaseGeometry.CullingLodCategory,
      description:
        'If above zero, the levels of detail: the smallest sphere drawn, in CSS pixels (a level of fewer, larger spheres where its spheres are smaller). Zero for none.',
    },
  );
}

/**
 * The levels by size on screen at `pixelSize` (none at 0) of `count` items, their strides adjusted for `sizeFactor`:
 * the levels' value (empty without levels). The items keep their order: a level draws every stride-th. `radius`: the sphere radius the pixel boundaries assume (`getLodPixelBoundary`); `base`: the items every `base`-th of a whole's (`getLodLevelsValue`)
 */
export function getLodPixelLevelsValue(
  count: number,
  pixelSize: number,
  sizeFactor: number,
  radius = LodNominalRadius,
  base = 1,
): LodLevelsValue {
  const lodLevels = getLodPixelLevels(pixelSize);
  if (lodLevels.length > 0) return getLodLevelsValue(lodLevels, sizeFactor, count, pixelSize, radius, base);
  // none, its base kept (coarse sites, every `base`-th of a whole): levels on again are of the whole still
  return base > 1 ? Object.assign([], { base }) : [];
}
