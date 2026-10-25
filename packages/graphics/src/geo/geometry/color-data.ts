/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ValueCell } from '@molstar/core/util';
import { type DataArray, type TextureImage, createDataArray } from '@molstar/graphics/gl/renderable/util';
import { Color } from '@molstar/core/util/color';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import type { LocationIterator } from '../util/location-iterator.js';
import { NullLocation } from '@molstar/model/model/location';
import type { LocationColor, ColorTheme, ColorVolume } from '@molstar/graphics/theme/color';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';

export type ColorTypeLocation = 'uniform' | 'instance' | 'group' | 'groupInstance' | 'vertex' | 'vertexInstance';
export type ColorTypeGrid = 'volume' | 'volumeInstance';
export type ColorTypeDirect = 'direct';
export type ColorType = ColorTypeLocation | ColorTypeGrid | ColorTypeDirect;

export type ColorData = {
  uColor: ValueCell<Vec3>;
  tColor: ValueCell<DataArray<Uint8Array>>;
  tColorGrid: ValueCell<Texture>;
  uPaletteDomain: ValueCell<Vec2>;
  uPaletteDefault: ValueCell<Vec3>;
  tPalette: ValueCell<TextureImage<Uint8Array>>;
  uColorGridDim: ValueCell<Vec3>;
  /** the size of the 2D layout of the grid (`tColorGrid`) */
  uColorGridTexDim: ValueCell<Vec2>;
  uColorGridTransform: ValueCell<Vec4>;
  uColorIndexType: ValueCell<number>;
  dColorType: ValueCell<string>;
  /** the `ColorType` */
  colorType: ValueCell<string>;
  dUsePalette: ValueCell<boolean>;
};

/**
 * How the color types that share a define read their data (`uColorIndexType`): one program for each of 'texture' (the
 * index of `tColor`) and 'grid' (the position `tColorGrid` is sampled at), the shaders select at runtime
 */
export const ColorTypeIndex: { readonly [T in ColorType]?: number } = {
  instance: 0,
  group: 1,
  groupInstance: 2,
  vertex: 3,
  vertexInstance: 4,
  volume: 0,
  volumeInstance: 1,
};

const ColorTypeDefines: { readonly [T in ColorType]: string } = {
  uniform: 'uniform',
  instance: 'texture',
  group: 'texture',
  groupInstance: 'texture',
  vertex: 'texture',
  vertexInstance: 'texture',
  volume: 'grid',
  volumeInstance: 'grid',
  direct: 'direct',
};

/**
 * The define value of the color type `type`: the types read from `tColor` by index share 'texture', the grid types
 * share 'grid' (see `ColorTypeIndex`)
 */
export function getColorTypeDefine(type: ColorType | string): string {
  return ColorTypeDefines[type as ColorType] ?? type;
}

/** Sets the color type `type` (a `ColorType`), its define and index (`ColorTypeIndex`) of `values` */
export function setColorType(
  values: Pick<ColorData, 'colorType' | 'dColorType' | 'uColorIndexType'>,
  type: ColorType | string,
) {
  ValueCell.updateIfChanged(values.colorType, type);
  ValueCell.updateIfChanged(values.dColorType, getColorTypeDefine(type));
  ValueCell.updateIfChanged(values.uColorIndexType, ColorTypeIndex[type as ColorType] ?? 0);
}

/** The values of the color type `type` (its define and index, see `setColorType`) */
function createColorType(type: ColorType) {
  return {
    uColorIndexType: ValueCell.create(ColorTypeIndex[type] ?? 0),
    dColorType: ValueCell.create(getColorTypeDefine(type)),
    colorType: ValueCell.create<string>(type),
  };
}

/** Whether colors of `type` are sampled at vertex positions (the position iterator of `createColors`) */
export function isPositionColorType(type: ColorType) {
  return type === 'vertex' || type === 'vertexInstance';
}

export function createColors(
  locationIt: LocationIterator,
  positionIt: LocationIterator,
  colorTheme: ColorTheme<any, any>,
  colorData?: ColorData,
): ColorData {
  const data = _createColors(locationIt, positionIt, colorTheme, colorData);
  if (colorTheme.palette) {
    ValueCell.updateIfChanged(data.dUsePalette, true);
    const [min, max] = colorTheme.palette.domain || [0, 1];
    ValueCell.update(data.uPaletteDomain, Vec2.set(data.uPaletteDomain.ref.value, min, max));
    ValueCell.update(
      data.uPaletteDefault,
      Color.toVec3Normalized(data.uPaletteDefault.ref.value, colorTheme.palette.defaultColor ?? Color(0xcccccc)),
    );
    updatePaletteTexture(colorTheme.palette, data.tPalette);
  } else {
    ValueCell.updateIfChanged(data.dUsePalette, false);
  }
  return data;
}

function _createColors(
  locationIt: LocationIterator,
  positionIt: LocationIterator,
  colorTheme: ColorTheme<any, any>,
  colorData?: ColorData,
): ColorData {
  switch (colorTheme.granularity) {
    case 'uniform':
      return createUniformColor(locationIt, colorTheme.color, colorData);
    case 'instance':
      return locationIt.nonInstanceable
        ? createGroupColor(locationIt, colorTheme.color, colorData)
        : createInstanceColor(locationIt, colorTheme.color, colorData);
    case 'group':
      return createGroupColor(locationIt, colorTheme.color, colorData);
    case 'groupInstance':
      return createGroupInstanceColor(locationIt, colorTheme.color, colorData);
    case 'vertex':
      return createVertexColor(positionIt, colorTheme.color, colorData);
    case 'vertexInstance':
      return createVertexInstanceColor(positionIt, colorTheme.color, colorData);
    case 'volume':
      return createGridColor(colorTheme.grid, 'volume', colorData);
    case 'volumeInstance':
      return createGridColor(colorTheme.grid, 'volumeInstance', colorData);
    case 'direct':
      return createDirectColor(colorData);
  }
}

function updatePaletteTexture(palette: ColorTheme.Palette, cell: ValueCell<TextureImage<Uint8Array>>) {
  let isSynced = true;
  const texture = cell.ref.value;
  if (palette.colors.length !== texture.width || texture.filter !== palette.filter) {
    isSynced = false;
  } else {
    const data = texture.array;
    let o = 0;
    for (const c of palette.colors) {
      const [r, g, b] = Color.toRgb(c);
      if (data[o++] !== r || data[o++] !== g || data[o++] !== b) {
        isSynced = false;
        break;
      }
      o += 1;
    }
  }

  if (isSynced) return;

  const array = new Uint8Array(palette.colors.length * 4);
  let o = 0;
  for (const c of palette.colors) {
    const [r, g, b] = Color.toRgb(c);
    array[o++] = r;
    array[o++] = g;
    array[o++] = b;
    array[o++] = 255;
  }

  ValueCell.update(cell, { array, height: 1, width: palette.colors.length, filter: palette.filter });
}

//

export function createValueColor(value: Color, colorData?: ColorData): ColorData {
  if (colorData) {
    ValueCell.update(colorData.uColor, Color.toVec3Normalized(colorData.uColor.ref.value, value));
    setColorType(colorData, 'uniform');
    return colorData;
  } else {
    return {
      uColor: ValueCell.create(Color.toVec3Normalized(Vec3(), value)),
      tColor: ValueCell.create({ array: new Uint8Array(4), count: 1 }),
      tColorGrid: ValueCell.create(createNullTexture()),
      uPaletteDomain: ValueCell.create(Vec2.create(0, 1)),
      uPaletteDefault: ValueCell.create(Vec3()),
      tPalette: ValueCell.create({ array: new Uint8Array(4), width: 1, height: 1 }),
      uColorGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uColorGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uColorGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      ...createColorType('uniform'),
      dUsePalette: ValueCell.create(false),
    };
  }
}

/** Creates color uniform */
function createUniformColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
  locationIt.reset();
  const loc = locationIt.hasNext ? locationIt.move() : { location: NullLocation, isSecondary: false };
  return createValueColor(color(loc.location, loc.isSecondary), colorData);
}

//

export function createTextureColor(colors: DataArray<Uint8Array>, type: ColorType, colorData?: ColorData): ColorData {
  if (colorData) {
    ValueCell.update(colorData.tColor, colors);
    setColorType(colorData, type);
    return colorData;
  } else {
    return {
      uColor: ValueCell.create(Vec3()),
      tColor: ValueCell.create(colors),
      tColorGrid: ValueCell.create(createNullTexture()),
      uPaletteDomain: ValueCell.create(Vec2.create(0, 1)),
      uPaletteDefault: ValueCell.create(Vec3()),
      tPalette: ValueCell.create({ array: new Uint8Array(4), width: 1, height: 1 }),
      uColorGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uColorGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uColorGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      ...createColorType(type),
      dUsePalette: ValueCell.create(false),
    };
  }
}

/** Writes `color` as opaque rgba */
function setColor(array: Uint8Array, offset: number, color: Color) {
  Color.toArray(color, array, offset);
  array[offset + 3] = 255;
}

/** Creates color texture with color for each instance */
function createInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
  const { instanceCount } = locationIt;
  const colors = createDataArray(
    Math.max(1, instanceCount),
    4,
    Uint8Array,
    colorData && colorData.tColor.ref.value.array,
  );
  // by index, not by `move`: a location per instance (its first group's), no iterator state to keep per step (large
  // models: millions of instances); reset as for the loop it replaces (no instances voided)
  locationIt.reset();
  if (locationIt.groupCount > 0) {
    for (let i = 0; i < instanceCount; ++i) {
      setColor(colors.array, i * 4, color(locationIt.getInstanceLocation(i), locationIt.isInstanceSecondary(i)));
    }
  }
  return createTextureColor(colors, 'instance', colorData);
}

/** Creates color texture with color for each group (i.e. shared across instances) */
function createGroupColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
  const { groupCount, hasLocation2 } = locationIt;
  const colors = createDataArray(
    Math.max(1, groupCount * (hasLocation2 ? 2 : 1)),
    4,
    Uint8Array,
    colorData && colorData.tColor.ref.value.array,
  );
  locationIt.reset();
  const indexMultiplier = hasLocation2 ? 8 : 4;
  while (locationIt.hasNext && !locationIt.isNextNewInstance) {
    const { location, location2, isSecondary, groupIndex } = locationIt.move();
    setColor(colors.array, groupIndex * indexMultiplier, color(location, isSecondary));
    if (hasLocation2) setColor(colors.array, groupIndex * indexMultiplier + 4, color(location2, isSecondary));
  }
  return createTextureColor(colors, 'group', colorData);
}

/** Creates color texture with color for each group in each instance */
function createGroupInstanceColor(
  locationIt: LocationIterator,
  color: LocationColor,
  colorData?: ColorData,
): ColorData {
  const { groupCount, instanceCount, hasLocation2 } = locationIt;
  const count = instanceCount * groupCount * (hasLocation2 ? 2 : 1);
  const colors = createDataArray(Math.max(1, count), 4, Uint8Array, colorData && colorData.tColor.ref.value.array);
  locationIt.reset();
  const indexMultiplier = hasLocation2 ? 8 : 4;
  while (locationIt.hasNext) {
    const { location, location2, isSecondary, index } = locationIt.move();
    setColor(colors.array, index * indexMultiplier, color(location, isSecondary));
    if (hasLocation2) setColor(colors.array, index * indexMultiplier + 4, color(location2, isSecondary));
  }
  return createTextureColor(colors, 'groupInstance', colorData);
}

/** Creates color texture with color for each vertex (i.e. shared across instances) */
function createVertexColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
  const { groupCount, stride } = locationIt;
  const colors = createDataArray(Math.max(1, groupCount), 4, Uint8Array, colorData && colorData.tColor.ref.value.array);
  locationIt.reset();
  locationIt.voidInstances();
  while (locationIt.hasNext && !locationIt.isNextNewInstance) {
    const { location, isSecondary, groupIndex } = locationIt.move();
    const c = color(location, isSecondary);
    for (let i = 0; i < stride; ++i) {
      setColor(colors.array, (groupIndex + i) * 4, c);
    }
  }
  return createTextureColor(colors, 'vertex', colorData);
}

/** Creates color texture with color for each vertex in each instance */
function createVertexInstanceColor(
  locationIt: LocationIterator,
  color: LocationColor,
  colorData?: ColorData,
): ColorData {
  const { groupCount, instanceCount, stride } = locationIt;
  const count = instanceCount * groupCount;
  const colors = createDataArray(Math.max(1, count), 4, Uint8Array, colorData && colorData.tColor.ref.value.array);
  locationIt.reset();
  while (locationIt.hasNext) {
    const { location, isSecondary, index } = locationIt.move();
    const c = color(location, isSecondary);
    for (let i = 0; i < stride; ++i) {
      setColor(colors.array, (index + i) * 4, c);
    }
  }
  return createTextureColor(colors, 'vertexInstance', colorData);
}

//

export function createGridColor(grid: ColorVolume, type: ColorType, colorData?: ColorData): ColorData {
  const { colors, dimension, transform } = grid;
  const width = colors.getWidth();
  const height = colors.getHeight();
  if (colorData) {
    ValueCell.update(colorData.tColorGrid, colors);
    ValueCell.update(colorData.uColorGridTexDim, Vec2.create(width, height));
    ValueCell.update(colorData.uColorGridDim, Vec3.clone(dimension));
    ValueCell.update(colorData.uColorGridTransform, Vec4.clone(transform));
    setColorType(colorData, type);
    return colorData;
  } else {
    return {
      uColor: ValueCell.create(Vec3()),
      tColor: ValueCell.create({ array: new Uint8Array(4), count: 1 }),
      tColorGrid: ValueCell.create(colors),
      uPaletteDomain: ValueCell.create(Vec2.create(0, 1)),
      uPaletteDefault: ValueCell.create(Vec3()),
      tPalette: ValueCell.create({ array: new Uint8Array(4), width: 1, height: 1 }),
      uColorGridTexDim: ValueCell.create(Vec2.create(width, height)),
      uColorGridDim: ValueCell.create(Vec3.clone(dimension)),
      uColorGridTransform: ValueCell.create(Vec4.clone(transform)),
      ...createColorType(type),
      dUsePalette: ValueCell.create(false),
    };
  }
}

//

/** Creates direct color */
function createDirectColor(colorData?: ColorData): ColorData {
  if (colorData) {
    setColorType(colorData, 'direct');
    return colorData;
  } else {
    return {
      uColor: ValueCell.create(Vec3()),
      tColor: ValueCell.create({ array: new Uint8Array(4), count: 1 }),
      tColorGrid: ValueCell.create(createNullTexture()),
      uPaletteDomain: ValueCell.create(Vec2.create(0, 1)),
      uPaletteDefault: ValueCell.create(Vec3()),
      tPalette: ValueCell.create({ array: new Uint8Array(4), width: 1, height: 1 }),
      uColorGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uColorGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uColorGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      ...createColorType('direct'),
      dUsePalette: ValueCell.create(false),
    };
  }
}
