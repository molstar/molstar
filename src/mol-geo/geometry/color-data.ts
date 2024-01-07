/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ValueCell } from '../../mol-util';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { Color } from '../../mol-util/color';
import { Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { LocationIterator } from '../util/location-iterator';
import { NullLocation } from '../../mol-model/location';
import { LocationColor, ColorTheme, ColorVolume } from '../../mol-theme/color';
import { createNullTexture, Texture } from '../../mol-gl/webgl/texture';

export type ColorTypeLocation = 'uniform' | 'instance' | 'group' | 'groupInstance' | 'vertex' | 'vertexInstance';
export type ColorTypeGrid = 'volume' | 'volumeInstance';
export type ColorTypeDirect = 'direct';
export type ColorType = ColorTypeLocation | ColorTypeGrid | ColorTypeDirect;

export type ColorData = {
    uColor: ValueCell<Vec3>,
    tColor: ValueCell<TextureImage<Uint8Array>>,
    tColorGrid: ValueCell<Texture>,
    tPalette: ValueCell<TextureImage<Uint8Array>>,
    uColorTexDim: ValueCell<Vec2>,
    uColorGridDim: ValueCell<Vec3>,
    uColorGridTransform: ValueCell<Vec4>,
    dColorType: ValueCell<string>,
    dUsePalette: ValueCell<boolean>,
}

export function createColors(locationIt: LocationIterator, positionIt: LocationIterator, colorTheme: ColorTheme<any, any>, colorData?: ColorData): ColorData {
    const data = _createColors(locationIt, positionIt, colorTheme, colorData);
    if (colorTheme.palette) {
        ValueCell.updateIfChanged(data.dUsePalette, true);
        updatePaletteTexture(colorTheme.palette, data.tPalette);
    } else {
        ValueCell.updateIfChanged(data.dUsePalette, false);
    }
    return data;
}

function _createColors(locationIt: LocationIterator, positionIt: LocationIterator, colorTheme: ColorTheme<any, any>, colorData?: ColorData): ColorData {
    switch (colorTheme.granularity) {
        case 'uniform': return createUniformColor(locationIt, colorTheme.color, colorData);
        case 'instance':
            return locationIt.nonInstanceable
                ? createGroupColor(locationIt, colorTheme.color, colorData)
                : createInstanceColor(locationIt, colorTheme.color, colorData);
        case 'group': return createGroupColor(locationIt, colorTheme.color, colorData);
        case 'groupInstance': return createGroupInstanceColor(locationIt, colorTheme.color, colorData);
        case 'vertex': return createVertexColor(positionIt, colorTheme.color, colorData);
        case 'vertexInstance': return createVertexInstanceColor(positionIt, colorTheme.color, colorData);
        case 'volume': return createGridColor(colorTheme.grid, 'volume', colorData);
        case 'volumeInstance': return createGridColor(colorTheme.grid, 'volumeInstance', colorData);
        case 'direct': return createDirectColor(colorData);
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
        }
    }

    if (isSynced) return;

    const array = new Uint8Array(palette.colors.length * 3);
    let o = 0;
    for (const c of palette.colors) {
        const [r, g, b] = Color.toRgb(c);
        array[o++] = r;
        array[o++] = g;
        array[o++] = b;
    }

    ValueCell.update(cell, { array, height: 1, width: palette.colors.length, filter: palette.filter });
}

//

export function createValueColor(value: Color, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.uColor, Color.toVec3Normalized(colorData.uColor.ref.value, value));
        ValueCell.updateIfChanged(colorData.dColorType, 'uniform');
        return colorData;
    } else {
        return {
            uColor: ValueCell.create(Color.toVec3Normalized(Vec3(), value)),
            tColor: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            tColorGrid: ValueCell.create(createNullTexture()),
            tPalette: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            uColorTexDim: ValueCell.create(Vec2.create(1, 1)),
            uColorGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uColorGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dColorType: ValueCell.create('uniform'),
            dUsePalette: ValueCell.create(false),
        };
    }
}

/** Creates color uniform */
function createUniformColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    return createValueColor(color(NullLocation, false), colorData);
}

//

export function createTextureColor(colors: TextureImage<Uint8Array>, type: ColorType, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.tColor, colors);
        ValueCell.update(colorData.uColorTexDim, Vec2.create(colors.width, colors.height));
        ValueCell.updateIfChanged(colorData.dColorType, type);
        return colorData;
    } else {
        return {
            uColor: ValueCell.create(Vec3()),
            tColor: ValueCell.create(colors),
            tColorGrid: ValueCell.create(createNullTexture()),
            tPalette: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            uColorTexDim: ValueCell.create(Vec2.create(colors.width, colors.height)),
            uColorGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uColorGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dColorType: ValueCell.create(type),
            dUsePalette: ValueCell.create(false),
        };
    }
}

/** Creates color texture with color for each instance */
function createInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { instanceCount } = locationIt;
    const colors = createTextureImage(Math.max(1, instanceCount), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext) {
        const { location, isSecondary, instanceIndex } = locationIt.move();
        Color.toArray(color(location, isSecondary), colors.array, instanceIndex * 3);
        locationIt.skipInstance();
    }
    return createTextureColor(colors, 'instance', colorData);
}

/** Creates color texture with color for each group (i.e. shared across instances) */
function createGroupColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount, hasLocation2 } = locationIt;
    const colors = createTextureImage(Math.max(1, groupCount * (hasLocation2 ? 2 : 1)), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    const indexMultiplier = hasLocation2 ? 6 : 3;
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const { location, location2, isSecondary, groupIndex } = locationIt.move();
        Color.toArray(color(location, isSecondary), colors.array, groupIndex * indexMultiplier);
        if (hasLocation2) Color.toArray(color(location2, isSecondary), colors.array, groupIndex * indexMultiplier + 3);
    }
    return createTextureColor(colors, 'group', colorData);
}

/** Creates color texture with color for each group in each instance */
function createGroupInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount, instanceCount, hasLocation2 } = locationIt;
    const count = instanceCount * groupCount * (hasLocation2 ? 2 : 1);
    const colors = createTextureImage(Math.max(1, count), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    const indexMultiplier = hasLocation2 ? 6 : 3;
    while (locationIt.hasNext) {
        const { location, location2, isSecondary, index } = locationIt.move();
        Color.toArray(color(location, isSecondary), colors.array, index * indexMultiplier);
        if (hasLocation2) Color.toArray(color(location2, isSecondary), colors.array, index * indexMultiplier + 3);
    }
    return createTextureColor(colors, 'groupInstance', colorData);
}

/** Creates color texture with color for each vertex (i.e. shared across instances) */
function createVertexColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount, stride } = locationIt;
    const colors = createTextureImage(Math.max(1, groupCount), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    locationIt.voidInstances();
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const { location, isSecondary, groupIndex } = locationIt.move();
        const c = color(location, isSecondary);
        for (let i = 0; i < stride; ++i) {
            Color.toArray(c, colors.array, (groupIndex + i) * 3);
        }
    }
    return createTextureColor(colors, 'vertex', colorData);
}

/** Creates color texture with color for each vertex in each instance */
function createVertexInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount, instanceCount, stride } = locationIt;
    const count = instanceCount * groupCount;
    const colors = createTextureImage(Math.max(1, count), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext) {
        const { location, isSecondary, index } = locationIt.move();
        const c = color(location, isSecondary);
        for (let i = 0; i < stride; ++i) {
            Color.toArray(c, colors.array, (index + i) * 3);
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
        ValueCell.update(colorData.uColorTexDim, Vec2.create(width, height));
        ValueCell.update(colorData.uColorGridDim, Vec3.clone(dimension));
        ValueCell.update(colorData.uColorGridTransform, Vec4.clone(transform));
        ValueCell.updateIfChanged(colorData.dColorType, type);
        return colorData;
    } else {
        return {
            uColor: ValueCell.create(Vec3()),
            tColor: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            tColorGrid: ValueCell.create(colors),
            tPalette: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            uColorTexDim: ValueCell.create(Vec2.create(width, height)),
            uColorGridDim: ValueCell.create(Vec3.clone(dimension)),
            uColorGridTransform: ValueCell.create(Vec4.clone(transform)),
            dColorType: ValueCell.create(type),
            dUsePalette: ValueCell.create(false),
        };
    }
}

//

/** Creates direct color */
function createDirectColor(colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.updateIfChanged(colorData.dColorType, 'direct');
        return colorData;
    } else {
        return {
            uColor: ValueCell.create(Vec3()),
            tColor: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            tColorGrid: ValueCell.create(createNullTexture()),
            tPalette: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            uColorTexDim: ValueCell.create(Vec2.create(1, 1)),
            uColorGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uColorGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dColorType: ValueCell.create('direct'),
            dUsePalette: ValueCell.create(false),
        };
    }
}
