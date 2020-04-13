/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { Color } from '../../mol-util/color';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { LocationIterator } from '../util/location-iterator';
import { NullLocation } from '../../mol-model/location';
import { LocationColor, ColorTheme } from '../../mol-theme/color';
import { Geometry } from './geometry';

export type ColorType = 'uniform' | 'instance' | 'group' | 'groupInstance'

export type ColorData = {
    uColor: ValueCell<Vec3>,
    tColor: ValueCell<TextureImage<Uint8Array>>,
    uColorTexDim: ValueCell<Vec2>,
    dColorType: ValueCell<string>,
}

export function createColors(locationIt: LocationIterator, colorTheme: ColorTheme<any>, colorData?: ColorData): ColorData {
    switch (Geometry.getGranularity(locationIt, colorTheme.granularity)) {
        case 'uniform': return createUniformColor(locationIt, colorTheme.color, colorData);
        case 'group': return createGroupColor(locationIt, colorTheme.color, colorData);
        case 'groupInstance': return createGroupInstanceColor(locationIt, colorTheme.color, colorData);
        case 'instance': return createInstanceColor(locationIt, colorTheme.color, colorData);
    }
}

export function createValueColor(value: Color, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.uColor, Color.toVec3Normalized(colorData.uColor.ref.value, value));
        if (colorData.dColorType.ref.value !== 'uniform') {
            ValueCell.update(colorData.dColorType, 'uniform');
        }
        return colorData;
    } else {
        return {
            uColor: ValueCell.create(Color.toVec3Normalized(Vec3(), value)),
            tColor: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            uColorTexDim: ValueCell.create(Vec2.create(1, 1)),
            dColorType: ValueCell.create('uniform'),
        };
    }
}

/** Creates color uniform */
export function createUniformColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    return createValueColor(color(NullLocation, false), colorData);
}

export function createTextureColor(colors: TextureImage<Uint8Array>, type: ColorType, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.tColor, colors);
        ValueCell.update(colorData.uColorTexDim, Vec2.create(colors.width, colors.height));
        if (colorData.dColorType.ref.value !== type) {
            ValueCell.update(colorData.dColorType, type);
        }
        return colorData;
    } else {
        return {
            uColor: ValueCell.create(Vec3.zero()),
            tColor: ValueCell.create(colors),
            uColorTexDim: ValueCell.create(Vec2.create(colors.width, colors.height)),
            dColorType: ValueCell.create(type),
        };
    }
}

/** Creates color texture with color for each instance/unit */
export function createInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
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

/** Creates color texture with color for each group (i.e. shared across instances/units) */
export function createGroupColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount } = locationIt;
    const colors = createTextureImage(Math.max(1, groupCount), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const { location, isSecondary, groupIndex } = locationIt.move();
        Color.toArray(color(location, isSecondary), colors.array, groupIndex * 3);
    }
    return createTextureColor(colors, 'group', colorData);
}

/** Creates color texture with color for each group in each instance (i.e. for each unit) */
export function createGroupInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount, instanceCount } = locationIt;
    const count = instanceCount * groupCount;
    const colors = createTextureImage(Math.max(1, count), 3, Uint8Array, colorData && colorData.tColor.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext) {
        const { location, isSecondary, index } = locationIt.move();
        Color.toArray(color(location, isSecondary), colors.array, index * 3);
    }
    return createTextureColor(colors, 'groupInstance', colorData);
}