/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';
import { Color } from 'mol-util/color';
import { Vec2, Vec3 } from 'mol-math/linear-algebra';
import { LocationIterator } from '../util/location-iterator';
import { NullLocation } from 'mol-model/location';
import { LocationColor, ColorTheme } from 'mol-theme/color';
import { RuntimeContext } from 'mol-task';
import { getGranularity } from './geometry';

export type ColorType = 'uniform' | 'instance' | 'group' | 'groupInstance'

export type ColorData = {
    uColor: ValueCell<Vec3>,
    aColor: ValueCell<Float32Array>,
    tColor: ValueCell<TextureImage<Uint8Array>>,
    uColorTexDim: ValueCell<Vec2>,
    dColorType: ValueCell<string>,
}

export function createColors(ctx: RuntimeContext, locationIt: LocationIterator, colorTheme: ColorTheme, colorData?: ColorData): Promise<ColorData> {
    switch (getGranularity(locationIt, colorTheme.granularity)) {
        case 'uniform': return createUniformColor(ctx, locationIt, colorTheme.color, colorData)
        case 'group': return createGroupColor(ctx, locationIt, colorTheme.color, colorData)
        case 'groupInstance': return createGroupInstanceColor(ctx, locationIt, colorTheme.color, colorData)
        case 'instance': return createInstanceColor(ctx, locationIt, colorTheme.color, colorData)
    }
}

export function createValueColor(value: Color, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.uColor, Color.toRgbNormalized(value) as Vec3)
        if (colorData.dColorType.ref.value !== 'uniform') {
            ValueCell.update(colorData.dColorType, 'uniform')
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Color.toRgbNormalized(value) as Vec3),
            aColor: ValueCell.create(new Float32Array(0)),
            tColor: ValueCell.create({ array: new Uint8Array(3), width: 1, height: 1 }),
            uColorTexDim: ValueCell.create(Vec2.create(1, 1)),
            dColorType: ValueCell.create('uniform'),
        }
    }
}

/** Creates color uniform */
export async function createUniformColor(ctx: RuntimeContext, locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): Promise<ColorData> {
    return createValueColor(color(NullLocation, false), colorData)
}

export function createTextureColor(colors: TextureImage<Uint8Array>, type: ColorType, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.tColor, colors)
        ValueCell.update(colorData.uColorTexDim, Vec2.create(colors.width, colors.height))
        if (colorData.dColorType.ref.value !== type) {
            ValueCell.update(colorData.dColorType, type)
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Vec3.zero()),
            aColor: ValueCell.create(new Float32Array(0)),
            tColor: ValueCell.create(colors),
            uColorTexDim: ValueCell.create(Vec2.create(colors.width, colors.height)),
            dColorType: ValueCell.create(type),
        }
    }
}

/** Creates color texture with color for each instance/unit */
export async function createInstanceColor(ctx: RuntimeContext, locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): Promise<ColorData> {
    const { instanceCount } = locationIt
    const colors = colorData && colorData.tColor.ref.value.array.length >= instanceCount * 3 ? colorData.tColor.ref.value : createTextureImage(instanceCount, 3)
    let i = 0
    while (locationIt.hasNext) {
        const { location, isSecondary, instanceIndex } = locationIt.move()
        Color.toArray(color(location, isSecondary), colors.array, instanceIndex * 3)
        locationIt.skipInstance()
        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Creating instance colors', current: i, max: instanceCount });
        }
        ++i
    }
    return createTextureColor(colors, 'instance', colorData)
}

/** Creates color texture with color for each group (i.e. shared across instances/units) */
export async function createGroupColor(ctx: RuntimeContext, locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): Promise<ColorData> {
    const { groupCount } = locationIt
    const colors = colorData && colorData.tColor.ref.value.array.length >= groupCount * 3 ? colorData.tColor.ref.value : createTextureImage(groupCount, 3)
    let i = 0
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const { location, isSecondary, groupIndex } = locationIt.move()
        Color.toArray(color(location, isSecondary), colors.array, groupIndex * 3)
        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Creating group colors', current: i, max: groupCount });
        }
        ++i
    }
    return createTextureColor(colors, 'group', colorData)
}

/** Creates color texture with color for each group in each instance (i.e. for each unit) */
export async function createGroupInstanceColor(ctx: RuntimeContext, locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): Promise<ColorData> {
    const { groupCount, instanceCount } = locationIt
    const count = instanceCount * groupCount
    const colors = colorData && colorData.tColor.ref.value.array.length >= count * 3 ? colorData.tColor.ref.value : createTextureImage(count, 3)
    let i = 0
    while (locationIt.hasNext) {
        const { location, isSecondary, index } = locationIt.move()
        Color.toArray(color(location, isSecondary), colors.array, index * 3)
        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Creating group instance colors', current: i, max: count });
        }
        ++i
    }
    return createTextureColor(colors, 'groupInstance', colorData)
}