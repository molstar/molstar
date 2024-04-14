/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { createNullTexture, Texture } from '../../mol-gl/webgl/texture';

export type EmissiveType = 'instance' | 'groupInstance' | 'volumeInstance';

export type EmissiveData = {
    tEmissive: ValueCell<TextureImage<Uint8Array>>
    uEmissiveTexDim: ValueCell<Vec2>
    dEmissive: ValueCell<boolean>,
    emissiveAverage: ValueCell<number>,

    tEmissiveGrid: ValueCell<Texture>,
    uEmissiveGridDim: ValueCell<Vec3>,
    uEmissiveGridTransform: ValueCell<Vec4>,
    dEmissiveType: ValueCell<string>,
    uEmissiveStrength: ValueCell<number>,
}

export function applyEmissiveValue(array: Uint8Array, start: number, end: number, value: number) {
    for (let i = start; i < end; ++i) {
        array[i] = value * 255;
    }
    return true;
}

export function getEmissiveAverage(array: Uint8Array, count: number): number {
    if (count === 0 || array.length < count) return 0;
    let sum = 0;
    for (let i = 0; i < count; ++i) {
        sum += array[i];
    }
    return sum / (255 * count);
}

export function clearEmissive(array: Uint8Array, start: number, end: number) {
    array.fill(0, start, end);
}

export function createEmissive(count: number, type: EmissiveType, emissiveData?: EmissiveData): EmissiveData {
    const emissive = createTextureImage(Math.max(1, count), 1, Uint8Array, emissiveData && emissiveData.tEmissive.ref.value.array);
    if (emissiveData) {
        ValueCell.update(emissiveData.tEmissive, emissive);
        ValueCell.update(emissiveData.uEmissiveTexDim, Vec2.create(emissive.width, emissive.height));
        ValueCell.updateIfChanged(emissiveData.dEmissive, count > 0);
        ValueCell.updateIfChanged(emissiveData.emissiveAverage, getEmissiveAverage(emissive.array, count));
        ValueCell.updateIfChanged(emissiveData.dEmissiveType, type);
        return emissiveData;
    } else {
        return {
            tEmissive: ValueCell.create(emissive),
            uEmissiveTexDim: ValueCell.create(Vec2.create(emissive.width, emissive.height)),
            dEmissive: ValueCell.create(count > 0),
            emissiveAverage: ValueCell.create(0),

            tEmissiveGrid: ValueCell.create(createNullTexture()),
            uEmissiveGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uEmissiveGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dEmissiveType: ValueCell.create(type),
            uEmissiveStrength: ValueCell.create(1),
        };
    }
}

const emptyEmissiveTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyEmissive(emissiveData?: EmissiveData): EmissiveData {
    if (emissiveData) {
        ValueCell.update(emissiveData.tEmissive, emptyEmissiveTexture);
        ValueCell.update(emissiveData.uEmissiveTexDim, Vec2.create(1, 1));
        return emissiveData;
    } else {
        return {
            tEmissive: ValueCell.create(emptyEmissiveTexture),
            uEmissiveTexDim: ValueCell.create(Vec2.create(1, 1)),
            dEmissive: ValueCell.create(false),
            emissiveAverage: ValueCell.create(0),

            tEmissiveGrid: ValueCell.create(createNullTexture()),
            uEmissiveGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uEmissiveGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dEmissiveType: ValueCell.create('groupInstance'),
            uEmissiveStrength: ValueCell.create(1),
        };
    }
}