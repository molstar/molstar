/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';

export type WiggleType = 'instance' | 'groupInstance';

export type WiggleData = {
    tWiggle: ValueCell<TextureImage<Uint8Array>>
    uWiggleTexDim: ValueCell<Vec2>
    dWiggle: ValueCell<boolean>,
    wiggleAverage: ValueCell<number>,
    dWiggleType: ValueCell<string>,
    uWiggleStrength: ValueCell<number>,
}

export function applyWiggleValue(array: Uint8Array, start: number, end: number, value: number) {
    for (let i = start; i < end; ++i) {
        array[i] = value * 255;
    }
    return true;
}

export function getWiggleAverage(array: Uint8Array, count: number): number {
    if (count === 0 || array.length < count) return 0;
    let sum = 0;
    for (let i = 0; i < count; ++i) {
        sum += array[i];
    }
    return sum / (255 * count);
}

export function clearWiggle(array: Uint8Array, start: number, end: number) {
    array.fill(0, start, end);
}

export function createWiggle(count: number, type: WiggleType, wiggleData?: WiggleData): WiggleData {
    const wiggle = createTextureImage(Math.max(1, count), 1, Uint8Array, wiggleData && wiggleData.tWiggle.ref.value.array);
    if (wiggleData) {
        ValueCell.update(wiggleData.tWiggle, wiggle);
        ValueCell.update(wiggleData.uWiggleTexDim, Vec2.create(wiggle.width, wiggle.height));
        ValueCell.updateIfChanged(wiggleData.dWiggle, count > 0);
        ValueCell.updateIfChanged(wiggleData.wiggleAverage, getWiggleAverage(wiggle.array, count));
        ValueCell.updateIfChanged(wiggleData.dWiggleType, type);
        return wiggleData;
    } else {
        return {
            tWiggle: ValueCell.create(wiggle),
            uWiggleTexDim: ValueCell.create(Vec2.create(wiggle.width, wiggle.height)),
            dWiggle: ValueCell.create(count > 0),
            wiggleAverage: ValueCell.create(0),
            dWiggleType: ValueCell.create(type),
            uWiggleStrength: ValueCell.create(1),
        };
    }
}

const emptyWiggleTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyWiggle(wiggleData?: WiggleData): WiggleData {
    if (wiggleData) {
        ValueCell.update(wiggleData.tWiggle, emptyWiggleTexture);
        ValueCell.update(wiggleData.uWiggleTexDim, Vec2.create(1, 1));
        return wiggleData;
    } else {
        return {
            tWiggle: ValueCell.create(emptyWiggleTexture),
            uWiggleTexDim: ValueCell.create(Vec2.create(1, 1)),
            dWiggle: ValueCell.create(false),
            wiggleAverage: ValueCell.create(0),
            dWiggleType: ValueCell.create('groupInstance'),
            uWiggleStrength: ValueCell.create(1),
        };
    }
}
