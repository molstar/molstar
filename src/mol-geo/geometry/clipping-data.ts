/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { Clipping } from '../../mol-theme/clipping';

export type ClippingData = {
    dClipObjectCount: ValueCell<number>,
    dClipVariant: ValueCell<string>,

    tClipping: ValueCell<TextureImage<Uint8Array>>
    uClippingTexDim: ValueCell<Vec2>
    dClipping: ValueCell<boolean>,
}

export function applyClippingGroups(array: Uint8Array, start: number, end: number, groups: Clipping.Groups) {
    for (let i = start; i < end; ++i) {
        array[i] = groups;
    }
    return true;
}

export function clearClipping(array: Uint8Array, start: number, end: number) {
    array.fill(0, start, end);
}

export function createClipping(count: number, clippingData?: ClippingData): ClippingData {
    const clipping = createTextureImage(Math.max(1, count), 1, Uint8Array, clippingData && clippingData.tClipping.ref.value.array);
    if (clippingData) {
        ValueCell.update(clippingData.tClipping, clipping);
        ValueCell.update(clippingData.uClippingTexDim, Vec2.create(clipping.width, clipping.height));
        ValueCell.update(clippingData.dClipping, count > 0);
        return clippingData;
    } else {
        return {
            dClipObjectCount: ValueCell.create(0),
            dClipVariant: ValueCell.create('instance'),

            tClipping: ValueCell.create(clipping),
            uClippingTexDim: ValueCell.create(Vec2.create(clipping.width, clipping.height)),
            dClipping: ValueCell.create(count > 0),
        };
    }
}

const emptyClippingTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyClipping(clippingData?: ClippingData): ClippingData {
    if (clippingData) {
        ValueCell.update(clippingData.tClipping, emptyClippingTexture);
        ValueCell.update(clippingData.uClippingTexDim, Vec2.create(1, 1));
        return clippingData;
    } else {
        return {
            dClipObjectCount: ValueCell.create(0),
            dClipVariant: ValueCell.create('instance'),

            tClipping: ValueCell.create(emptyClippingTexture),
            uClippingTexDim: ValueCell.create(Vec2.create(1, 1)),
            dClipping: ValueCell.create(false),
        };
    }
}