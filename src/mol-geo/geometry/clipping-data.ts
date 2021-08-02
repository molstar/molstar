/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

    uClipObjectType: ValueCell<number[]>,
    uClipObjectInvert: ValueCell<boolean[]>,
    uClipObjectPosition: ValueCell<number[]>,
    uClipObjectRotation: ValueCell<number[]>,
    uClipObjectScale: ValueCell<number[]>,
}

export function applyClippingGroups(array: Uint8Array, start: number, end: number, groups: Clipping.Groups) {
    array.fill(groups, start, end);
    return true;
}

export function clearClipping(array: Uint8Array, start: number, end: number) {
    array.fill(0, start, end);
}

export function createClipping(count: number, variant: Clipping.Variant, objects: Clipping.Objects, clippingData?: ClippingData): ClippingData {
    const clipping = createTextureImage(Math.max(1, count), 1, Uint8Array, clippingData && clippingData.tClipping.ref.value.array);
    if (clippingData) {
        ValueCell.update(clippingData.dClipObjectCount, objects.count);
        ValueCell.update(clippingData.dClipVariant, variant);

        ValueCell.update(clippingData.tClipping, clipping);
        ValueCell.update(clippingData.uClippingTexDim, Vec2.create(clipping.width, clipping.height));
        ValueCell.updateIfChanged(clippingData.dClipping, count > 0);

        ValueCell.update(clippingData.uClipObjectType, objects.type);
        ValueCell.update(clippingData.uClipObjectInvert, objects.invert);
        ValueCell.update(clippingData.uClipObjectPosition, objects.position);
        ValueCell.update(clippingData.uClipObjectRotation, objects.rotation);
        ValueCell.update(clippingData.uClipObjectScale, objects.scale);
        return clippingData;
    } else {
        return {
            dClipObjectCount: ValueCell.create(objects.count),
            dClipVariant: ValueCell.create(variant),

            tClipping: ValueCell.create(clipping),
            uClippingTexDim: ValueCell.create(Vec2.create(clipping.width, clipping.height)),
            dClipping: ValueCell.create(count > 0),

            uClipObjectType: ValueCell.create(objects.type),
            uClipObjectInvert: ValueCell.create(objects.invert),
            uClipObjectPosition: ValueCell.create(objects.position),
            uClipObjectRotation: ValueCell.create(objects.rotation),
            uClipObjectScale: ValueCell.create(objects.scale),
        };
    }
}

const emptyClippingTexture = { array: new Uint8Array(1), width: 1, height: 1 };
export function createEmptyClipping(clippingData?: ClippingData): ClippingData {
    if (clippingData) {
        ValueCell.update(clippingData.dClipObjectCount, 0);
        ValueCell.update(clippingData.tClipping, emptyClippingTexture);
        ValueCell.update(clippingData.uClippingTexDim, Vec2.create(1, 1));
        ValueCell.updateIfChanged(clippingData.dClipping, false);
        return clippingData;
    } else {
        const { objects, variant } = Clipping.Empty;
        return {
            dClipObjectCount: ValueCell.create(0),
            dClipVariant: ValueCell.create(variant),

            tClipping: ValueCell.create(emptyClippingTexture),
            uClippingTexDim: ValueCell.create(Vec2.create(1, 1)),
            dClipping: ValueCell.create(false),

            uClipObjectType: ValueCell.create(objects.type),
            uClipObjectInvert: ValueCell.create(objects.invert),
            uClipObjectPosition: ValueCell.create(objects.position),
            uClipObjectRotation: ValueCell.create(objects.rotation),
            uClipObjectScale: ValueCell.create(objects.scale),
        };
    }
}