/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { Color } from '../../mol-util/color';

export type OverpaintData = {
    tOverpaint: ValueCell<TextureImage<Uint8Array>>
    uOverpaintTexDim: ValueCell<Vec2>
    dOverpaint: ValueCell<boolean>,
}

export function applyOverpaintColor(array: Uint8Array, start: number, end: number, color: Color) {
    for (let i = start; i < end; ++i) {
        Color.toArray(color, array, i * 4);
        array[i * 4 + 3] = 255;
    }
    return true;
}

export function clearOverpaint(array: Uint8Array, start: number, end: number) {
    array.fill(0, start * 4, end * 4);
    return true;
}

export function createOverpaint(count: number, overpaintData?: OverpaintData): OverpaintData {
    const overpaint = createTextureImage(Math.max(1, count), 4, Uint8Array, overpaintData && overpaintData.tOverpaint.ref.value.array);
    if (overpaintData) {
        ValueCell.update(overpaintData.tOverpaint, overpaint);
        ValueCell.update(overpaintData.uOverpaintTexDim, Vec2.create(overpaint.width, overpaint.height));
        ValueCell.update(overpaintData.dOverpaint, count > 0);
        return overpaintData;
    } else {
        return {
            tOverpaint: ValueCell.create(overpaint),
            uOverpaintTexDim: ValueCell.create(Vec2.create(overpaint.width, overpaint.height)),
            dOverpaint: ValueCell.create(count > 0),
        };
    }
}

const emptyOverpaintTexture = { array: new Uint8Array(4), width: 1, height: 1 };
export function createEmptyOverpaint(overpaintData?: OverpaintData): OverpaintData {
    if (overpaintData) {
        ValueCell.update(overpaintData.tOverpaint, emptyOverpaintTexture);
        ValueCell.update(overpaintData.uOverpaintTexDim, Vec2.create(1, 1));
        return overpaintData;
    } else {
        return {
            tOverpaint: ValueCell.create(emptyOverpaintTexture),
            uOverpaintTexDim: ValueCell.create(Vec2.create(1, 1)),
            dOverpaint: ValueCell.create(false),
        };
    }
}