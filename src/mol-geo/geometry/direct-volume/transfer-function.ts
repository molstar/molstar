/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TextureImage } from '../../../mol-gl/renderable/util';
import { spline } from '../../../mol-math/interpolate';
import { UUID, ValueCell } from '../../../mol-util';
import { ControlPoint } from '../../../mol-plugin-ui/controls/line-graph/line-graph-component';
import { Color } from '../../../mol-util/color';

export function createTransferFunctionTexture(controlPoints: ControlPoint[], texture?: ValueCell<TextureImage<Uint8Array>>): ValueCell<TextureImage<Uint8Array>> {
    const blackColor = Color.fromHexStyle('#000000');
    const cpStart: ControlPoint = {
        data: {
            x: 0, alpha: 0
        },
        color: blackColor,
        id: UUID.create22(),
        // could have issues with index, for now 0 everywhere
        index: 0,
    };
    const cpEnd: ControlPoint = {
        data: {
            x: 1, alpha: 0
        },
        color: blackColor,
        id: UUID.create22(),
        // could have issues with index, for now 0 everywhere
        index: 0,
    };
    const cp = [
        cpStart,
        cpStart,
        ...controlPoints,
        cpEnd,
        cpEnd
    ];

    const n = 256;
    const array = texture ? texture.ref.value.array : new Uint8Array(n);
    let k = 0;
    let x1: number, x2: number;
    let a0: number, a1: number, a2: number, a3: number;
    const il = controlPoints.length + 1;
    for (let i = 0; i < il; ++i) {

        x1 = cp[i + 1].data.x;
        x2 = cp[i + 2].data.x;

        a0 = cp[i].data.alpha;
        a1 = cp[i + 1].data.alpha;
        a2 = cp[i + 2].data.alpha;
        a3 = cp[i + 3].data.alpha;

        const jl = Math.round((x2 - x1) * n);
        for (let j = 0; j < jl; ++j) {
            const t = j / jl;
            array[k] = Math.max(0, spline(a0, a1, a2, a3, t, 0.5) * 255);
            ++k;
        }
    }
    const textureImage = { array, width: 256, height: 1 };
    if (texture) {
        ValueCell.update(texture, textureImage);
        return texture;
    } else {
        return ValueCell.create(textureImage);
    }
}