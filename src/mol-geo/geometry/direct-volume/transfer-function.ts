/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TextureImage } from '../../../mol-gl/renderable/util';
import { spline } from '../../../mol-math/interpolate';
import { ColorScale, Color } from '../../../mol-util/color';
import { ValueCell } from '../../../mol-util';
import { Vec2 } from '../../../mol-math/linear-algebra';
import { ColorListName } from '../../../mol-util/color/lists';

export interface ControlPoint { x: number, alpha: number }

export function getControlPointsFromString(s: string): ControlPoint[] {
    return s.split(/\s*,\s*/).map(p => {
        const ps = p.split(/\s*:\s*/);
        return { x: parseFloat(ps[0]), alpha: parseFloat(ps[1]) };
    });
}

export function getControlPointsFromVec2Array(array: Vec2[]): ControlPoint[] {
    return array.map(v => ({ x: v[0], alpha: v[1] }));
}

export function createTransferFunctionTexture(controlPoints: ControlPoint[], listOrName: Color[] | ColorListName, texture?: ValueCell<TextureImage<Uint8Array>>): ValueCell<TextureImage<Uint8Array>> {
    const cp = [
        { x: 0, alpha: 0 },
        { x: 0, alpha: 0 },
        ...controlPoints,
        { x: 1, alpha: 0 },
        { x: 1, alpha: 0 },
    ];
    const scale = ColorScale.create({ domain: [0, 1], listOrName });

    const n = 256;
    const array = texture ? texture.ref.value.array : new Uint8Array(n * 4);

    let k = 0;
    let x1: number, x2: number;
    let a0: number, a1: number, a2: number, a3: number;

    const il = controlPoints.length + 1;
    for (let i = 0; i < il; ++i) {
        x1 = cp[i + 1].x;
        x2 = cp[i + 2].x;

        a0 = cp[i].alpha;
        a1 = cp[i + 1].alpha;
        a2 = cp[i + 2].alpha;
        a3 = cp[i + 3].alpha;

        const jl = Math.round((x2 - x1) * n);
        for (let j = 0; j < jl; ++j) {
            const t = j / jl;
            array[k * 4 + 3] = Math.max(0, spline(a0, a1, a2, a3, t, 0.5) * 255);
            scale.colorToArray(k / 255, array, k * 4);
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