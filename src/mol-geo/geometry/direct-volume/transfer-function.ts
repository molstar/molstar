/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TextureImage } from '../../../mol-gl/renderable/util';
import { spline } from '../../../mol-math/interpolate';
import { ValueCell } from '../../../mol-util';
import { Vec2 } from '../../../mol-math/linear-algebra';

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

const TF_WIDTH = 256;

/** Build the 1D alpha curve (length 256) from sorted control points using a Catmull-Rom spline. */
function build1DAlphaArray(controlPoints: ControlPoint[], out?: Uint8Array): Uint8Array {
    const cp = [
        { x: 0, alpha: 0 },
        { x: 0, alpha: 0 },
        ...controlPoints,
        { x: 1, alpha: 0 },
        { x: 1, alpha: 0 },
    ];

    const array = out ?? new Uint8Array(TF_WIDTH);
    if (out) array.fill(0);

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

        const jl = Math.round((x2 - x1) * TF_WIDTH);
        for (let j = 0; j < jl; ++j) {
            const t = j / jl;
            if (k < TF_WIDTH) {
                array[k] = Math.max(0, spline(a0, a1, a2, a3, t, 0.5) * 255);
                ++k;
            }
        }
    }
    return array;
}

export function createTransferFunctionTexture(controlPoints: ControlPoint[], texture?: ValueCell<TextureImage<Uint8Array>>): ValueCell<TextureImage<Uint8Array>> {
    // If updating an existing 2D texture, fall back to creating a fresh image so the size matches.
    const reuse = texture && texture.ref.value.height === 1 && texture.ref.value.array.length === TF_WIDTH;
    const array = reuse ? texture!.ref.value.array : new Uint8Array(TF_WIDTH);
    build1DAlphaArray(controlPoints, array);

    const textureImage = { array, width: TF_WIDTH, height: 1 };
    if (texture) {
        ValueCell.update(texture, textureImage);
        return texture;
    } else {
        return ValueCell.create(textureImage);
    }
}
