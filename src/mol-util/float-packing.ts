/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clamp } from '../mol-math/interpolate';
import { Vec3, Vec4 } from '../mol-math/linear-algebra';
import { NumberArray } from './type-helpers';

/** encode positive integer as rgb byte triplet into array at offset */
export function packIntToRGBArray(value: number, array: NumberArray, offset: number) {
    value = clamp(Math.round(value), 0, 16777216 - 1) + 1;
    array[offset + 2] = value % 256;
    value = Math.floor(value / 256);
    array[offset + 1] = value % 256;
    value = Math.floor(value / 256);
    array[offset] = value % 256;
    return array;
}

/** decode positive integer encoded as rgb byte triplet */
export function unpackRGBToInt(r: number, g: number, b: number) {
    return (Math.floor(r) * 256 * 256 + Math.floor(g) * 256 + Math.floor(b)) - 1;
}

const UnpackDownscale = 255 / 256; // 0..1 -> fraction (excluding 1)
const PackFactors = Vec3.create(256 * 256 * 256, 256 * 256, 256);
const UnpackFactors = Vec4.create(
    UnpackDownscale / PackFactors[0],
    UnpackDownscale / PackFactors[1],
    UnpackDownscale / PackFactors[2],
    UnpackDownscale / 1
);

const tmpDepthRGBA = Vec4();
export function unpackRGBAToDepth(r: number, g: number, b: number, a: number) {
    Vec4.set(tmpDepthRGBA, r / 255, g / 255, b / 255, a / 255);
    return Vec4.dot(tmpDepthRGBA, UnpackFactors);
}