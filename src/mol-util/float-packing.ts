/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clamp } from '../mol-math/interpolate';
import { fasterExp, fasterLog } from '../mol-math/approx';
import { Vec3, Vec4 } from '../mol-math/linear-algebra';
import { NumberArray } from './type-helpers';

const maxFloat = 10000.0; // NOTE same constant is set in shaders
const floatLogFactor = fasterLog(maxFloat + 1);

/** encode float logarithmically */
export function encodeFloatLog(value: number) { return fasterLog(value + 1) / floatLogFactor; }

/** decode logarithmically encoded float */
export function decodeFloatLog(value: number) { return fasterExp(value * floatLogFactor) - 1; }

/** encode float as rgb triplet into array at offset */
export function encodeFloatRGBtoArray(value: number, array: NumberArray, offset: number) {
    value = clamp(value, 0, 16777216 - 1) + 1;
    array[offset + 2] = value % 256;
    value = Math.floor(value / 256);
    array[offset + 1] = value % 256;
    value = Math.floor(value / 256);
    array[offset] = value % 256;
    return array;
}

/** decode float encoded as rgb triplet */
export function decodeFloatRGB(r: number, g: number, b: number) {
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