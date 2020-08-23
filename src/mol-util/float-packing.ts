/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clamp } from '../mol-math/interpolate';
import { fasterExp, fasterLog } from '../mol-math/approx';

const maxFloat = 10000.0; // NOTE same constant is set in shaders
const floatLogFactor = fasterLog(maxFloat + 1.0);

/** encode float logarithmically */
export function encodeFloatLog(value: number) { return fasterLog(value + 1.0) / floatLogFactor; }

/** decode logarithmically encoded float */
export function decodeFloatLog(value: number) { return fasterExp(value * floatLogFactor) - 1.0; }

/** encode float as normalized rgb triplet */
export function encodeFloatRGB(value: number) {
    value = clamp(value, 0.0, 16777216.0 - 1.0) + 1.0;
    const b = (value % 256) / 255.0;
    value = Math.floor(value / 256.0);
    const g = (value % 256) / 255.0;
    value = Math.floor(value / 256.0);
    const r = (value % 256) / 255.0;
    return [r, g, b];
}

/** decode float encoded as rgb triplet */
export function decodeFloatRGB(r: number, g: number, b: number) {
    return (Math.floor(r) * 256 * 256 + Math.floor(g) * 256 + Math.floor(b)) - 1;
}