/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

float decodeFloatRGB(const in vec3 rgb) {
    return rgb.r * 256.0 * 256.0 * 255.0 + rgb.g * 256.0 * 255.0 + rgb.b * 255.0;
}

#pragma glslify: export(decodeFloatRGB)