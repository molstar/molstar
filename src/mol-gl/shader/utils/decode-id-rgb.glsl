/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#pragma glslify: decodeFloatRGB = require(../utils/decode-float-rgb.glsl)

float decodeIdRGB(const in vec3 v) {
	return decodeFloatRGB(v) - 1.0;
}

#pragma glslify: export(decodeIdRGB)