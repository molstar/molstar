/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#pragma glslify: encodeFloatRGB = require(../utils/encode-float-rgb.glsl)

vec3 encodeIdRGB(const in float v) {
	return encodeFloatRGB(v + 1.0);
}

#pragma glslify: export(encodeIdRGB)