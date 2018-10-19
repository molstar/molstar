/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#pragma glslify: encodeFloatRGBA = require(../utils/encode-float-rgba.glsl)

vec4 encodeIdRGBA(const in float v) {
	return encodeFloatRGBA(v + 1.0);
}

#pragma glslify: export(encodeIdRGBA)