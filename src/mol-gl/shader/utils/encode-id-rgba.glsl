#pragma glslify: encodeFloatRGBA = require(../utils/encode-float-rgba.glsl)

vec4 encodeIdRGBA(const in float v) {
	return encodeFloatRGBA(v + 1.0);
}

#pragma glslify: export(encodeIdRGBA)