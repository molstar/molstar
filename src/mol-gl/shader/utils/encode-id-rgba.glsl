#pragma glslify: encodeFloatRGBA = require(../utils/encode-float-rgba.glsl)

#define MAX_ID 16777216.0

vec4 encodeIdRGBA( const in float v ) {
	return encodeFloatRGBA(((v + 1.0) / MAX_ID));
}

#pragma glslify: export(encodeIdRGBA)