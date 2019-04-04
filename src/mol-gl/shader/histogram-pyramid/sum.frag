/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp sampler2D;

uniform sampler2D tTexture;

// uv coordinate of each fragment obtained from vertex shader.
varying vec2 vCoordinate;

#pragma glslify: encodeFloatRGB = require(../utils/encode-float-rgb.glsl)

void main(void) {
    gl_FragColor = vec4(encodeFloatRGB(texture2D(tTexture, vCoordinate).r), 1.0);
}