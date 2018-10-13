/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision mediump float;

varying vec3 position;
varying float radius;

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform float uCurrentSlice;
uniform float uCurrentX;
uniform float uCurrentY;
uniform float uAlpha;

// #if dDrawBuffers == 8
    // layout(location = 0) out vec4 out0;
    layout(location = 1) out vec4 out1;
    layout(location = 2) out vec4 out2;
    layout(location = 3) out vec4 out3;
    layout(location = 4) out vec4 out4;
    layout(location = 5) out vec4 out5;
    layout(location = 6) out vec4 out6;
    layout(location = 7) out vec4 out7;
// #endif

float calcDensity(float x, float y, float z) {
    vec3 fragPos = vec3(x, y, z) / uGridDim;
    float dist = length(fragPos * uBboxSize - position * uBboxSize);
    float density = exp(-uAlpha * ((dist * dist) / (radius * radius)));
    return density;
}

void main() {
    vec2 tmpVec = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    out_FragColor = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 0.0));
    out1 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 1.0));
    out2 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 2.0));
    out3 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 3.0));
    out4 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 4.0));
    out5 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 5.0));
    out6 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 6.0));
    out7 = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice + 7.0));
    // gl_FragColor = vec4(1, 1, 1, calcDensity(tmpVec.x, tmpVec.y, uCurrentSlice));
}