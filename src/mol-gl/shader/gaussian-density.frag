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

#if dDrawBuffers >= 4
    layout(location = 1) out vec4 out1;
    layout(location = 2) out vec4 out2;
    layout(location = 3) out vec4 out3;
#endif
#if dDrawBuffers >= 8
    layout(location = 4) out vec4 out4;
    layout(location = 5) out vec4 out5;
    layout(location = 6) out vec4 out6;
    layout(location = 7) out vec4 out7;
#endif

float calcDensity(float x, float y, float z) {
    vec3 fragPos = vec3(x, y, z) / uGridDim;
    float dist = length(fragPos * uBboxSize - position * uBboxSize);
    float density = exp(-uAlpha * ((dist * dist) / (radius * radius)));
    return density;
}

const vec3 color = vec3(1.0, 1.0, 1.0);

void main() {
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    gl_FragColor = vec4(color, calcDensity(v.x, v.y, uCurrentSlice));
    #if dDrawBuffers >= 4
        out1 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 1.0));
        out2 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 2.0));
        out3 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 3.0));
    #endif
    #if dDrawBuffers >= 8
        out4 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 4.0));
        out5 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 5.0));
        out6 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 6.0));
        out7 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 7.0));
    #endif
}