/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision highp float;

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

float calcDensity(float x, float y, float z, float radiusSq) {
    vec3 fragPos = vec3(x, y, z) / uGridDim;
    float dist = distance(fragPos * uBboxSize, position * uBboxSize);
    float density = exp(-uAlpha * ((dist * dist) / radiusSq));
    return density;
}

const vec3 color = vec3(1.0, 1.0, 1.0);

void main() {
    float radiusSq = radius * radius;
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    out_FragColor = vec4(color, calcDensity(v.x, v.y, uCurrentSlice, radiusSq));
    #if dDrawBuffers >= 4
        out1 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 1.0, radiusSq));
        out2 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 2.0, radiusSq));
        out3 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 3.0, radiusSq));
    #endif
    #if dDrawBuffers >= 8
        out4 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 4.0, radiusSq));
        out5 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 5.0, radiusSq));
        out6 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 6.0, radiusSq));
        out7 = vec4(color, calcDensity(v.x, v.y, uCurrentSlice + 7.0, radiusSq));
    #endif
}