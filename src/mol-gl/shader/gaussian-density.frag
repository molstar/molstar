/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision highp float;

varying vec3 vPosition;
varying float vRadius;

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform float uCurrentSlice;
uniform float uCurrentX;
uniform float uCurrentY;
uniform float uAlpha;

vec4 calc(float x, float y, float z, float radiusSq) {
    vec3 fragPos = vec3(x, y, z) / uGridDim;
    float dist = distance(fragPos * uBboxSize, vPosition * uBboxSize);
    float density = exp(-uAlpha * ((dist * dist) / radiusSq));
    return vec4(density);
}

void main() {
    float radiusSq = vRadius * vRadius;
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    gl_FragColor = calc(v.x, v.y, uCurrentSlice, radiusSq);
}