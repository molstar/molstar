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

void main() {
    vec3 tmpVec = gl_FragCoord.xyz;
    tmpVec.x = tmpVec.x - uCurrentX;
    tmpVec.y = tmpVec.y - uCurrentY;
    vec3 fragPos = vec3(
        (tmpVec.x - 0.5) / uGridDim.x,
        (tmpVec.y - 0.5) / uGridDim.y,
        (uCurrentSlice) / uGridDim.z
    );
    float dist = length(fragPos * uBboxSize - position * uBboxSize);
    float density = 1.0 - smoothstep( 0.0, radius * 2.0, dist);
    gl_FragColor = vec4(1, 1, 1, density);
    // density = 1.0 - clamp((dist - (radius + 1.4)) + 0.5, 0.0, 1.0);				
    // gl_FragColor = vec4(vec3(density), 1.0);
}