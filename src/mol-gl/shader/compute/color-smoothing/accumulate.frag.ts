/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const accumulate_frag = `
precision highp float;

varying vec3 vPosition;
varying vec4 vColor;

uniform float uCurrentSlice;
uniform float uCurrentX;
uniform float uCurrentY;
uniform float uResolution;

const float p = 2.0;

void main() {
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    vec3 fragPos = vec3(v.x, v.y, uCurrentSlice);
    float dist = distance(fragPos, vPosition);
    if (dist > p) discard;

    float f = p - dist;
    gl_FragColor = vColor * f;
    gl_FragData[1] = vec4(f);
}
`;