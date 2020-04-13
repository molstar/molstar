/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export default `
precision highp float;

varying vec3 vPosition;
varying float vRadiusSqInv;
#if defined(dCalcType_groupId)
    precision highp sampler2D;
    uniform sampler2D tMinDistanceTex;
    uniform vec3 uGridTexDim;
    varying float vGroup;
#endif

#include common

uniform vec3 uGridDim;
uniform vec2 uGridTexScale;
uniform float uCurrentSlice;
uniform float uCurrentX;
uniform float uCurrentY;
uniform float uAlpha;
uniform float uResolution;

void main() {
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    vec3 fragPos = vec3(v.x, v.y, uCurrentSlice);
    float dist = distance(fragPos, vPosition) * uResolution;

    #if defined(dCalcType_density)
        float density = exp(-uAlpha * ((dist * dist) * vRadiusSqInv));
        gl_FragColor.a = density;
    #elif defined(dCalcType_minDistance)
        gl_FragColor.a = 10000.0 - dist;
    #elif defined(dCalcType_groupId)
        float minDistance = 10000.0 - texture2D(tMinDistanceTex, (gl_FragCoord.xy) / (uGridTexDim.xy / uGridTexScale)).a;
        if (dist > minDistance + uResolution * 0.05)
            discard;
        gl_FragColor.rgb = encodeFloatRGB(vGroup);
    #endif
}
`;