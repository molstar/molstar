/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export const gaussianDensity_frag = `
precision highp float;

varying vec3 vPosition;
varying float vRadiusSqInv;
#if defined(dCalcType_groupId)
    #if defined(dGridTexType_2d)
        precision highp sampler2D;
        uniform sampler2D tMinDistanceTex;
        uniform vec3 uGridTexDim;
    #elif defined(dGridTexType_3d)
        precision highp sampler3D;
        uniform sampler3D tMinDistanceTex;
    #endif
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
uniform float uRadiusFactorInv;

void main() {
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    vec3 fragPos = vec3(v.x, v.y, uCurrentSlice);
    float dist = distance(fragPos, vPosition) * uResolution;

    #if defined(dCalcType_density)
        float density = exp(-uAlpha * ((dist * dist) * vRadiusSqInv));
        gl_FragColor.a = density * uRadiusFactorInv;
    #elif defined(dCalcType_minDistance)
        gl_FragColor.a = 1.0 - dist * uRadiusFactorInv;
    #elif defined(dCalcType_groupId)
        #if defined(dGridTexType_2d)
            float minDistance = 1.0 - texture2D(tMinDistanceTex, (gl_FragCoord.xy) / (uGridTexDim.xy / uGridTexScale)).a;
        #elif defined(dGridTexType_3d)
            float minDistance = 1.0 - texelFetch(tMinDistanceTex, ivec3(gl_FragCoord.xy, uCurrentSlice), 0).a;
        #endif
        if (dist * uRadiusFactorInv > minDistance + uResolution * 0.05)
            discard;
        gl_FragColor.rgb = packIntToRGB(vGroup);
    #endif
}
`;