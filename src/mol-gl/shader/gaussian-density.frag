/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision highp float;

varying vec3 vPosition;
varying float vRadius;
#if defined(dCalcType_groupId)
    precision highp sampler3D;
    uniform sampler3D tMinDistanceTex;
    varying float vGroup;
#endif

#pragma glslify: encodeIdRGBA = require(./utils/encode-id-rgba.glsl)

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform float uCurrentSlice;
uniform float uCurrentX;
uniform float uCurrentY;
uniform float uAlpha;

// encode distance logarithmically with given maxDistance
const float maxDistance = 10000.0;
const float distLogFactor = log(maxDistance + 1.0);
float encodeDistLog(float dist) { return log(dist + 1.0) / distLogFactor; }
float decodeDistLog(float logDist) { return exp(logDist * distLogFactor) - 1.0; }

void main() {
    vec2 v = gl_FragCoord.xy - vec2(uCurrentX, uCurrentY) - 0.5;
    vec3 fragPos = vec3(v.x, v.y, uCurrentSlice) / uGridDim;
    float dist = distance(fragPos * uBboxSize, vPosition * uBboxSize);

    #if defined(dCalcType_density)
        float radiusSq = vRadius * vRadius;
        float density = exp(-uAlpha * ((dist * dist) / radiusSq));
        gl_FragColor = vec4(density);
    #elif defined(dCalcType_minDistance)
        gl_FragColor.r = 1.0 - encodeDistLog(dist);
    #elif defined(dCalcType_groupId)
        float minDistance = decodeDistLog(1.0 - texture(tMinDistanceTex, fragPos).r);
        if (dist > minDistance + log(minDistance) / 2.0)
            discard;
        gl_FragColor = encodeIdRGBA(vGroup);
    #endif
}