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
    #if defined(dGridTexType_2d)
        precision mediump sampler2D;
        uniform sampler2D tMinDistanceTex;
        uniform vec2 uGridTexDim;
    #elif defined(dGridTexType_3d)
        precision highp sampler3D;
        uniform sampler3D tMinDistanceTex;
    #endif
    varying float vGroup;
#endif

#pragma glslify: encodeIdRGBA = require(./utils/encode-id-rgba.glsl)
#pragma glslify: texture3dFrom2dNearest = require(./utils/texture3d-from-2d-nearest.glsl)

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform float uCurrentSlice;
uniform float uCurrentX;
uniform float uCurrentY;
uniform float uAlpha;

#if defined(dCalcType_groupId)
    #if defined(dGridTexType_2d)
        vec4 textureMinDist(vec3 pos) {
            return texture3dFrom2dNearest(tMinDistanceTex, pos, uGridDim, uGridTexDim);
        }
    #elif defined(dGridTexType_3d)
        vec4 textureMinDist(vec3 pos) {
            return texture(tMinDistanceTex, pos);
        }
    #endif
#endif

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
        gl_FragColor.a = 1.0 - encodeDistLog(dist);
    #elif defined(dCalcType_groupId)
        float minDistance = decodeDistLog(1.0 - textureMinDist(fragPos).a);
        if (dist > minDistance + length(uBboxSize / uGridDim) / 1.5)
            discard;
        gl_FragColor = encodeIdRGBA(vGroup);
    #endif
}