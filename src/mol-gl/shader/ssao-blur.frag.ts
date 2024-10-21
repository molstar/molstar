/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

export const ssaoBlur_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tSsaoDepth;
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform float uKernel[dOcclusionKernelSize];
uniform float uBlurDepthBias;
uniform float uBlurNormalBias;
uniform float uBlurStepSize;
uniform vec2 uBlurStepOffset;

uniform float uBlurDirectionX;
uniform float uBlurDirectionY;

uniform mat4 uInvProjection;
uniform float uNear;
uniform float uFar;

#include common

float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

bool isBackground(const in float depth) {
    return depth > 0.999;
}

bool isNearClip(const in float depth) {
    return depth == 0.0;
}

bool outsideBounds(const in vec2 p) {
    return p.x < uBounds.x || p.y < uBounds.y || p.x > uBounds.z || p.y > uBounds.w;
}

float getPixelSize(const in vec2 coords, const in float depth) {
    vec3 viewPos0 = screenSpaceToViewSpace(vec3(coords, depth), uInvProjection);
    vec3 viewPos1 = screenSpaceToViewSpace(vec3(coords + vec2(1.0, 0.0) / uTexSize, depth), uInvProjection);
    return distance(viewPos0, viewPos1);
}

#ifdef dBlurNormalBias
    float getDepth(const in vec2 coords) {
        return unpackRGToUnitInterval(texture2D(tSsaoDepth, coords).zw);
    }

    // adapted from https://gist.github.com/bgolus/a07ed65602c009d5e2f753826e8078a0
    vec3 viewNormalAtPixelPositionAccurate(const in vec2 vpos) {
        // current pixel's depth
        float c = getDepth(vpos);

        // get current pixel's view space position
        vec3 viewSpacePos_c = screenSpaceToViewSpace(vec3(vpos, c), uInvProjection);

        // get view space position at 1 pixel offsets in each major direction
        vec3 viewSpacePos_l = screenSpaceToViewSpace(vec3(vpos + vec2(-1.0, 0.0) / uTexSize, getDepth(vpos + vec2(-1.0, 0.0) / uTexSize)), uInvProjection);
        vec3 viewSpacePos_r = screenSpaceToViewSpace(vec3(vpos + vec2( 1.0, 0.0) / uTexSize, getDepth(vpos + vec2( 1.0, 0.0) / uTexSize)), uInvProjection);
        vec3 viewSpacePos_d = screenSpaceToViewSpace(vec3(vpos + vec2( 0.0,-1.0) / uTexSize, getDepth(vpos + vec2( 0.0,-1.0) / uTexSize)), uInvProjection);
        vec3 viewSpacePos_u = screenSpaceToViewSpace(vec3(vpos + vec2( 0.0, 1.0) / uTexSize, getDepth(vpos + vec2( 0.0, 1.0) / uTexSize)), uInvProjection);

        // get the difference between the current and each offset position
        vec3 l = viewSpacePos_c - viewSpacePos_l;
        vec3 r = viewSpacePos_r - viewSpacePos_c;
        vec3 d = viewSpacePos_c - viewSpacePos_d;
        vec3 u = viewSpacePos_u - viewSpacePos_c;

        // get depth values at 1 & 2 pixels offsets from current along the horizontal axis
        vec4 H = vec4(
            getDepth(vpos + vec2(-1.0, 0.0) / uTexSize),
            getDepth(vpos + vec2( 1.0, 0.0) / uTexSize),
            getDepth(vpos + vec2(-2.0, 0.0) / uTexSize),
            getDepth(vpos + vec2( 2.0, 0.0) / uTexSize)
        );

        // get depth values at 1 & 2 pixels offsets from current along the vertical axis
        vec4 V = vec4(
            getDepth(vpos + vec2(0.0,-1.0) / uTexSize),
            getDepth(vpos + vec2(0.0, 1.0) / uTexSize),
            getDepth(vpos + vec2(0.0,-2.0) / uTexSize),
            getDepth(vpos + vec2(0.0, 2.0) / uTexSize)
        );

        // current pixel's depth difference from slope of offset depth samples
        // differs from original article because we're using non-linear depth values
        // see article's comments
        vec2 he = abs((2.0 * H.xy - H.zw) - c);
        vec2 ve = abs((2.0 * V.xy - V.zw) - c);

        // pick horizontal and vertical diff with the smallest depth difference from slopes
        vec3 hDeriv = he.x < he.y ? l : r;
        vec3 vDeriv = ve.x < ve.y ? d : u;

        // get view space normal from the cross product of the best derivatives
        vec3 viewNormal = normalize(cross(hDeriv, vDeriv));

        return viewNormal;
    }
#endif

void main(void) {
    vec2 coords = gl_FragCoord.xy / uTexSize;

    vec2 packedDepth = texture2D(tSsaoDepth, coords).zw;

    if (outsideBounds(coords)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packedDepth);
        return;
    }

    float selfDepth = unpackRGToUnitInterval(packedDepth);
    if (isBackground(selfDepth) || isNearClip(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packedDepth);
        return;
    }

    float selfViewZ = getViewZ(selfDepth);
    float pixelSize = getPixelSize(coords, selfDepth);

    #ifdef dBlurNormalBias
        vec3 selfNormal = viewNormalAtPixelPositionAccurate(coords);
    #endif

    vec2 offset = vec2(uBlurDirectionX, uBlurDirectionY) / uTexSize;
    coords += uBlurStepOffset * vec2(uBlurDirectionX, uBlurDirectionY);
    offset *= uBlurStepSize;

    float sum = 0.0;
    float kernelSum = 0.0;
    // only if kernelSize is odd
    for (int i = -dOcclusionKernelSize / 2; i <= dOcclusionKernelSize / 2; i++) {
        if (abs(float(i)) > 1.0 && abs(float(i)) * pixelSize > 0.8) continue;

        vec2 sampleCoords = coords + float(i) * offset;
        if (outsideBounds(sampleCoords)) {
            continue;
        }

        vec4 sampleSsaoDepth = texture2D(tSsaoDepth, sampleCoords);

        float sampleDepth = unpackRGToUnitInterval(sampleSsaoDepth.zw);
        if (isBackground(sampleDepth) || isNearClip(sampleDepth)) {
            continue;
        }

        #ifdef dBlurNormalBias
            vec3 sampleNormal = viewNormalAtPixelPositionAccurate(sampleCoords);
            if (saturate(dot(selfNormal, sampleNormal)) < uBlurNormalBias) {
                continue;
            }
        #endif

        float sampleViewZ = getViewZ(sampleDepth);
        if (abs(selfViewZ - sampleViewZ) >= uBlurDepthBias) {
            continue;
        }

        float kernel = uKernel[int(abs(float(i)))]; // abs is not defined for int in webgl1
        float sampleValue = unpackRGToUnitInterval(sampleSsaoDepth.xy);

        sum += kernel * sampleValue;
        kernelSum += kernel;
    }
    gl_FragColor = vec4(packUnitIntervalToRG(sum / kernelSum), packedDepth);
}
`;