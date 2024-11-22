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
    // checking for 1.0 is not enough, because of precision issues
    return depth >= 0.999;
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

void main(void) {
    vec2 coords = gl_FragCoord.xy / uTexSize;

    vec2 packedDepth = texture2D(tSsaoDepth, coords).zw;

    if (outsideBounds(coords)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packedDepth);
        return;
    }

    float selfDepth = unpackRGToUnitInterval(packedDepth);
    // (if background and if second pass) or if near clip
    if ((isBackground(selfDepth) && uBlurDirectionY != 0.0) || isNearClip(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packedDepth);
        return;
    }

    float selfViewZ = getViewZ(selfDepth);
    float pixelSize = getPixelSize(coords, selfDepth);

    vec2 offset = vec2(uBlurDirectionX, uBlurDirectionY) / uTexSize;

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