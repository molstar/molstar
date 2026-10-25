/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { common } from './chunks/common.glsl.js';

export const ssaoBlur_frag = `
precision highp float;
precision highp int;

${common}

// the packed depth is linear and reversed (see the SSAO pass)
float getViewZ(const in float depth) {
    return depth * (uFar - uNear) - uFar;
}

// the background (far plane) at 0, packed as (0, 0); the near plane at 1, which
// packUnitIntervalToRG stores as (0, 1) in the uint8 target and unpacks to 255/256, never 1.0, so the packed bytes
// are tested
bool isBackgroundPacked(const in vec2 packedDepth) {
    return packedDepth.x == 0.0 && packedDepth.y == 0.0;
}

bool isNearClipPacked(const in vec2 packedDepth) {
    return packedDepth.x == 0.0 && packedDepth.y == 1.0;
}

bool outsideBounds(const in vec2 p) {
    return p.x < uBounds.x || p.y < uBounds.y || p.x > uBounds.z || p.y > uBounds.w;
}

// view space size of a pixel at view z
float getPixelSize(const in float viewZ) {
    #if dOrthographic == 1
        return 2.0 * uInvProjection[0][0] / uTexSize.x;
    #else
        return 2.0 * uInvProjection[0][0] * -viewZ / uTexSize.x;
    #endif
}

void main(void) {
    vec2 coords = gl_FragCoord.xy / uTexSize;

    vec4 selfSsaoDepth = texture2D(tSsaoDepth, coords);
    vec2 packedDepth = selfSsaoDepth.zw;

    if (outsideBounds(coords)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packedDepth);
        return;
    }

    float selfDepth = unpackRGToUnitInterval(packedDepth);
    // (if background and if second pass) or if near clip
    if ((isBackgroundPacked(packedDepth) && uBlurDirectionY != 0.0) || isNearClipPacked(packedDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packedDepth);
        return;
    }

    float selfViewZ = getViewZ(selfDepth);
    float pixelSize = getPixelSize(selfViewZ);

    vec2 offset = vec2(uBlurDirectionX, uBlurDirectionY) / uTexSize;

    float sum = 0.0;
    float kernelSum = 0.0;
    const int halfKernelSize = dOcclusionKernelSize / 2;
    // only if kernelSize is odd
    for (int i = -halfKernelSize; i <= halfKernelSize; i++) {
        if (abs(float(i)) > 1.0 && abs(float(i)) * pixelSize > 0.8) continue;

        vec2 sampleCoords = coords + float(i) * offset;
        if (outsideBounds(sampleCoords)) {
            continue;
        }

        vec4 sampleSsaoDepth = texture2D(tSsaoDepth, sampleCoords);

        float sampleDepth = unpackRGToUnitInterval(sampleSsaoDepth.zw);
        if (isBackgroundPacked(sampleSsaoDepth.zw) || isNearClipPacked(sampleSsaoDepth.zw)) {
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
    // no sample within the depth bias (a background pixel of the first pass): the pixel unblurred
    gl_FragColor = kernelSum == 0.0 ? selfSsaoDepth : vec4(packUnitIntervalToRG(sum / kernelSum), packedDepth);
}
`;
