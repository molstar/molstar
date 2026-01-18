/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

export const outlines_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tDepthOpaque;
uniform sampler2D tDepthTransparent;
uniform vec2 uTexSize;

uniform float uNear;
uniform float uFar;
uniform mat4 uInvProjection;

uniform float uOutlineThreshold;

#include common

float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

float getDepthOpaque(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepthOpaque, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepthOpaque, coords));
    #endif
}

vec2 getDepthTransparentWithAlpha(const in vec2 coords) {
    #ifdef dTransparentOutline
        return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords));
    #else
        return vec2(1.0, 0.0);
    #endif
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

float getPixelSize(const in vec2 coords, const in float depth) {
    vec3 viewPos0 = screenSpaceToViewSpace(vec3(coords, depth), uInvProjection);
    vec3 viewPos1 = screenSpaceToViewSpace(vec3(coords + vec2(1.0, 0.0) / uTexSize, depth), uInvProjection);
    return distance(viewPos0, viewPos1);
}

void main(void) {
    float backgroundViewZ = 2.0 * uFar;

    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec2 invTexSize = 1.0 / uTexSize;

    float selfDepthOpaque = getDepthOpaque(coords);
    float selfViewZOpaque = isBackground(selfDepthOpaque) ? backgroundViewZ : getViewZ(selfDepthOpaque);
    float pixelSizeOpaque = getPixelSize(coords, selfDepthOpaque) * uOutlineThreshold;

    vec2 selfDepthTransparentWithAlpha = getDepthTransparentWithAlpha(coords);
    float selfDepthTransparent = selfDepthTransparentWithAlpha.x;
    float selfViewZTransparent = isBackground(selfDepthTransparent) ? backgroundViewZ : getViewZ(selfDepthTransparent);
    float pixelSizeTransparent = getPixelSize(coords, selfDepthTransparent) * uOutlineThreshold;

    float bestOpaqueDepth = 1.0;
    float bestTransparentDepth = 1.0;
    float bestTransparentAlpha = 0.0;

    float opaqueOutlineFlag = 0.0;
    float transparentOutlineFlag = 0.0;

    for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
            vec2 sampleCoords = coords + vec2(float(x), float(y)) * invTexSize;

            // Opaque
            float sampleDepthOpaque = getDepthOpaque(sampleCoords);
            float sampleViewZOpaque = isBackground(sampleDepthOpaque) ? backgroundViewZ : getViewZ(sampleDepthOpaque);
            if (abs(selfViewZOpaque - sampleViewZOpaque) > pixelSizeOpaque && selfDepthOpaque > sampleDepthOpaque && sampleDepthOpaque <= bestOpaqueDepth) {
                bestOpaqueDepth = sampleDepthOpaque;
                opaqueOutlineFlag = 1.0;
            }

            // Transparent
            vec2 sampleDepthTransparentWithAlpha = getDepthTransparentWithAlpha(sampleCoords);
            float sampleDepthTransparent = sampleDepthTransparentWithAlpha.x;
            float sampleAlphaTransparent = sampleDepthTransparentWithAlpha.y;
            float sampleViewZTransparent = isBackground(sampleDepthTransparent) ? backgroundViewZ : getViewZ(sampleDepthTransparent);
            if (abs(selfViewZTransparent - sampleViewZTransparent) > pixelSizeTransparent && selfDepthTransparent > sampleDepthTransparent && sampleDepthTransparent <= bestTransparentDepth) {
                bestTransparentDepth = sampleDepthTransparent;
                bestTransparentAlpha = sampleAlphaTransparent;
                transparentOutlineFlag = 1.0;
            }
        }
    }

    // lazy curvature veto: only compute if an outline was found and center is not background
    const float kCurvatureGate = 0.75;
    vec2 dx = vec2(invTexSize.x, 0.0);
    vec2 dy = vec2(0.0, invTexSize.y);

    if (opaqueOutlineFlag > 0.0 && !isBackground(selfDepthOpaque)) {
        float dL = getDepthOpaque(coords - dx);
        float dR = getDepthOpaque(coords + dx);
        float dU = getDepthOpaque(coords + dy);
        float dD = getDepthOpaque(coords - dy);

        float vzL = isBackground(dL) ? backgroundViewZ : getViewZ(dL);
        float vzR = isBackground(dR) ? backgroundViewZ : getViewZ(dR);
        float vzU = isBackground(dU) ? backgroundViewZ : getViewZ(dU);
        float vzD = isBackground(dD) ? backgroundViewZ : getViewZ(dD);

        float ddx = abs(vzL + vzR - 2.0 * selfViewZOpaque);
        float ddy = abs(vzU + vzD - 2.0 * selfViewZOpaque);
        float curvOpaque = max(ddx, ddy);

        if (curvOpaque < pixelSizeOpaque * kCurvatureGate) {
            opaqueOutlineFlag = 0.0;
            bestOpaqueDepth = 1.0;
        }
    }

    if (transparentOutlineFlag > 0.0 && !isBackground(selfDepthTransparent)) {
        vec2 daL = getDepthTransparentWithAlpha(coords - dx);
        vec2 daR = getDepthTransparentWithAlpha(coords + dx);
        vec2 daU = getDepthTransparentWithAlpha(coords + dy);
        vec2 daD = getDepthTransparentWithAlpha(coords - dy);

        float dL = daL.x;
        float dR = daR.x;
        float dU = daU.x;
        float dD = daD.x;

        float vzL = isBackground(dL) ? backgroundViewZ : getViewZ(dL);
        float vzR = isBackground(dR) ? backgroundViewZ : getViewZ(dR);
        float vzU = isBackground(dU) ? backgroundViewZ : getViewZ(dU);
        float vzD = isBackground(dD) ? backgroundViewZ : getViewZ(dD);

        float ddx = abs(vzL + vzR - 2.0 * selfViewZTransparent);
        float ddy = abs(vzU + vzD - 2.0 * selfViewZTransparent);
        float curvTransparent = max(ddx, ddy);

        if (curvTransparent < pixelSizeTransparent * kCurvatureGate) {
            transparentOutlineFlag = 0.0;
            bestTransparentDepth = 1.0;
            bestTransparentAlpha = 0.0;
        }
    }

    if (transparentOutlineFlag > 0.0 && bestOpaqueDepth < 1.0 && bestTransparentDepth > bestOpaqueDepth) {
        transparentOutlineFlag = 0.0;
        bestTransparentAlpha = 0.0;
    }

    vec2 depthPacked; // Pack depth in G/B channels
    float outlineTypeFlag = 0.0;
    if (opaqueOutlineFlag > 0.0 && transparentOutlineFlag > 0.0) {
        outlineTypeFlag = 0.75; // Both
        depthPacked = packUnitIntervalToRG(bestOpaqueDepth);
    } else if (transparentOutlineFlag > 0.0) {
        outlineTypeFlag = 0.5;  // Transparent only
        depthPacked = packUnitIntervalToRG(bestTransparentDepth);
    } else if (opaqueOutlineFlag > 0.0) {
        outlineTypeFlag = 0.25; // Opaque only
        depthPacked = packUnitIntervalToRG(bestOpaqueDepth);
    }

    float alpha = clamp(bestTransparentAlpha, 0.0, 0.5) * 2.0; // limiting to range [0.0, 0.5] to improve alpha precision since we don't need a wider range
    float packedFlagWithAlpha = pack2x4(vec2(outlineTypeFlag, alpha)); // pack outlineType with alpha
    gl_FragColor = vec4(packedFlagWithAlpha, depthPacked.x, depthPacked.y, bestTransparentDepth);
}
`;