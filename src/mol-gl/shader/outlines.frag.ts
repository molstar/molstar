/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 @author Alexander Rose <alexander.rose@weirdbyte.de>
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

uniform float uMaxPossibleViewZDiff;

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

float getDepthTransparent(const in vec2 coords) {
    return unpackRGBAToDepth(texture2D(tDepthTransparent, coords));
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

void main(void) {
    float backgroundViewZ = uFar + 3.0 * uMaxPossibleViewZDiff;

    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec2 invTexSize = 1.0 / uTexSize;

    float selfDepthOpaque = getDepthOpaque(coords);
    float selfViewZOpaque = isBackground(selfDepthOpaque) ? backgroundViewZ : getViewZ(selfDepthOpaque);

    float selfDepthTransparent = getDepthTransparent(coords);
    float selfViewZTransparent = isBackground(selfDepthTransparent) ? backgroundViewZ : getViewZ(selfDepthTransparent);

    float outline = 1.0;
    float bestDepth = 1.0;

    for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
            vec2 sampleCoords = coords + vec2(float(x), float(y)) * invTexSize;

            float sampleDepthOpaque = getDepthOpaque(sampleCoords);
            float sampleDepthTransparent = getDepthTransparent(sampleCoords);

            float sampleViewZOpaque = isBackground(sampleDepthOpaque) ? backgroundViewZ : getViewZ(sampleDepthOpaque);
            if (abs(selfViewZOpaque - sampleViewZOpaque) > uMaxPossibleViewZDiff && selfDepthOpaque > sampleDepthOpaque && sampleDepthOpaque <= bestDepth) {
                outline = 0.0;
                bestDepth = sampleDepthOpaque;
            }

            if (sampleDepthTransparent < sampleDepthOpaque) {
                float sampleViewZTransparent = isBackground(sampleDepthTransparent) ? backgroundViewZ : getViewZ(sampleDepthTransparent);
                if (abs(selfViewZTransparent - sampleViewZTransparent) > uMaxPossibleViewZDiff && selfDepthTransparent > sampleDepthTransparent && sampleDepthTransparent <= bestDepth) {
                    outline = 0.0;
                    bestDepth = sampleDepthTransparent;
                }
            }
        }
    }

    gl_FragColor = vec4(outline, packUnitIntervalToRG(bestDepth), 0.0);
}
`;