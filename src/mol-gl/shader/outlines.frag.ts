/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const outlines_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tDepth;
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

float getDepth(const in vec2 coords) {
    return unpackRGBAToDepth(texture2D(tDepth, coords));
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

void main(void) {
    float backgroundViewZ = uFar + 3.0 * uMaxPossibleViewZDiff;

    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec2 invTexSize = 1.0 / uTexSize;

    float selfDepth = getDepth(coords);
    float selfViewZ = isBackground(selfDepth) ? backgroundViewZ : getViewZ(getDepth(coords));

    float outline = 1.0;
    float bestDepth = 1.0;

    for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
            vec2 sampleCoords = coords + vec2(float(x), float(y)) * invTexSize;
            float sampleDepth = getDepth(sampleCoords);
            float sampleViewZ = isBackground(sampleDepth) ? backgroundViewZ : getViewZ(sampleDepth);

            if (abs(selfViewZ - sampleViewZ) > uMaxPossibleViewZDiff && selfDepth > sampleDepth && sampleDepth <= bestDepth) {
                outline = 0.0;
                bestDepth = sampleDepth;
            }
        }
    }

    gl_FragColor = vec4(outline, packUnitIntervalToRG(bestDepth), 0.0);
}
`;