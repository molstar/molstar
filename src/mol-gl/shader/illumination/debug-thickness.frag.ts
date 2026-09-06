/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const debugThickness_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tColor;
uniform sampler2D tThickness;
uniform sampler2D tDepth;
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform float uNear;
uniform float uFar;

uniform float uMinThickness;
uniform float uThicknessFactor;

#include common

float getDepth(const in vec2 coords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    return texture2D(tDepth, c).r;
}

float getBackDepth(const in vec2 coords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    return unpackRGBAToDepth(texture2D(tThickness, c));
}

float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

float linearizeDepth(const in float depth) {
    return saturate((abs(getViewZ(depth)) - uNear) / (uFar - uNear));
}

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;

    // three panels: front depth | back depth | thickness
    float panel = floor(coords.x * 3.0);
    vec2 sampleCoords = vec2(fract(coords.x * 3.0), coords.y);

    float frontDepth = getDepth(sampleCoords);
    float backDepth = getBackDepth(sampleCoords);

    vec3 color;
    if (panel < 0.5) {
        color = vec3(linearizeDepth(frontDepth));
    } else if (panel < 1.5) {
        color = vec3(linearizeDepth(backDepth));
    } else {
        if (frontDepth == 1.0) {
            color = vec3(0.0);
        } else {
            // same estimation as the auto thickness mode in trace.frag
            float thickness = max(uMinThickness, (getViewZ(frontDepth) - getViewZ(backDepth)) * uThicknessFactor * texture2D(tColor, sampleCoords).a);
            // tone-map to handle a wide range of thickness values
            color = vec3(1.0 - exp(-thickness * 0.1));
        }
    }

    // panel separators
    if (abs(gl_FragCoord.x - uTexSize.x / 3.0) < 1.0 || abs(gl_FragCoord.x - 2.0 * uTexSize.x / 3.0) < 1.0) {
        color = vec3(1.0, 0.5, 0.0);
    }

    gl_FragColor = vec4(color, 1.0);
}
`;
