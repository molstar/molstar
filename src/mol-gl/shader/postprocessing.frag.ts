/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const postprocessing_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tSsaoDepth;
uniform sampler2D tColor;
uniform sampler2D tDepthOpaque;
uniform sampler2D tDepthTransparent;
uniform sampler2D tShadows;
uniform sampler2D tOutlines;
uniform vec2 uTexSize;

uniform float uNear;
uniform float uFar;
uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;
uniform vec3 uOutlineColor;
uniform bool uTransparentBackground;

uniform vec2 uOcclusionOffset;

uniform float uMaxPossibleViewZDiff;
uniform mat4 uInvProjection;

const float outlineDistanceFactor = 5.0;
const vec3 occlusionColor = vec3(0.0);

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
    #ifdef dTransparentOutline
        return unpackRGBAToDepth(texture2D(tDepthTransparent, coords));
    #else
        return 1.0;
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

float getOutline(const in vec2 coords, const in float opaqueDepth, out float closestTexel) {
    float backgroundViewZ = uFar + 3.0 * uMaxPossibleViewZDiff;
    vec2 invTexSize = 1.0 / uTexSize;

    float selfDepth = min(opaqueDepth, getDepthTransparent(coords));
    float selfViewZ = isBackground(selfDepth) ? backgroundViewZ : getViewZ(selfDepth);
    float pixelSize = getPixelSize(coords, selfDepth);

    float outline = 1.0;
    closestTexel = 1.0;
    for (int y = -dOutlineScale; y <= dOutlineScale; y++) {
        for (int x = -dOutlineScale; x <= dOutlineScale; x++) {
            if (x * x + y * y > dOutlineScale * dOutlineScale) {
                continue;
            }

            vec2 sampleCoords = coords + vec2(float(x), float(y)) * invTexSize;

            vec4 sampleOutlineCombined = texture2D(tOutlines, sampleCoords);
            float sampleOutline = sampleOutlineCombined.r;
            float sampleOutlineDepth = unpackRGToUnitInterval(sampleOutlineCombined.gb);
            float sampleOutlineViewZ = isBackground(sampleOutlineDepth) ? backgroundViewZ : getViewZ(sampleOutlineDepth);

            if (sampleOutline == 0.0 && sampleOutlineDepth < closestTexel && abs(selfViewZ - sampleOutlineViewZ) > uMaxPossibleViewZDiff + (pixelSize * outlineDistanceFactor)) {
                outline = 0.0;
                closestTexel = sampleOutlineDepth;
            }
        }
    }
    return closestTexel < opaqueDepth ? outline : 1.0;
}

float getSsao(vec2 coords) {
    float rawSsao = unpackRGToUnitInterval(texture2D(tSsaoDepth, coords).xy);
    if (rawSsao > 0.999) {
        return 1.0;
    } else if (rawSsao > 0.001) {
        return rawSsao;
    }
    // treat values close to 0.0 as errors and return no occlusion
    return 1.0;
}

void main(void) {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec4 color = texture2D(tColor, coords);

    float viewDist;
    float fogFactor;
    float opaqueDepth = getDepthOpaque(coords);

    #ifdef dOcclusionEnable
        if (!isBackground(opaqueDepth)) {
            viewDist = abs(getViewZ(opaqueDepth));
            fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
            float occlusionFactor = getSsao(coords + uOcclusionOffset);
            if (!uTransparentBackground) {
                color.rgb = mix(mix(occlusionColor, uFogColor, fogFactor), color.rgb, occlusionFactor);
            } else {
                color.rgb = mix(occlusionColor * (1.0 - fogFactor), color.rgb, occlusionFactor);
            }
        }
    #endif

    #ifdef dShadowEnable
        if (!isBackground(opaqueDepth)) {
            viewDist = abs(getViewZ(opaqueDepth));
            fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
            vec4 shadow = texture2D(tShadows, coords);
            if (!uTransparentBackground) {
                color.rgb = mix(mix(vec3(0), uFogColor, fogFactor), color.rgb, shadow.a);
            } else {
                color.rgb = mix(vec3(0) * (1.0 - fogFactor), color.rgb, shadow.a);
            }
        }
    #endif

    // outline needs to be handled after occlusion and shadow to keep them clean
    #ifdef dOutlineEnable
        float closestTexel;
        float outline = getOutline(coords, opaqueDepth, closestTexel);
        if (outline == 0.0) {
            viewDist = abs(getViewZ(closestTexel));
            fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
            if (!uTransparentBackground) {
                color.rgb = mix(uOutlineColor, uFogColor, fogFactor);
            } else {
                color.a = 1.0 - fogFactor;
                color.rgb = mix(uOutlineColor, color.rgb, fogFactor);
            }
        }
    #endif

    gl_FragColor = color;
}
`;