/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { common } from './chunks/common.glsl.js';

export const ssao_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

${common}

float smootherstep(float edge0, float edge1, float x) {
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

// ordered dither of a 2x2 and a 4x4 block (values k / 4 and k / 16), \`pixel\` in whole pixels
float bayer2(vec2 pixel) {
    pixel = floor(pixel);
    return fract(dot(pixel, vec2(0.5, pixel.y * 0.75)));
}

float bayer4(const in vec2 pixel) {
    return bayer2(0.5 * pixel) * 0.25 + bayer2(pixel);
}

// rotation of the sample kernel at \`pixel\`: the 16 angles of a 4x4 block in ordered dither, so that the blur
// (wider than a block) averages all of them
vec2 getKernelRotation(const in vec2 pixel) {
    float angle = 2.0 * PI * (bayer4(pixel) + 0.5 / 16.0);
    return vec2(cos(angle), sin(angle));
}

// depth is reversed: the background (far plane) at 0
bool isBackground(const in float depth) {
    return depth == 0.0;
}

// inset by half a texel to avoid sampling texels outside the viewport
vec2 clampToBounds(const in vec2 coords) {
    vec2 halfTexel = 0.5 / uTexSize;
    return clamp(coords, uBounds.xy + halfTexel, uBounds.zw - halfTexel);
}

float getDepth(const in vec2 coords, const in int transparentFlag) {
    if (transparentFlag == 1){
        #if defined(dIncludeTransparent)
            return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords)).x;
        #else
            return 0.0;
        #endif
    } else {
        #ifdef dPackedDepth
            return unpackRGBAToDepth(texture2D(tDepth, coords));
        #else
            return texture2D(tDepth, coords).r;
        #endif
    }
}

#if defined(dIncludeTransparent)
    vec2 getDepthTransparentWithAlpha(const in vec2 coords){
        return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords));
    }
#endif

#define dQuarterThreshold 0.1
#define dHalfThreshold 0.05

float getLevelDepth(const in sampler2D tex, const in vec2 coords) {
    #ifdef dPackedDepthLevels
        return unpackRGBAToDepth(texture2D(tex, coords));
    #else
        return texture2D(tex, coords).r;
    #endif
}

float getMappedDepth(const in vec2 coords, const in vec2 selfCoords) {
    float d = distance(coords, selfCoords);
    if (d > dQuarterThreshold) {
        return getLevelDepth(tDepthQuarter, coords);
    } else if (d > dHalfThreshold) {
        return getLevelDepth(tDepthHalf, coords);
    } else {
        return getDepth(coords, 0);
    }
}

#if defined(dIncludeTransparent)
    vec2 getMappedDepthTransparentWithAlpha(const in vec2 coords, const in vec2 selfCoords) {
        float d = distance(coords, selfCoords);
        if (d > dQuarterThreshold) {
            return unpackRGBAToDepthWithAlpha(texture2D(tDepthQuarterTransparent, coords));
        } else if (d > dHalfThreshold) {
            return unpackRGBAToDepthWithAlpha(texture2D(tDepthHalfTransparent, coords));
        } else {
            return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords));
        }
    }
#endif

// adapted from https://gist.github.com/bgolus/a07ed65602c009d5e2f753826e8078a0
vec3 viewNormalAtPixelPositionAccurate(const in vec2 vpos, const in int transparentFlag) {
    // current pixel's depth
    float c = getDepth(vpos, transparentFlag);

    // get current pixel's view space position
    vec3 viewSpacePos_c = screenSpaceToViewSpace(vec3(vpos, c), uInvProjection);

    // get view space position at 1 pixel offsets in each major direction
    vec3 viewSpacePos_l = screenSpaceToViewSpace(vec3(vpos + vec2(-1.0, 0.0) / uTexSize, getDepth(vpos + vec2(-1.0, 0.0) / uTexSize, transparentFlag)), uInvProjection);
    vec3 viewSpacePos_r = screenSpaceToViewSpace(vec3(vpos + vec2( 1.0, 0.0) / uTexSize, getDepth(vpos + vec2( 1.0, 0.0) / uTexSize, transparentFlag)), uInvProjection);
    vec3 viewSpacePos_d = screenSpaceToViewSpace(vec3(vpos + vec2( 0.0,-1.0) / uTexSize, getDepth(vpos + vec2( 0.0,-1.0) / uTexSize, transparentFlag)), uInvProjection);
    vec3 viewSpacePos_u = screenSpaceToViewSpace(vec3(vpos + vec2( 0.0, 1.0) / uTexSize, getDepth(vpos + vec2( 0.0, 1.0) / uTexSize, transparentFlag)), uInvProjection);

    // get the difference between the current and each offset position
    vec3 l = viewSpacePos_c - viewSpacePos_l;
    vec3 r = viewSpacePos_r - viewSpacePos_c;
    vec3 d = viewSpacePos_c - viewSpacePos_d;
    vec3 u = viewSpacePos_u - viewSpacePos_c;

    // get depth values at 1 & 2 pixels offsets from current along the horizontal axis
    vec4 H = vec4(
        getDepth(vpos + vec2(-1.0, 0.0) / uTexSize, transparentFlag),
        getDepth(vpos + vec2( 1.0, 0.0) / uTexSize, transparentFlag),
        getDepth(vpos + vec2(-2.0, 0.0) / uTexSize, transparentFlag),
        getDepth(vpos + vec2( 2.0, 0.0) / uTexSize, transparentFlag)
    );

    // get depth values at 1 & 2 pixels offsets from current along the vertical axis
    vec4 V = vec4(
        getDepth(vpos + vec2(0.0,-1.0) / uTexSize, transparentFlag),
        getDepth(vpos + vec2(0.0, 1.0) / uTexSize, transparentFlag),
        getDepth(vpos + vec2(0.0,-2.0) / uTexSize, transparentFlag),
        getDepth(vpos + vec2(0.0, 2.0) / uTexSize, transparentFlag)
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

float getPixelSize(const in vec2 coords, const in float depth) {
    vec3 viewPos0 = screenSpaceToViewSpace(vec3(coords, depth), uInvProjection);
    vec3 viewPos1 = screenSpaceToViewSpace(vec3(coords + vec2(1.0, 0.0) / uTexSize, depth), uInvProjection);
    return distance(viewPos0, viewPos1);
}

// StarCraft II Ambient Occlusion by [Filion and McNaughton 2008]
void main(void) {
    vec2 invTexSize = 1.0 / uTexSize;
    vec2 selfCoords = gl_FragCoord.xy * invTexSize;
    float selfDepth = getDepth(selfCoords, uTransparencyFlag);
    if (isBackground(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), packUnitIntervalToRG(0.0));
        return;
    }

    vec3 selfViewNormal = viewNormalAtPixelPositionAccurate(selfCoords, uTransparencyFlag);
    vec3 selfViewPos = screenSpaceToViewSpace(vec3(selfCoords, selfDepth), uInvProjection);
    // linear depth for the blur, reversed (the near plane at 1, the far plane and the background at 0): 16-bit fixed
    // point keeps a constant view z precision, unlike the window depth
    vec2 selfPackedDepth = packUnitIntervalToRG(clamp((uFar + selfViewPos.z) / (uFar - uNear), 0.0, 1.0));

    vec3 randomVec = vec3(getKernelRotation(gl_FragCoord.xy), 0.0);
    vec3 tangent = normalize(randomVec - selfViewNormal * dot(randomVec, selfViewNormal));
    vec3 bitangent = cross(selfViewNormal, tangent);
    mat3 TBN = mat3(tangent, bitangent, selfViewNormal);

    float occlusion = 0.0;
    #ifdef dMultiScale
        float pixelSize = getPixelSize(selfCoords, selfDepth);

        for(int l = 0; l < dLevels; l++) {
            // TODO: smooth transition
            if (pixelSize * uNearThreshold > uLevelRadius[l]) continue;
            if (pixelSize * uFarThreshold < uLevelRadius[l]) continue;

            float levelOcclusion = 0.0;
            for(int i = 0; i < dNSamples; i++) {
                // get sample position:
                vec3 sampleViewPos = TBN * uSamples[i];
                sampleViewPos = selfViewPos + sampleViewPos * uLevelRadius[l];

                // project sample position:
                vec4 offset = vec4(sampleViewPos, 1.0);
                offset = uProjection * offset;
                offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;
                offset.xy = clampToBounds(offset.xy);

                // get sample depth:
                float sampleOcc = 0.0;
                #ifdef dIllumination
                    if (uTransparencyFlag == 1) {
                #endif
                        float sampleDepth = getMappedDepth(offset.xy, selfCoords);
                        if (!isBackground(sampleDepth)) {
                            float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, sampleDepth), uInvProjection).z;
                            sampleOcc = step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, uLevelRadius[l] / abs(selfViewPos.z - sampleViewZ)) * uLevelBias[l];
                        }
                #ifdef dIllumination
                    }
                #endif
                #if defined(dIncludeTransparent)
                    vec2 sampleDepthWithAlpha = getMappedDepthTransparentWithAlpha(offset.xy, selfCoords);
                    if (!isBackground(sampleDepthWithAlpha.x))  {
                        float sampleViewZWithAlpha = screenSpaceToViewSpace(vec3(offset.xy, sampleDepthWithAlpha.x), uInvProjection).z;
                        sampleOcc = max(sampleOcc, step(sampleViewPos.z + 0.025, sampleViewZWithAlpha) * smootherstep(0.0, 1.0, uLevelRadius[l] / abs(selfViewPos.z - sampleViewZWithAlpha)) * uLevelBias[l] * sampleDepthWithAlpha.y);
                    }
                #endif

                levelOcclusion += sampleOcc;
            }
            levelOcclusion /= float(dNSamples);
            occlusion = max(occlusion, levelOcclusion);
        }
    #else
        for(int i = 0; i < dNSamples; i++) {
            vec3 sampleViewPos = TBN * uSamples[i];
            sampleViewPos = selfViewPos + sampleViewPos * uRadius;

            vec4 offset = vec4(sampleViewPos, 1.0);
            offset = uProjection * offset;
            offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;
            offset.xy = clampToBounds(offset.xy);

            float sampleOcc = 0.0;
            #ifdef dIllumination
                if (uTransparencyFlag == 1) {
            #endif
                    // NOTE: using getMappedDepth here causes issues on some mobile devices
                    float sampleDepth = getDepth(offset.xy, 0);
                    if (!isBackground(sampleDepth)) {
                        float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, sampleDepth), uInvProjection).z;
                        sampleOcc = step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, uRadius / abs(selfViewPos.z - sampleViewZ));
                    }
            #ifdef dIllumination
                }
            #endif
            #if defined(dIncludeTransparent)
                vec2 sampleDepthWithAlpha = getDepthTransparentWithAlpha(offset.xy);
                if (!isBackground(sampleDepthWithAlpha.x)) {
                    float sampleViewZWithAlpha = screenSpaceToViewSpace(vec3(offset.xy, sampleDepthWithAlpha.x), uInvProjection).z;
                    sampleOcc = max(sampleOcc, step(sampleViewPos.z + 0.025, sampleViewZWithAlpha) * smootherstep(0.0, 1.0, uRadius / abs(selfViewPos.z - sampleViewZWithAlpha)) * sampleDepthWithAlpha.y);
                }
            #endif

            occlusion += sampleOcc;
        }

        occlusion /= float(dNSamples);
    #endif
    occlusion = 1.0 - (uBias * occlusion);

    vec2 packedOcclusion = packUnitIntervalToRG(clamp(occlusion, 0.01, 1.0));

    gl_FragColor = vec4(packedOcclusion, selfPackedDepth);
}
`;
