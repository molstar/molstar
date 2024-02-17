/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const ssao_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

#include common

uniform sampler2D tDepth;
uniform sampler2D tDepthHalf;
uniform sampler2D tDepthQuarter;
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform vec3 uSamples[dNSamples];

uniform mat4 uProjection;
uniform mat4 uInvProjection;

#ifdef dMultiScale
    uniform float uLevelRadius[dLevels];
    uniform float uLevelBias[dLevels];
    uniform float uNearThreshold;
    uniform float uFarThreshold;
#else
    uniform float uRadius;
#endif
uniform float uBias;

float smootherstep(float edge0, float edge1, float x) {
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

float noise(const in vec2 coords) {
    float a = 12.9898;
    float b = 78.233;
    float c = 43758.5453;
    float dt = dot(coords, vec2(a,b));
    float sn = mod(dt, PI);
    return abs(fract(sin(sn) * c)); // is abs necessary?
}

vec2 getNoiseVec2(const in vec2 coords) {
    return vec2(noise(coords), noise(coords + vec2(PI, 2.71828)));
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

float getDepth(const in vec2 coords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    #ifdef depthTextureSupport
        return texture2D(tDepth, c).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepth, c));
    #endif
}

#define dQuarterThreshold 0.1
#define dHalfThreshold 0.05

float getMappedDepth(const in vec2 coords, const in vec2 selfCoords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    float d = distance(coords, selfCoords);
    #ifdef depthTextureSupport
        if (d > dQuarterThreshold) {
            return texture2D(tDepthQuarter, c).r;
        } else if (d > dHalfThreshold) {
            return texture2D(tDepthHalf, c).r;
        } else {
            return texture2D(tDepth, c).r;
        }
    #else
        if (d > dQuarterThreshold) {
            return unpackRGBAToDepth(texture2D(tDepthQuarter, c));
        } else if (d > dHalfThreshold) {
            return unpackRGBAToDepth(texture2D(tDepthHalf, c));
        } else {
            return unpackRGBAToDepth(texture2D(tDepth, c));
        }
    #endif
}

vec3 normalFromDepth(const in float depth, const in float depth1, const in float depth2, vec2 offset1, vec2 offset2) {
    vec3 p1 = vec3(offset1, depth1 - depth);
    vec3 p2 = vec3(offset2, depth2 - depth);

    vec3 normal = cross(p1, p2);
    normal.z = -normal.z;

    return normalize(normal);
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

    float selfDepth = getDepth(selfCoords);
    vec2 selfPackedDepth = packUnitIntervalToRG(selfDepth);

    if (isBackground(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(1.0), selfPackedDepth);
        return;
    }

    vec2 offset1 = vec2(0.0, invTexSize.y);
    vec2 offset2 = vec2(invTexSize.x, 0.0);

    float selfDepth1 = getDepth(selfCoords + offset1);
    float selfDepth2 = getDepth(selfCoords + offset2);

    vec3 selfViewNormal = normalFromDepth(selfDepth, selfDepth1, selfDepth2, offset1, offset2);
    vec3 selfViewPos = screenSpaceToViewSpace(vec3(selfCoords, selfDepth), uInvProjection);

    vec3 randomVec = normalize(vec3(getNoiseVec2(selfCoords) * 2.0 - 1.0, 0.0));
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
                vec3 sampleViewPos = TBN * uSamples[i];
                sampleViewPos = selfViewPos + sampleViewPos * uLevelRadius[l];

                vec4 offset = vec4(sampleViewPos, 1.0);
                offset = uProjection * offset;
                offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;

                float sampleDepth = getMappedDepth(offset.xy, selfCoords);
                float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, sampleDepth), uInvProjection).z;

                levelOcclusion += step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, uLevelRadius[l] / abs(selfViewPos.z - sampleViewZ)) * uLevelBias[l];
            }
            occlusion = max(occlusion, levelOcclusion);
        }
    #else
        for(int i = 0; i < dNSamples; i++) {
            vec3 sampleViewPos = TBN * uSamples[i];
            sampleViewPos = selfViewPos + sampleViewPos * uRadius;

            vec4 offset = vec4(sampleViewPos, 1.0);
            offset = uProjection * offset;
            offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;

            float sampleDepth = getMappedDepth(offset.xy, selfCoords);
            float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, sampleDepth), uInvProjection).z;

            occlusion += step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, uRadius / abs(selfViewPos.z - sampleViewZ));
        }
    #endif
    occlusion = 1.0 - (uBias * occlusion / float(dNSamples));

    vec2 packedOcclusion = packUnitIntervalToRG(clamp(occlusion, 0.01, 1.0));

    gl_FragColor = vec4(packedOcclusion, selfPackedDepth);
}
`;