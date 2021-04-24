/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform vec3 uSamples[dNSamples];

uniform mat4 uProjection;
uniform mat4 uInvProjection;

uniform float uRadius;
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

bool outsideBounds(const in vec2 p) {
    return p.x < uBounds.x || p.y < uBounds.y || p.x > uBounds.z || p.y > uBounds.w;
}

float getDepth(const in vec2 coords) {
    return outsideBounds(coords) ? 1.0 : unpackRGBAToDepth(texture2D(tDepth, coords));
}

vec3 normalFromDepth(const in float depth, const in float depth1, const in float depth2, vec2 offset1, vec2 offset2) {
    vec3 p1 = vec3(offset1, depth1 - depth);
    vec3 p2 = vec3(offset2, depth2 - depth);

    vec3 normal = cross(p1, p2);
    normal.z = -normal.z;

    return normalize(normal);
}

// StarCraft II Ambient Occlusion by [Filion and McNaughton 2008]
void main(void) {
    vec2 invTexSize = 1.0 / uTexSize;
    vec2 selfCoords = gl_FragCoord.xy * invTexSize;

    float selfDepth = getDepth(selfCoords);
    vec2 selfPackedDepth = packUnitIntervalToRG(selfDepth);

    if (isBackground(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(0.0), selfPackedDepth);
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
    for(int i = 0; i < dNSamples; i++){
        vec3 sampleViewPos = TBN * uSamples[i];
        sampleViewPos = selfViewPos + sampleViewPos * uRadius;

        vec4 offset = vec4(sampleViewPos, 1.0);
        offset = uProjection * offset;
        offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;

        float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, getDepth(offset.xy)), uInvProjection).z;

        occlusion += step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, uRadius / abs(selfViewPos.z - sampleViewZ));
    }
    occlusion = 1.0 - (uBias * occlusion / float(dNSamples));

    vec2 packedOcclusion = packUnitIntervalToRG(occlusion);

    gl_FragColor = vec4(packedOcclusion, selfPackedDepth);
}
`;