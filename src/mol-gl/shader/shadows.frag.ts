/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const shadows_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

#include common

uniform sampler2D tDepth;
uniform vec2 uTexSize;

uniform float uNear;
uniform float uFar;

#if dLightCount != 0
    uniform vec3 uLightDirection[dLightCount];
    uniform vec3 uLightColor[dLightCount];
#endif

uniform mat4 uProjection;
uniform mat4 uInvProjection;

uniform float uMaxDistance;
uniform float uTolerance;
uniform float uBias;

bool isBackground(const in float depth) {
    return depth == 1.0;
}

float getViewZ(in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

float getDepth(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepth, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #endif
}

// based on https://panoskarabelas.com/posts/screen_space_shadows/
float screenSpaceShadow(in vec2 coords, in vec3 position, in vec3 lightDirection, in float stepLength) {
    // Ray position and direction (in view-space)
    vec3 rayPos = position;
    vec3 rayDir = -lightDirection;

    // Compute ray step
    vec3 rayStep = rayDir * stepLength;

    // Ray march towards the light
    float occlusion = 0.0;
    vec4 rayCoords = vec4(0.0);
    for (int i = 0; i < dSteps; ++i) {
        // Step the ray
        rayPos += rayStep;

        // Compute the difference between the ray's and the camera's depth
        rayCoords = uProjection * vec4(rayPos, 1.0);
        rayCoords.xyz = (rayCoords.xyz / rayCoords.w) * 0.5 + 0.5;
        float depth = getDepth(rayCoords.xy);
        float viewZ = getViewZ(depth);
        float zDelta = rayPos.z - viewZ;

        if (zDelta < uTolerance) {
            occlusion = 1.0;
            break;
        }
    }

    // Fade out as we approach the edges of the screen
    vec2 fade = max(12.0 * abs(rayCoords.xy - 0.5) - 5.0, vec2(0.0));
    occlusion *= saturate(1.0 - dot(fade, fade));

    return 1.0 - (uBias * occlusion);
}

void main(void) {
    vec2 invTexSize = 1.0 / uTexSize;
    vec2 selfCoords = gl_FragCoord.xy * invTexSize;

    float selfDepth = getDepth(selfCoords);

    if (isBackground(selfDepth)) {
        gl_FragColor = vec4(0.0);
        return;
    }

    vec3 selfViewPos = screenSpaceToViewSpace(vec3(selfCoords, selfDepth), uInvProjection);
    float stepLength = uMaxDistance / float(dSteps);

    float o = 1.0;
    #if dLightCount != 0
        float sh[dLightCount];
        #pragma unroll_loop_start
        for (int i = 0; i < dLightCount; ++i) {
            sh[i] = screenSpaceShadow(selfCoords, selfViewPos, uLightDirection[i], stepLength);
            o = min(o, sh[i]);
        }
        #pragma unroll_loop_end
    #endif

    gl_FragColor = vec4(o);
}
`;