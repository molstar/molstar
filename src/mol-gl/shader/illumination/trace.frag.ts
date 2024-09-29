/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const trace_frag = `
precision highp int;
precision highp float;
precision highp sampler2D;

uniform sampler2D tColor;
uniform sampler2D tNormal;
uniform sampler2D tShaded;
uniform sampler2D tThickness;
uniform sampler2D tAccumulate;
uniform sampler2D tDepth;
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform float uNear;
uniform float uFar;
uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

#if dLightCount != 0
    uniform vec3 uLightDirection[dLightCount];
    uniform vec3 uLightColor[dLightCount];
#endif
uniform vec3 uAmbientColor;
uniform vec3 uLightStrength;

uniform int uFrameNo;

uniform float uRayDistance;
uniform float uMinThickness;
uniform float uThicknessFactor;
uniform float uThickness;

uniform float uShadowSoftness;
uniform float uShadowThickness;

uniform mat4 uProjection;
uniform mat4 uInvProjection;

#include common

// parts adapted from
// - https://blog.demofox.org/2020/05/25/casual-shadertoy-path-tracing-1-basic-camera-diffuse-emissive/
// - https://github.com/0beqz/realism-effects/blob/v2-debug/src/ssgi/shader/ssgi.frag

//

// after a hit, it moves the ray this far along the normal away from a surface.
// Helps prevent incorrect intersections when rays bounce off of objects.
#define RayPosNormalNudge 0.0001

#if __VERSION__ == 100
    #define StateType float

    // from https://www.shadertoy.com/view/4djSRW
    float hash14(vec4 p4) {
        p4 = fract(p4  * vec4(0.1031, 0.1030, 0.0973, 0.1099));
        p4 += dot(p4, p4.wzxy + 33.33);
        return fract((p4.x + p4.y) * (p4.z + p4.w));
    }

    float randomFloat(inout float state) {
        state += 0.06711056;
        return 1.0 - hash14(vec4(gl_FragCoord.xy, float(uFrameNo), state));
    }
#else
    #define StateType uint

    // https://www.pcg-random.org/
    // https://jcgt.org/published/0009/03/02/
    uint pcg(inout uint seed) {
        seed = seed * 747796405u + 2891336453u;
        uint word = ((seed >> ((seed >> 28u) + 4u)) ^ seed) * 277803737u;
        return (word >> 22u) ^ word;
    }

    float randomFloat(inout uint state) {
        return float(pcg(state)) / 4294967296.0;
    }
#endif

vec3 randomUnitVector(inout StateType state) {
    float z = randomFloat(state) * 2.0 - 1.0;
    float a = randomFloat(state) * TWO_PI;
    float r = sqrt(1.0 - z * z);
    float x = r * cos(a);
    float y = r * sin(a);
    return vec3(x, y, z);
}

struct RayHitInfo {
    bool missed;
    vec3 position;
    vec3 normal;
    vec3 color;
    vec3 emissive;
};

//

float getDepth(const in vec2 coords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    return texture2D(tDepth, c).r;
}

float getThickness(const in vec2 coords) {
    vec2 c = vec2(clamp(coords.x, uBounds.x, uBounds.z), clamp(coords.y, uBounds.y, uBounds.w));
    return unpackRGBAToDepth(texture2D(tThickness, c));
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

vec2 viewSpaceToScreenSpace(const vec3 position) {
    vec4 projectedCoord = uProjection * vec4(position, 1.0);
    projectedCoord.xy /= projectedCoord.w;
    // [-1, 1] --> [0, 1] (NDC to screen position)
    projectedCoord.xy = projectedCoord.xy * 0.5 + 0.5;
    return projectedCoord.xy;
}

vec2 binarySearch(inout vec3 dir, inout vec3 hitPos) {
    float rayHitDepthDifference;
    vec2 coords;

    dir *= 0.5;
    hitPos -= dir;

    for (int i = 0; i < dRefineSteps; i++) {
        coords = viewSpaceToScreenSpace(hitPos);
        float depth = getDepth(coords);
        float z = getViewZ(depth);
        rayHitDepthDifference = z - hitPos.z;

        dir *= 0.5;
        if (rayHitDepthDifference >= 0.0) {
            hitPos -= dir;
        } else {
            hitPos += dir;
        }
    }

    coords = viewSpaceToScreenSpace(hitPos);

    return coords;
}

float calculateGrowthFactor(float begin, float end, float steps) {
    return pow(end / begin, 1.0 / steps);
}

vec2 rayMarch(in vec3 dir, in float thickness, inout vec3 hitPos, out bool missed) {
    float rayHitDepthDifference;
    vec2 coords;

    float begin = float(dFirstStepSize);
    dir *= begin;
    missed = false;
    float gf = calculateGrowthFactor(begin, uRayDistance, float(dSteps));

    for (int i = 1; i < dSteps; i++) {
        hitPos += dir;
        dir *= gf;

        coords = viewSpaceToScreenSpace(hitPos);
        float depth = getDepth(coords);
        float z = getViewZ(depth);
        rayHitDepthDifference = z - hitPos.z;

        if (thickness == 0.0) {
            #ifdef dThicknessMode_auto
                thickness = max(uMinThickness, (getViewZ(getThickness(coords)) - z) * uThicknessFactor * texture2D(tColor, coords).a);
            #else
                thickness = uThickness;
            #endif
        }

        if (rayHitDepthDifference >= 0.0 && rayHitDepthDifference < thickness) {
            if (dRefineSteps == 0) {
                return coords;
            } else {
                return binarySearch(dir, hitPos);
            }
        }
    }

    missed = true;

    return coords;
}

//

void trace(in vec3 rayPos, in vec3 rayDir, inout RayHitInfo hitInfo) {
    vec3 hitPos = vec3(rayPos);
    bool missed;
    vec2 coords;
    coords = rayMarch(rayDir, 0.0, hitPos, missed);

    hitInfo.missed = missed;
    hitInfo.position = hitPos;
    hitInfo.normal = -texture2D(tNormal, coords).rgb;
    hitInfo.color = texture2D(tColor, coords).rgb;
    hitInfo.emissive = texture2D(tColor, coords).rgb * texture2D(tNormal, coords).a * 2.0;

    float depth = getDepth(coords);
    if (isBackground(depth)) {
        hitInfo.emissive = vec3(0.0);
    }
}

vec3 viewPos;

vec3 colorForRay(in vec3 startRayPos, in vec3 startRayDir, inout StateType rngState) {
    vec3 ret = vec3(0.0, 0.0, 0.0);

    vec3 throughput = vec3(1.0, 1.0, 1.0);
    vec3 rayPos = startRayPos;
    vec3 rayDir = startRayDir;

    RayHitInfo hitInfo;
    RayHitInfo prevHitInfo;

    for (int bounceIndex = 0; bounceIndex <= dBounces; ++bounceIndex) {
        // shoot a ray out into the world
        if (bounceIndex == 0) {
            vec2 coords = gl_FragCoord.xy / uTexSize;
            float depth = getDepth(coords);

            hitInfo.missed = false;
            hitInfo.position = screenSpaceToViewSpace(vec3(coords, depth), uInvProjection);
            hitInfo.normal = -texture2D(tNormal, coords).rgb;
            hitInfo.color = texture2D(tShaded, coords).rgb;
            hitInfo.emissive = texture2D(tColor, coords).rgb * texture2D(tNormal, coords).a;

            // shadow
            #ifdef dShadowEnable
                #if dLightCount != 0
                    vec3 directLight = vec3(uAmbientColor);
                    #pragma unroll_loop_start
                    bool missed;
                    vec3 hitPos;
                    for (int i = 0; i < dLightCount; ++i) {
                        missed = false;
                        hitPos = viewPos + hitInfo.normal * RayPosNormalNudge;
                        hitPos += -uLightDirection[i] * (randomFloat(rngState));
                        rayMarch(-uLightDirection[i] + randomUnitVector(rngState) * uShadowSoftness, uShadowThickness, hitPos, missed);
                        if (missed) directLight += uLightColor[i];
                    }
                    #pragma unroll_loop_end
                    hitInfo.color *= directLight / uLightStrength;
                #endif
            #endif

            if (hitInfo.normal == vec3(0.0)) {
                hitInfo.missed = true;
            }
        } else {
            prevHitInfo = hitInfo;
            trace(rayPos, rayDir, hitInfo);
        }

        // if the ray missed, we are done
        if (hitInfo.missed) {
            vec3 accIrradiance = vec3(1.0);
            #ifdef dGlow
                if (bounceIndex > 1) {
                    accIrradiance = uLightStrength;
                }
            #else
                if (bounceIndex > 1) {
                    accIrradiance = uAmbientColor;
                    #if dLightCount != 0
                        #pragma unroll_loop_start
                        float dotNL;
                        vec3 irradiance;
                        for (int i = 0; i < dLightCount; ++i) {
                            dotNL = saturate(dot(prevHitInfo.normal, -uLightDirection[i]));
                            irradiance = dotNL * uLightColor[i];
                            accIrradiance += irradiance;
                        }
                        #pragma unroll_loop_end
                    #endif
                }
            #endif
            ret += prevHitInfo.color * accIrradiance * throughput;
            break;
        }

        // add emissive light
        ret += hitInfo.emissive * throughput;

        // update the ray position
        rayPos = hitInfo.position + hitInfo.normal * RayPosNormalNudge;

        // new ray direction from normal oriented cosine weighted hemisphere sample
        rayDir = normalize(hitInfo.normal + randomUnitVector(rngState));

        if (bounceIndex == 0) {
            continue;
        }

        // update the colorMultiplier.
        throughput *= hitInfo.color;

        // Russian Roulette
        // As the throughput gets smaller, the ray is more likely to get terminated early.
        // Survivors have their value boosted to make up for fewer samples being in the average.
        {
            float p = max(throughput.r, max(throughput.g, throughput.b));
            if (randomFloat(rngState) > p)
                break;

            // Add the energy we 'lose' by randomly terminating paths
            throughput *= 1.0 / p;
        }
    }

    // return pixel color
    return ret;
}

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    float depth = getDepth(coords);

    if (isBackground(depth)) {
        gl_FragColor = texture2D(tColor, coords);
        return;
    }

    #if __VERSION__ == 100
        float rngState = 26699.0;
    #else
        // initialize a random number state based on gl_FragCoord and uFrameNo
        uint rngState = uint(uint(gl_FragCoord.x) * 1973u + uint(gl_FragCoord.y) * 9277u + uint(uFrameNo) * 26699u) | 1u;
    #endif

    vec3 cameraPos = vec3(0.0, 0.0, 0.0);
    viewPos = screenSpaceToViewSpace(vec3(coords, depth), uInvProjection);
    vec3 rayDir = normalize(viewPos);

    // raytrace for this pixel
    vec3 color = vec3(0.0, 0.0, 0.0);
    for (int index = 0; index < int(dRendersPerFrame); ++index) {
        color += colorForRay(cameraPos, rayDir, rngState) / float(dRendersPerFrame);
    }

    // average the frames together
    vec4 lastFrameColor = texture2D(tAccumulate, coords);
    float blend = (uFrameNo < 1 || lastFrameColor.a == 0.0) ? 1.0 : 1.0 / (1.0 + (1.0 / lastFrameColor.a));
    color = mix(lastFrameColor.rgb, color, blend);

    // show the result
    gl_FragColor = vec4(color, blend);
}
`;