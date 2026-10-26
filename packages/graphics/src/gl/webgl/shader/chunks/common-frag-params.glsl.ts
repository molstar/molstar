/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

export const common_frag_params = `
#if defined(dClipObjects)
    #if defined(dClipping)
        #if __VERSION__ == 100 || !defined(dVaryingGroup)
            varying float vClipping;
        #else
            flat in float vClipping;
        #endif
    #endif
#endif

#if defined(dNeedsMarker)
    #if __VERSION__ == 100 || !defined(dVaryingGroup)
        varying float vMarker;
    #else
        flat in float vMarker;
    #endif
#endif

#if defined(dRenderVariant_colorDpoit)
    #define MAX_DPOIT_DEPTH 99999.0 // NOTE constant also set in TypeScript
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

bool interior;

float getDepthPacked(const in vec2 coords) {
    return unpackRGBAToDepth(texture2D(tDepth, coords));
}

float getDepth(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepth, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #endif
}

float calcDepth(const in vec3 pos) {
    vec2 clipZW = pos.z * uProjection[2].zw + uProjection[3].zw;
    return clipZToDepth(clipZW.x / clipZW.y);
}

vec3 perturbNormal(in vec3 modelPosition, in vec3 viewPosition, in vec3 normal, in float frequency, in float scale) {
    vec3 t1 = normalize(cross(normal, abs(normal.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0)));
    vec3 t2 = cross(normal, t1);
    float e = max((uIsOrtho == 1.0 ? 2.0 : 2.0 * abs(viewPosition.z)) / (uProjection[1][1] * uDrawingBufferSize.y), 0.01 / frequency);
    vec3 p = modelPosition * frequency;
    float h0 = fbm(p);
    float h1 = fbm(p + (mat3(uInvView) * t1) * (e * frequency));
    float h2 = fbm(p + (mat3(uInvView) * t2) * (e * frequency));
    return normalize(normal - scale * ((h1 - h0) * t1 + (h2 - h0) * t2) / e);
}

#ifdef dXrayShaded
    float calcXrayShadedAlpha(in float alpha, const in vec3 normal) {
        float d = abs(dot(normal, vec3(0.0, 0.0, 1.0)));
        // d^falloff, at the silhouette (d = 0) its limit: pow(0, y) is undefined for y = 0 (the falloff's minimum)
        float p = d > 0.0 ? pow(d, uXrayEdgeFalloff) : (uXrayEdgeFalloff == 0.0 ? 1.0 : 0.0);
        #if defined(dXrayShaded_on)
            alpha *= 1.0 - p;
        #elif defined(dXrayShaded_inverted)
            alpha *= p;
        #endif
        return clamp(alpha, 0.001, 0.999);
    }
#endif
`;
