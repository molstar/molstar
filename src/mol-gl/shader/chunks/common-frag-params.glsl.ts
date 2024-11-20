export const common_frag_params = `
uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

uniform int uPickType;
uniform int uMarkingType;

uniform vec4 uCameraPlane;
uniform vec4 uLod;

#if dClipObjectCount != 0
    uniform int uClipObjectType[dClipObjectCount];
    uniform bool uClipObjectInvert[dClipObjectCount];
    uniform vec3 uClipObjectPosition[dClipObjectCount];
    uniform vec4 uClipObjectRotation[dClipObjectCount];
    uniform vec3 uClipObjectScale[dClipObjectCount];

    #if defined(dClipping)
        #if __VERSION__ == 100 || defined(dClippingType_instance) || !defined(dVaryingGroup)
            varying float vClipping;
        #else
            flat in float vClipping;
        #endif
    #endif
#endif

#if defined(dColorMarker)
    uniform vec3 uHighlightColor;
    uniform vec3 uSelectColor;
    uniform vec3 uDimColor;
    uniform float uHighlightStrength;
    uniform float uSelectStrength;
    uniform float uDimStrength;
    uniform int uMarkerPriority;
    uniform float uMarkerAverage;
#endif

#if defined(dNeedsMarker)
    uniform float uMarker;
    #if __VERSION__ == 100 || defined(dMarkerType_instance) || !defined(dVaryingGroup)
        varying float vMarker;
    #else
        flat in float vMarker;
    #endif
#endif

#if defined(dRenderVariant_colorDpoit)
    #define MAX_DPOIT_DEPTH 99999.0 // NOTE constant also set in TypeScript
    uniform sampler2D tDpoitDepth;
    uniform sampler2D tDpoitFrontColor;
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

uniform vec2 uViewOffset;

uniform float uNear;
uniform float uFar;
uniform float uIsOrtho;

uniform bool uFog;
uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform bool uTransparentBackground;

uniform bool uDoubleSided;
uniform float uInteriorDarkening;
uniform bool uInteriorColorFlag;
uniform vec3 uInteriorColor;
bool interior;

uniform float uXrayEdgeFalloff;
uniform float uCelSteps;
uniform float uExposure;

uniform mat4 uProjection;

uniform int uRenderMask;
uniform bool uMarkingDepthTest;

uniform sampler2D tDepth;
uniform vec2 uDrawingBufferSize;

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
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}

// "Bump Mapping Unparametrized Surfaces on the GPU" Morten S. Mikkelsen
// https://mmikk.github.io/papers3d/mm_sfgrad_bump.pdf
vec3 perturbNormal(in vec3 position, in vec3 normal, in float height, in float scale) {
    vec3 sigmaS = dFdx(position);
    vec3 sigmaT = dFdy(position);

    vec3 r1 = cross(sigmaT, normal);
    vec3 r2 = cross(normal, sigmaS);
    float det = dot(sigmaS, r1);
    if (det == 0.0) return normal;

    float bs = dFdx(height);
    float bt = dFdy(height);

    vec3 surfGrad = sign(det) * (bs * r1 + bt * r2);
    return normalize(abs(det) * normal - scale * surfGrad);
}

float hash(in float h) {
    return fract(sin(h) * 43758.5453123);
}

float noise(in vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);

    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(
        mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
            mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
        mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
            mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

float fbm(in vec3 p) {
    float f = 0.0;
    f += 0.5 * noise(p);
    p *= 2.01;
    f += 0.25 * noise(p);
    p *= 2.02;
    f += 0.125 * noise(p);

    return f;
}

#ifdef dXrayShaded
    float calcXrayShadedAlpha(in float alpha, const in vec3 normal) {
        #if defined(dXrayShaded_on)
            alpha *= 1.0 - pow(abs(dot(normal, vec3(0.0, 0.0, 1.0))), uXrayEdgeFalloff);
        #elif defined(dXrayShaded_inverted)
            alpha *= pow(abs(dot(normal, vec3(0.0, 0.0, 1.0))), uXrayEdgeFalloff);
        #endif
        return clamp(alpha, 0.001, 0.999);
    }
#endif
`;