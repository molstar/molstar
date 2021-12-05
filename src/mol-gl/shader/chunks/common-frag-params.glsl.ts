export const common_frag_params = `
uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

#if dClipObjectCount != 0
    uniform int uClipObjectType[dClipObjectCount];
    uniform bool uClipObjectInvert[dClipObjectCount];
    uniform vec3 uClipObjectPosition[dClipObjectCount];
    uniform vec4 uClipObjectRotation[dClipObjectCount];
    uniform vec3 uClipObjectScale[dClipObjectCount];

    #if defined(dClipping)
        #if __VERSION__ == 100
            varying float vClipping;
        #else
            flat in float vClipping;
        #endif
    #endif
#endif

uniform vec3 uHighlightColor;
uniform vec3 uSelectColor;
uniform float uHighlightStrength;
uniform float uSelectStrength;
uniform int uMarkerPriority;

#if defined(dMarkerType_uniform)
    uniform float uMarker;
#elif defined(dMarkerType_groupInstance)
    #if __VERSION__ == 100
        varying float vMarker;
    #else
        flat in float vMarker;
    #endif
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

uniform vec2 uViewOffset;

uniform float uNear;
uniform float uFar;
uniform float uIsOrtho;

uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform bool uTransparentBackground;

uniform float uInteriorDarkening;
uniform bool uInteriorColorFlag;
uniform vec3 uInteriorColor;
bool interior;

uniform float uXrayEdgeFalloff;

uniform mat4 uProjection;

uniform bool uRenderWboit;
uniform bool uMarkingDepthTest;

uniform sampler2D tDepth;
uniform vec2 uDrawingBufferSize;

float getDepth(const in vec2 coords) {
    // always packed due to merged depth from primitives and volumes
    return unpackRGBAToDepth(texture2D(tDepth, coords));
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
`;