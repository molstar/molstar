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
`;