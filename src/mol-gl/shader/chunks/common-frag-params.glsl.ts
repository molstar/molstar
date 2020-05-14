export default `
uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

#if dClipObjectCount != 0
    uniform int uClipObjectType[dClipObjectCount];
    uniform vec3 uClipObjectPosition[dClipObjectCount];
    uniform vec4 uClipObjectRotation[dClipObjectCount];
    uniform vec3 uClipObjectScale[dClipObjectCount];

    #if defined(dClipping)
        #if __VERSION__ != 300
            varying float vClipping;
        #else
            flat in float vClipping;
        #endif
    #endif
#endif

uniform vec3 uHighlightColor;
uniform vec3 uSelectColor;
#if __VERSION__ != 300
    varying float vMarker;
#else
    flat in float vMarker;
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

uniform vec2 uViewOffset;

uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform int uTransparentBackground;

uniform float uInteriorDarkening;
uniform int uInteriorColorFlag;
uniform vec3 uInteriorColor;
bool interior;
`;