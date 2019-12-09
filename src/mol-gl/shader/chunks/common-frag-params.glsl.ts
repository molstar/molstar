export default `
uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

uniform vec3 uHighlightColor;
uniform vec3 uSelectColor;
#if __VERSION__ != 300
    varying float vMarker;
#else
    flat in float vMarker;
#endif

varying vec3 vViewPosition;

uniform vec2 uViewOffset;

uniform int uFogFlag;
uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform int uPickable;
uniform int uTransparentBackground;

uniform float uInteriorDarkening;
uniform int uInteriorColorFlag;
uniform vec3 uInteriorColor;
bool interior;
`