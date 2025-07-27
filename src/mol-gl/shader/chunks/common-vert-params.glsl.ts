export const common_vert_params = `
uniform mat4 uProjection, uModel, uView;
uniform vec3 uCameraPosition;
uniform vec4 uCameraPlane;

uniform int uObjectId;
uniform int uVertexCount;
uniform int uInstanceCount;
uniform int uGroupCount;
uniform vec4 uInvariantBoundingSphere;
uniform vec4 uLod;

uniform bool uDoubleSided;
uniform int uPickType;

#if dClipObjectCount != 0
    uniform int uClipObjectType[dClipObjectCount];
    uniform bool uClipObjectInvert[dClipObjectCount];
    uniform vec3 uClipObjectPosition[dClipObjectCount];
    uniform vec4 uClipObjectRotation[dClipObjectCount];
    uniform vec3 uClipObjectScale[dClipObjectCount];
    uniform mat4 uClipObjectTransform[dClipObjectCount];

    #if defined(dClipping)
        uniform vec2 uClippingTexDim;
        uniform sampler2D tClipping;
        #if __VERSION__ == 100 || defined(dClippingType_instance) || !defined(dVaryingGroup)
            varying float vClipping;
        #else
            flat out float vClipping;
        #endif
    #endif
#endif

#if defined(dNeedsMarker)
    uniform float uMarker;
    uniform vec2 uMarkerTexDim;
    uniform sampler2D tMarker;
    #if __VERSION__ == 100 || defined(dMarkerType_instance) || !defined(dVaryingGroup)
        varying float vMarker;
    #else
        flat out float vMarker;
    #endif
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

uniform float uModelScale;

#if defined(noNonInstancedActiveAttribs)
    // int() is needed for some Safari versions
    // see https://bugs.webkit.org/show_bug.cgi?id=244152
    #define VertexID int(gl_VertexID)
#else
    attribute float aVertex;
    #define VertexID int(aVertex)
#endif

#if defined(enabledMultiDraw)
    #define DrawID gl_DrawID
#else
    #define DrawID uDrawId
#endif
`;