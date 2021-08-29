export const common_vert_params = `
uniform mat4 uProjection, uModel, uView;
uniform vec3 uCameraPosition;

uniform int uObjectId;
uniform int uVertexCount;
uniform int uInstanceCount;
uniform int uGroupCount;
uniform vec4 uInvariantBoundingSphere;

#if dClipObjectCount != 0
    uniform int uClipObjectType[dClipObjectCount];
    uniform bool uClipObjectInvert[dClipObjectCount];
    uniform vec3 uClipObjectPosition[dClipObjectCount];
    uniform vec4 uClipObjectRotation[dClipObjectCount];
    uniform vec3 uClipObjectScale[dClipObjectCount];

    #if defined(dClipping)
        uniform vec2 uClippingTexDim;
        uniform sampler2D tClipping;
        #if __VERSION__ == 100
            varying float vClipping;
        #else
            flat out float vClipping;
        #endif
    #endif
#endif

#if defined(dMarkerType_uniform)
    uniform float uMarker;
#elif defined(dMarkerType_groupInstance)
    uniform vec2 uMarkerTexDim;
    uniform sampler2D tMarker;
    #if __VERSION__ == 100
        varying float vMarker;
    #else
        flat out float vMarker;
    #endif
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

#if __VERSION__ == 100
    attribute float aVertex;
    #define VertexID int(aVertex)
#else
    // not using gl_VertexID but aVertex to ensure there is an active attribute with divisor 0
    // since FF 85 this is not needed anymore but lets keep it for backwards compatibility
    // https://bugzilla.mozilla.org/show_bug.cgi?id=1679693
    // see also note in src/mol-gl/webgl/render-item.ts
    attribute float aVertex;
    #define VertexID int(aVertex)
    // #define VertexID gl_VertexID
#endif
`;