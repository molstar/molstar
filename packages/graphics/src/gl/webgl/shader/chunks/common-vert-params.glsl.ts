export const common_vert_params = `
#if dClipObjectCount != 0
    #if defined(dClipping)
        #if __VERSION__ == 100 || !defined(dVaryingGroup)
            varying float vClipping;
        #else
            flat out float vClipping;
        #endif
    #endif
#endif

#if defined(dNeedsMarker)
    #if __VERSION__ == 100 || !defined(dVaryingGroup)
        varying float vMarker;
    #else
        flat out float vMarker;
    #endif
#endif

varying vec3 vModelPosition;
varying vec3 vViewPosition;

#if defined(emulateVertexId)
    #define VertexID int(aVertex)
#else
    // int() is needed for some Safari versions
    // see https://bugs.webkit.org/show_bug.cgi?id=244152
    #define VertexID int(gl_VertexID)
#endif

#if defined(enabledMultiDraw)
    #define DrawID gl_DrawID
#else
    #define DrawID uDrawId
#endif

#if defined(dSegmented)
    // per-instance segment data for merged render items, see SegmentSchema
    #define groupInstanceIndex(inst, grp) (aSegment.y + (inst) * aSegment.x + (grp))
    #define groupIndex(grp) (aSegment.z + (grp))
    #define InvariantBoundingSphere aSegmentSphere
#else
    #define groupInstanceIndex(inst, grp) ((inst) * float(uGroupCount) + (grp))
    #define groupIndex(grp) (grp)
    #define InvariantBoundingSphere uInvariantBoundingSphere
#endif
// the index of the marker and layer data, per instance or group instance (uInstanceGranularity)
#define granularityIndex(inst, grp) (uInstanceGranularity ? (inst) : groupInstanceIndex(inst, grp))
`;
