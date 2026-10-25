/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const common_vert_params = `
#if defined(dClipObjects)
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

#if defined(dInstanceData)
    // the instance data of merged render items on WebGL2 without base instances (see \`getInstanceDataSchema\`): read
    // from a texture at the draw's first instance plus gl_InstanceID, instead of instanced attributes re-pointed per draw
    // (those stay declared, unread); 7 texels per instance: the transform's columns, (instance, segment), the segment's
    // sphere and level scales
    vec4 instanceTexel(const in int k) {
        int t = (uInstanceBase + gl_InstanceID) * 7 + k;
        int w = textureSize(tInstanceData, 0).x;
        return texelFetch(tInstanceData, ivec2(t % w, t / w), 0);
    }
    #define aTransform mat4(instanceTexel(0), instanceTexel(1), instanceTexel(2), instanceTexel(3))
    #define aInstance instanceTexel(4).x
    #define aSegment instanceTexel(4).yzw
    #define aSegmentSphere instanceTexel(5)
    #define aSegmentLod instanceTexel(6)
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
