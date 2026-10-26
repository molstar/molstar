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
    // sphere and level scales. Rows of 1 << uInstanceDataShift instances, addressed by a shift and a mask (an integer
    // division is emulated); the transform and (instance, segment), read by most chunks, once per vertex into globals
    // (\`readInstanceData\`, first in main by the group chunks), the sphere and level scales where read
    ivec2 instanceDataCoord;
    mat4 instanceTransform;
    vec4 instanceSegment;
    void readInstanceData() {
        int i = uInstanceBase + gl_InstanceID;
        instanceDataCoord = ivec2((i & ((1 << uInstanceDataShift) - 1)) * 7, i >> uInstanceDataShift);
        instanceTransform = mat4(
            texelFetch(tInstanceData, instanceDataCoord, 0),
            texelFetch(tInstanceData, instanceDataCoord + ivec2(1, 0), 0),
            texelFetch(tInstanceData, instanceDataCoord + ivec2(2, 0), 0),
            texelFetch(tInstanceData, instanceDataCoord + ivec2(3, 0), 0)
        );
        instanceSegment = texelFetch(tInstanceData, instanceDataCoord + ivec2(4, 0), 0);
    }
    #define aTransform instanceTransform
    #define aInstance instanceSegment.x
    #define aSegment instanceSegment.yzw
    #define aSegmentSphere texelFetch(tInstanceData, instanceDataCoord + ivec2(5, 0), 0)
    #define aSegmentLod texelFetch(tInstanceData, instanceDataCoord + ivec2(6, 0), 0)
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
