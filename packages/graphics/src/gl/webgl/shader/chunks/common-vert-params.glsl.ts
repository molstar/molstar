/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceDataVec4s } from '../../../renderable/merged.js';

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
    // (those stay declared, unread); 5 texels per instance: the transform's rows, (instance, segment groups), (segment
    // index). Rows of 1 << uInstanceDataShift instances, addressed by a shift and a mask (an integer division is
    // emulated); read by most chunks, once per vertex into globals (\`readInstanceData\`, first in main by the group
    // chunks)
    ivec2 instanceDataCoord;
    mat4 instanceTransform;
    vec4 instanceSegment;
    float instanceSegmentIndex;
    void readInstanceData() {
        int i = uInstanceBase + gl_InstanceID;
        instanceDataCoord = ivec2((i & ((1 << uInstanceDataShift) - 1)) * ${InstanceDataVec4s}, i >> uInstanceDataShift);
        instanceTransform = affineFromRows(
            texelFetch(tInstanceData, instanceDataCoord, 0),
            texelFetch(tInstanceData, instanceDataCoord + ivec2(1, 0), 0),
            texelFetch(tInstanceData, instanceDataCoord + ivec2(2, 0), 0)
        );
        instanceSegment = texelFetch(tInstanceData, instanceDataCoord + ivec2(3, 0), 0);
        instanceSegmentIndex = texelFetch(tInstanceData, instanceDataCoord + ivec2(4, 0), 0).x;
    }
    #define aTransform instanceTransform
    #define aInstance instanceSegment.x
    #define aSegment vec4(instanceSegment.yzw, instanceSegmentIndex)
#endif

#if defined(dSegmented)
    // the segment data of merged render items: the groups per instance, the sphere and level scales by the instance's
    // segment (\`aSegment.w\`), see SegmentSchema
    #define aSegmentSphere read_tSegment(int(aSegment.w) * 2)
    #define aSegmentLod read_tSegment(int(aSegment.w) * 2 + 1)
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
