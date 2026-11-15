/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceDataVec4s } from '../../../renderable/merged.js';
import { SegmentTexels } from '../../../renderable/schema.js';

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

// the position where the instance is placed, without the transient transform: what grids of 'volumeInstance' data are
// looked up by (built from \`aTransform\`); set by \`assign_position\`
vec3 gridModelPosition;

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
    // (those stay declared, unread); 4 texels per instance: the transform's rows, (instance, segment index); the
    // segments' groups in a slot each after the instances (\`uInstanceDataSegments\`); in pages of instances, a texture
    // each, bound per draw (\`uInstanceBase\` the page's). Rows of 1 << uInstanceDataShift
    // slots, addressed by a shift and a mask (an integer division is emulated); read by most chunks, once per vertex
    // into globals (\`readInstanceData\`, first in main by the group chunks)
    mat4 instanceTransform;
    uint instanceId;
    ivec4 instanceSegment;
    ivec2 instanceDataSlot(int i) {
        return ivec2((i & ((1 << uInstanceDataShift) - 1)) * ${InstanceDataVec4s}, i >> uInstanceDataShift);
    }
    void readInstanceData() {
        ivec2 c = instanceDataSlot(uInstanceBase + gl_InstanceID);
        instanceTransform = affineFromRows(
            texelFetch(tInstanceData, c, 0),
            texelFetch(tInstanceData, c + ivec2(1, 0), 0),
            texelFetch(tInstanceData, c + ivec2(2, 0), 0)
        );
        vec4 ids = texelFetch(tInstanceData, c + ivec2(3, 0), 0);
        instanceId = floatBitsToUint(ids.x);
        int segment = floatBitsToInt(ids.y);
        vec4 groups = texelFetch(tInstanceData, instanceDataSlot(uInstanceDataSegments + segment), 0);
        instanceSegment = ivec4(floatBitsToInt(groups.xyz), segment);
    }
    #define aTransform instanceTransform
    // the integers as their bits (\`packInstanceData\`)
    #define aInstance instanceId
    #define aSegment instanceSegment
#endif

#if defined(dSegmented)
    // the segment data of merged render items: the groups per instance, the sphere, level scales and strides and first item
    // (spheres', points') by the instance's
    // segment (\`aSegment.w\`), see SegmentSchema
    #define aSegmentSphere read_tSegment(int(aSegment.w) * ${SegmentTexels})
    #define aSegmentLod read_tSegment(int(aSegment.w) * ${SegmentTexels} + 1)
    #define aSegmentStride read_tSegment(int(aSegment.w) * ${SegmentTexels} + 2)
    #define aSegmentItem read_tSegment(int(aSegment.w) * ${SegmentTexels} + 3)
    #define aSegmentSites read_tSegment(int(aSegment.w) * ${SegmentTexels} + 4)
    #define groupInstanceIndex(inst, grp) (int(aSegment.y) + int(inst) * int(aSegment.x) + int(grp))
    #define groupIndex(grp) (int(aSegment.z) + int(grp))
    #define InvariantBoundingSphere aSegmentSphere
#else
    #define groupInstanceIndex(inst, grp) (int(inst) * uGroupCount + int(grp))
    #define groupIndex(grp) int(grp)
    #define InvariantBoundingSphere uInvariantBoundingSphere
#endif
// the index of the marker and layer data, per instance or group instance (uInstanceGranularity); integers: the instance
// ids (\`aInstance\`) and offsets past 2^24 exact (WebGL2; floats on WebGL1)
#define granularityIndex(inst, grp) (uInstanceGranularity ? int(inst) : groupInstanceIndex(inst, grp))
`;
