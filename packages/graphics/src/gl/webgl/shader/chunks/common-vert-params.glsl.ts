/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceDataVec4s } from '../../../renderable/merged.js';
import { SegmentTexels } from '../../../renderable/schema.js';
import {
  CopyLayerBlock,
  CopyLayerKinds,
  CopyMarkerTexels,
  CopyOrderPerRow,
  CopyPosesPerRow,
  CopyPoseTexels,
} from '../../../instance-expansion.js';

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
    #if defined(dInstanceExpansion)
        // an expanded object's instance (\`instance-expansion.ts\`): its object instance \`k\` placed by the pose \`v\`
        // of the list culled on the CPU (the draw's copies in \`tCopyOrder\` from \`uCopyBase\`, their poses in
        // \`tCopyPoses\`, ${CopyPoseTexels} texels each: the affine rows, ${CopyPosesPerRow} a row); its draw slot what
        // it is picked by (its pose's and its own instance's: \`v n + aInstance\`)
        uint expansionSlot;
        // its copy's marker (\`tCopyMarkers\`: two words a copy from \`uCopyMarkerBase\`, its own marker or its block's
        // offset, then the blocks' bytes by the template particle's index, ${CopyMarkerTexels} bytes a row), or-ed with
        // its template particle's (\`assign_marker_varying\`)
        uint expansionCopyMarker;
        // its copy (\`uCopyMarkerBase\` and on) and template particle: its copy's own layers (\`copyLayerOr\`)
        int expansionCopy;
        uint expansionLocal;
        uvec4 copyMarkerTexel(int i) {
            return uvec4(texelFetch(tCopyMarkers, ivec2(i % ${CopyMarkerTexels}, i / ${CopyMarkerTexels}), 0) * 255.0 + 0.5);
        }
        uint copyMarker(int copy, uint local) {
            uvec4 w = copyMarkerTexel(2 * copy);
            // \`CopyMarkerBlock\`: the word's top bit
            if (w.a < 128u) return w.r;
            uvec4 o = copyMarkerTexel(2 * copy + 1);
            int b = int(o.r | (o.g << 8u) | (o.b << 16u) | (o.a << 24u)) + int(local);
            return copyMarkerTexel(b / 4)[b % 4];
        }
    #endif
    ivec2 instanceDataSlot(int i) {
        return ivec2((i & ((1 << uInstanceDataShift) - 1)) * ${InstanceDataVec4s}, i >> uInstanceDataShift);
    }
    void readInstanceData() {
        #if defined(dInstanceExpansion)
            int k = gl_InstanceID % uExpansionCount;
            int v = uCopyBase + gl_InstanceID / uExpansionCount;
            ivec2 c = instanceDataSlot(uInstanceBase + k);
        #else
            ivec2 c = instanceDataSlot(uInstanceBase + gl_InstanceID);
        #endif
        instanceTransform = affineFromRows(
            texelFetch(tInstanceData, c, 0),
            texelFetch(tInstanceData, c + ivec2(1, 0), 0),
            texelFetch(tInstanceData, c + ivec2(2, 0), 0)
        );
        vec4 ids = texelFetch(tInstanceData, c + ivec2(3, 0), 0);
        instanceId = floatBitsToUint(ids.x);
        #if defined(dInstanceExpansion)
            int t = v / 4;
            int copy = int(texelFetch(tCopyOrder, ivec2(t % ${CopyOrderPerRow / 4}, t / ${CopyOrderPerRow / 4}), 0)[v % 4]);
            ivec2 p = ivec2((copy % ${CopyPosesPerRow}) * ${CopyPoseTexels}, copy / ${CopyPosesPerRow});
            instanceTransform = affineFromRows(
                texelFetch(tCopyPoses, p, 0),
                texelFetch(tCopyPoses, p + ivec2(1, 0), 0),
                texelFetch(tCopyPoses, p + ivec2(2, 0), 0)
            ) * instanceTransform;
            expansionSlot = uint(v) * uint(uExpansionCount) + instanceId;
            expansionCopy = uCopyMarkerBase + copy;
            expansionLocal = floatBitsToUint(ids.z);
            expansionCopyMarker = copyMarker(expansionCopy, expansionLocal);
        #endif
        #if defined(dSegmented)
            int segment = floatBitsToInt(ids.y);
            vec4 groups = texelFetch(tInstanceData, instanceDataSlot(uInstanceDataSegments + segment), 0);
            instanceSegment = ivec4(floatBitsToInt(groups.xyz), segment);
        #endif
    }
    #if defined(dInstanceExpansion) && defined(dCopyLayers)
        // layer \`k\` of the instance's copy (\`tCopyLayers\`: ${CopyLayerKinds} words a copy, none, a value's first byte
        // or a block's, by template particle, \`size\` bytes an item), else \`own\`, its object's
        vec4 copyLayerOr(int k, int size, vec4 own) {
            int i = expansionCopy * ${CopyLayerKinds} + k;
            uvec4 w = uvec4(texelFetch(tCopyLayers, ivec2(i % ${CopyMarkerTexels}, i / ${CopyMarkerTexels}), 0) * 255.0 + 0.5);
            uint word = w.r | (w.g << 8u) | (w.b << 16u) | (w.a << 24u);
            if (word == 0u) return own;
            int b = int(word & 0x3fffffffu);
            if ((word & ${CopyLayerBlock}u) != 0u) b += int(expansionLocal) * size;
            int t = b / 4;
            vec4 texel = texelFetch(tCopyLayers, ivec2(t % ${CopyMarkerTexels}, t / ${CopyMarkerTexels}), 0);
            return size == 4 ? texel : vec4(texel[b % 4]);
        }
        #define copyLayerOr4(k, own) copyLayerOr(k, 4, own)
        #define copyLayerOr1(k, own) copyLayerOr(k, 1, vec4(own)).r
    #endif
    #define aTransform instanceTransform
    // the integers as their bits (\`packInstanceData\`)
    #define aInstance instanceId
    #define aSegment instanceSegment
#endif

#if !defined(copyLayerOr4)
    // an object's own layers alone: not drawn by copies with layers of their own
    #define copyLayerOr4(k, own) own
    #define copyLayerOr1(k, own) own
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
