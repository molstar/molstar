/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SegmentTexels } from '../../../renderable/schema.js';

/**
 * Segment data of merged render items (`dSegmented`, see `SegmentSchema`): the group count and offsets per instance
 * (`aSegment.xyz`), the invariant bounding sphere (the object's otherwise) and the level scales by the instance's segment
 * (`aSegment.w`, `SegmentTexels` elements of `tSegment` per segment)
 */
export const segment_wgsl = /* wgsl */ `
fn invariantBoundingSphere(in: VertexInput) -> vec4<f32> {
    #if defined(dSegmented)
        return read_tSegment(i32(in.aSegment.w) * ${SegmentTexels});
    #else
        return uObject.uInvariantBoundingSphere;
    #endif
}

#if defined(dSegmented)
    /** the segment's group count, groupInstance offset and group offset */
    fn segmentGroups(in: VertexInput) -> vec3<i32> {
        return in.aSegment.xyz;
    }

    /** the segment's radius scale per lod level */
    fn segmentLod(in: VertexInput) -> vec4<f32> {
        return read_tSegment(i32(in.aSegment.w) * ${SegmentTexels} + 1);
    }

    /** the segment's stride per lod level (spheres', points': quad \`p\` draws item \`p * stride\`) */
    fn segmentStride(in: VertexInput) -> vec4<f32> {
        return read_tSegment(i32(in.aSegment.w) * ${SegmentTexels} + 2);
    }

    /** the segment's first item in the merged items (spheres', points'), its high and low 16 bits exact as floats */
    fn segmentItem(in: VertexInput) -> i32 {
        let b = read_tSegment(i32(in.aSegment.w) * ${SegmentTexels} + 3);
        return i32(b.x) * 65536 + i32(b.y);
    }

    /** the segment's first entry of the merged size table (spheres' \`siteType\` sizes) */
    fn segmentSizeBase(in: VertexInput) -> i32 {
        return i32(read_tSegment(i32(in.aSegment.w) * ${SegmentTexels} + 3).z);
    }

    /** the segment's first entry of the merged second color table (spheres' \`dSiteColor\` 'siteType') */
    fn segmentSiteColorBase(in: VertexInput) -> i32 {
        return i32(read_tSegment(i32(in.aSegment.w) * ${SegmentTexels} + 3).w);
    }

    /** the segment's quantized centers' origin (xyz) and scale (spheres' \`dSites\`) */
    fn segmentSites(in: VertexInput) -> vec4<f32> {
        return read_tSegment(i32(in.aSegment.w) * ${SegmentTexels} + 4);
    }
#endif
`;
