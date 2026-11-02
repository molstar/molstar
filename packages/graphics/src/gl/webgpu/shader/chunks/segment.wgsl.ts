/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Segment data of merged render items (`dSegmented`, see `SegmentSchema`): the group count and offsets per instance
 * (`aSegment.xyz`), the invariant bounding sphere (the object's otherwise) and the level scales by the instance's segment
 * (`aSegment.w`, two elements of `tSegment` per segment)
 */
export const segment_wgsl = /* wgsl */ `
fn invariantBoundingSphere(in: VertexInput) -> vec4<f32> {
    #if defined(dSegmented)
        return read_tSegment(i32(in.aSegment.w) * 2);
    #else
        return uObject.uInvariantBoundingSphere;
    #endif
}

#if defined(dSegmented)
    /** the segment's group count, groupInstance offset and group offset */
    fn segmentGroups(in: VertexInput) -> vec3<f32> {
        return in.aSegment.xyz;
    }

    /** the segment's radius scale per lod level */
    fn segmentLod(in: VertexInput) -> vec4<f32> {
        return read_tSegment(i32(in.aSegment.w) * 2 + 1);
    }
#endif
`;
