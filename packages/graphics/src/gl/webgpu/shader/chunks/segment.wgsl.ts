/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Per-instance data of merged render items (`dSegmented`, see `SegmentSchema`): the invariant bounding sphere of the
 * instance's segment, the object's otherwise
 */
export const segment_wgsl = /* wgsl */ `
fn invariantBoundingSphere(in: VertexInput) -> vec4<f32> {
    #if defined(dSegmented)
        return in.aSegmentSphere;
    #else
        return uObject.uInvariantBoundingSphere;
    #endif
}
`;
