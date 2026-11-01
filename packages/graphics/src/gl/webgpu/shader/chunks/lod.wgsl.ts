/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Levels of detail of the items drawn as quads (spheres, points), free of any binding: shared by the render stages
 * (`assign_lod_factor_wgsl`) and compute kernels that draw the same items
 */
export const lod_at_wgsl = /* wgsl */ `
struct LodAt {
    /** the level's scale faded within the overlap at the ends of its range */
    factor: f32,
    /** beyond the range */
    culled: bool,
};

/** The level of range \`lod\` (min, max, overlap) and scale \`scale\` at distance \`d\` */
fn lodAt(d: f32, lod: vec3<f32>, scale: f32) -> LodAt {
    let factor = min(smoothstep(lod.x, lod.x + lod.z, d), 1.0 - smoothstep(lod.y - lod.z, lod.y, d)) * scale;
    return LodAt(factor, d < lod.x || d > lod.y);
}
`;
