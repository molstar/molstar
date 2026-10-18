/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * View space z of window depth values, depth in [0, 1] as the fragment stage's position.z and depth textures; depth is
 * reversed (near 1, far 0, see `DepthRangeTag`), which keeps the precision of small (far) depths
 */
export const depth_wgsl = /* wgsl */ `
fn perspectiveDepthToViewZ(depth: f32, near: f32, far: f32) -> f32 {
    return -(near * far) / (near + (far - near) * depth);
}

fn orthographicDepthToViewZ(depth: f32, near: f32, far: f32) -> f32 {
    return depth * (far - near) - far;
}

fn depthToViewZ(isOrtho: f32, depth: f32, near: f32, far: f32) -> f32 {
    if (isOrtho == 1.0) {
        return orthographicDepthToViewZ(depth, near, far);
    }
    return perspectiveDepthToViewZ(depth, near, far);
}
`;
