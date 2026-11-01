/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Levels of detail of the items drawn as quads (spheres, points), after `vModelPosition`: `lodFactor` the level's scale
 * (a coarser level's fewer items larger) faded at the level's range (1 without levels), `lodCulled` beyond the range
 */
export const assign_lod_factor = `
float lodFactor = 1.0;
bool lodCulled = false;
if (uLod.w != 0.0 && (uLod.x != 0.0 || uLod.y != 0.0)) {
    // members of a merged item can have different scales for the same lod level (item-count dependent)
    #if defined(dSegmented)
        float lodScale = uLodLevel == 0 ? aSegmentLod.x : (uLodLevel == 1 ? aSegmentLod.y : (uLodLevel == 2 ? aSegmentLod.z : aSegmentLod.w));
    #else
        float lodScale = uLod.w;
    #endif
    if (uModelScale != 1.0) {
        lodFactor = lodScale;
    } else {
        float d = (dot(uCameraPlane.xyz, vModelPosition) + uCameraPlane.w) / uModelScale;
        lodFactor = min(
            smoothstep(uLod.x, uLod.x + uLod.z, d),
            1.0 - smoothstep(uLod.y - uLod.z, uLod.y, d)
        ) * lodScale;
        lodCulled = d < uLod.x || d > uLod.y;
    }
}
`;
