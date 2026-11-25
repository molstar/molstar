/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The overlay's occlusion (`dOverlay`): a fragment behind the scene's opaque depth (`tDepth`), or its front transparent
 * surface (`tDepthAlpha`) where that is opaque enough to be picked (`uPickingAlphaThreshold`), has its alpha multiplied
 * by `uOccludedAlpha`; in picking it is dropped where that alpha is below
 * the picking threshold, so that the scene keeps the pick. Behind by more than a thousandth of the distance: a clip
 * object's walls lie on the caps they cut, which would flip per pixel. Compared in view space, the overlay has near and far planes
 * of its own (`uOccluderNear`, `uOccluderFar` are the scene's). Expects `material` and `fragmentDepth`.
 */
export const apply_occlusion = `
#if defined(dOverlay) && (defined(dRenderVariant_color) || defined(dRenderVariant_pick) || defined(dRenderVariant_depth))
    if (uOccludedAlpha < 1.0) {
        vec2 occluderCoord = gl_FragCoord.xy / uDrawingBufferSize;
        float occluderDepth = getDepth(occluderCoord);
        float fragmentViewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
        float occluderViewZ = depthToViewZ(uIsOrtho, occluderDepth, uOccluderNear, uOccluderFar);
        // depth is reversed: 0 where the scene drew nothing
        float occlusionTolerance = 0.001 * abs(fragmentViewZ);
        bool occluded = occluderDepth > 0.0 && fragmentViewZ < occluderViewZ - occlusionTolerance;
        vec2 occluderFront = unpackRGBAToDepthWithAlpha(texture2D(tDepthAlpha, occluderCoord));
        if (occluderFront.y >= uPickingAlphaThreshold) {
            occluded = occluded || fragmentViewZ < depthToViewZ(uIsOrtho, occluderFront.x, uOccluderNear, uOccluderFar) - occlusionTolerance;
        }
        if (occluded) {
            #if defined(dRenderVariant_color)
                material.a *= uOccludedAlpha;
            #else
                if (uAlpha * uOccludedAlpha < uPickingAlphaThreshold) discard;
            #endif
        }
    }
#endif
`;
