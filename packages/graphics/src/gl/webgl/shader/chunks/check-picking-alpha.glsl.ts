/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const check_picking_alpha = `
float viewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
float fogFactor = smoothstep(uFogNear, uFogFar, abs(viewZ));
float fogAlpha = (1.0 - fogFactor) * uAlpha;
// widgets are picked at any opacity (uPickingAlphaMin 1)
float alpha = max(uAlpha, uPickingAlphaMin);
#ifdef dXrayShaded
    // add bias to make picking xray shaded elements easier
    alpha = calcXrayShadedAlpha(alpha, normal) + (0.3 * uPickingAlphaThreshold);
#endif
// if not opaque enough ignore so the element below can be picked
if (alpha < uPickingAlphaThreshold || fogAlpha < 0.1) {
    #ifdef dTransparentBackfaces_opaque
        if (!interior) discard;
    #else
        discard;
    #endif
}
`;
