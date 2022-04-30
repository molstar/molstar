export const check_picking_alpha = `
float viewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
float fogFactor = smoothstep(uFogNear, uFogFar, abs(viewZ));
float alpha = (1.0 - fogFactor) * uAlpha;
// if not opaque enough ignore so the element below can be picked
if (uAlpha < uPickingAlphaThreshold || alpha < 0.1) {
    #ifdef dTransparentBackfaces_opaque
        if (!interior) discard;
    #else
        discard;
    #endif
}
`;