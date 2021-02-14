export const check_picking_alpha = `
float viewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
float fogFactor = smoothstep(uFogNear, uFogFar, abs(viewZ));
float alpha = (1.0 - fogFactor) * uAlpha;
if (uAlpha < uPickingAlphaThreshold || alpha < 0.1)
    discard; // ignore so the element below can be picked
`;