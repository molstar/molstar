export default `
float depth = length(vViewPosition);
float fogFactor = smoothstep(uFogNear, uFogFar, depth);
float alpha = (1.0 - fogFactor) * uAlpha;
if (uAlpha < uPickingAlphaThreshold || alpha < 0.1)
    discard; // ignore so the element below can be picked
`;