export const apply_marking = `
float markingOpacity = clamp(material.a, 0.0, 1.0);
if (uMarkingType == 1) {
    if (marker > 0.0 || markingOpacity == 0.0)
        discard;
    material = packDepthWithAlphaToRGBA(fragmentDepth, markingOpacity);
} else {
    if (marker == 0.0)
        discard;
    float viewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
    float depthTest = 1.0;
    if (uMarkingDepthTest) {
        float markingDepth = unpackRGBAToDepthWithAlpha(texture2D(tDepthAlpha, gl_FragCoord.xy / uDrawingBufferSize)).x;
        // hidden behind it by more than a ten-thousandth of the distance and four steps of the packing's 24 bits (far
        // away a step is the larger): the depth may be the scene's, which holds the marked fragments themselves (rounded
        // by the packing, the sphere raster's arithmetic not the impostors')
        float markingViewZ = depthToViewZ(uIsOrtho, max(markingDepth - 4.0 / 16777216.0, 0.0), uNear, uFar);
        depthTest = markingDepth > 0.0 && viewZ < markingViewZ - 0.0001 * abs(viewZ) ? 1.0 : 0.0;
    }
    bool isHighlight = (uMarkerPriority == 1 && marker != 2.0) || (uMarkerPriority != 1 && marker == 1.0);
    float fogFactor = smoothstep(uFogNear, uFogFar, abs(viewZ));
    if (fogFactor == 1.0)
        discard;
    material = packMarkingMask(depthTest, isHighlight, 1.0 - fogFactor, markingOpacity);
}
`;
