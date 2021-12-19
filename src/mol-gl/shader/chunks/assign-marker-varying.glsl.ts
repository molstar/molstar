export const assign_marker_varying = `
#if defined(dRenderVariant_color) || defined(dRenderVariant_marking)
    vMarker = readFromTexture(tMarker, aInstance * float(uGroupCount) + group, uMarkerTexDim).a;
#endif
`;