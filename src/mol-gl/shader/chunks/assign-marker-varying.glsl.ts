export const assign_marker_varying = `
#if defined(dMarkerType_groupInstance)
    vMarker = readFromTexture(tMarker, aInstance * float(uGroupCount) + group, uMarkerTexDim).a;
#endif
`;