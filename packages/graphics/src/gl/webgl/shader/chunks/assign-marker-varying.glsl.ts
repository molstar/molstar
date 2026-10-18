export const assign_marker_varying = `
#if defined(dNeedsMarker)
    vMarker = readFromTexture(tMarker, granularityIndex(aInstance, group), uMarkerTexDim).r;
#endif
`;
