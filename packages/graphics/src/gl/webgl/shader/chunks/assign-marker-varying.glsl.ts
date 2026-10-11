export const assign_marker_varying = `
#if defined(dNeedsMarker)
    #if defined(dMarkerType_instance)
        vMarker = readFromTexture(tMarker, aInstance, uMarkerTexDim).r;
    #elif defined(dMarkerType_groupInstance)
        vMarker = readFromTexture(tMarker, groupInstanceIndex(aInstance, group), uMarkerTexDim).r;
    #endif
#endif
`;
