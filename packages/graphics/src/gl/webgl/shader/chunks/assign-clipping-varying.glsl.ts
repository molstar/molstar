export const assign_clipping_varying = `
#if dClipObjectCount != 0 && defined(dClipping)
    #if defined(dClippingType_instance)
        vClipping = readFromTexture(tClipping, aInstance, uClippingTexDim).r;
    #elif defined(dMarkerType_groupInstance)
        vClipping = readFromTexture(tClipping, groupInstanceIndex(aInstance, group), uClippingTexDim).r;
    #endif
#endif
`;
