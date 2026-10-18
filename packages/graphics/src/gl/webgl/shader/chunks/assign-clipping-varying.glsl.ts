export const assign_clipping_varying = `
#if dClipObjectCount != 0 && defined(dClipping)
    vClipping = readFromTexture(tClipping, granularityIndex(aInstance, group), uClippingTexDim).r;
#endif
`;
