export default `
#if dClipObjectCount != 0 && defined(dClipping)
    vClipping = readFromTexture(tClipping, aInstance * float(uGroupCount) + group, uClippingTexDim).a;
#endif
`;