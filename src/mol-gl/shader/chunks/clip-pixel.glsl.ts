export default `
#if defined(dClipVariant_pixel) && dClipObjectCount != 0
    #if defined(dClipping)
        int clippingFlag = int(floor(vClipping * 255.0 + 0.5));
    #else
        int clippingFlag = 0;
    #endif

    if (clipTest(vec4(vModelPosition, 0.0), clippingFlag))
        discard;
#endif
`;