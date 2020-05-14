export default `
#if defined(dClipVariant_pixel) && dClipObjectCount != 0
    int flag = 0;
    #if defined(dClipping)
        flag = int(floor(vClipping * 255.0 + 0.5));
    #endif

    if (clipTest(vec4(vModelPosition, 0.0), flag))
        discard;
#endif
`;