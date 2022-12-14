export const clip_pixel = `
#if defined(dClipVariant_pixel) && dClipObjectCount != 0
    if (clipTest(vec4(vModelPosition, 0.0)))
        discard;
#endif
`;