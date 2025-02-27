export const clip_pixel = `
#if defined(dClipVariant_pixel) && dClipObjectCount != 0
    if (clipTest(vModelPosition))
        discard;
#endif
`;