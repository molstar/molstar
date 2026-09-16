export const apply_interior_color = `
if (interior) {
    #if defined(dInteriorColorType_uniform)
        #ifdef dUsePalette
            material.rgb = texture2D(tPalette, vec2(decodePaletteV(uInteriorThemeColor), 0.5)).rgb;
        #else
            material.rgb = uInteriorThemeColor;
        #endif
    #elif defined(dInteriorColorType_texture)
        #ifdef dUsePalette
            material.rgb = texture2D(tPalette, vec2(vInteriorPaletteV, 0.5)).rgb;
        #else
            material.rgb = vInteriorColor;
        #endif
    #endif
    material.rgb = mix(material.rgb, uInteriorColor.rgb, uInteriorColor.a);

    float isf = clamp(uInteriorSubstance.a, 0.0, 0.99); // clamp to avoid artifacts
    metalness = mix(metalness, uInteriorSubstance.r, isf);
    roughness = mix(roughness, uInteriorSubstance.g, isf);
    bumpiness = mix(bumpiness, uInteriorSubstance.b, isf);

    #ifdef dTransparentBackfaces_opaque
        material.a = 1.0;
    #endif
}
`;