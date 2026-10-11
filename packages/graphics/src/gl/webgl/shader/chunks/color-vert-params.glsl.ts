export const color_vert_params = `
#if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
    #if defined(dColorType_texture) || defined(dColorType_grid) || defined(dColorType_direct)
        varying vec4 vColor;
    #endif

    #ifdef dUsePalette
        varying float vPaletteV;
    #endif

    #ifdef dOverpaint
        #if defined(dOverpaintType_instance) || defined(dOverpaintType_groupInstance) || defined(dOverpaintType_vertexInstance) || defined(dOverpaintType_volumeInstance)
            varying vec4 vOverpaint;
        #endif
    #endif

    #ifdef dEmissive
        #if defined(dEmissiveType_instance) || defined(dEmissiveType_groupInstance) || defined(dEmissiveType_vertexInstance) || defined(dEmissiveType_volumeInstance)
            varying float vEmissive;
        #endif
    #endif

    #ifdef dSubstance
        #if defined(dSubstanceType_instance) || defined(dSubstanceType_groupInstance) || defined(dSubstanceType_vertexInstance) || defined(dSubstanceType_volumeInstance)
            varying vec4 vSubstance;
        #endif
    #endif
#elif defined(dRenderVariant_emissive)
    #ifdef dEmissive
        #if defined(dEmissiveType_instance) || defined(dEmissiveType_groupInstance) || defined(dEmissiveType_vertexInstance) || defined(dEmissiveType_volumeInstance)
            varying float vEmissive;
        #endif
    #endif
#elif defined(dRenderVariant_pick)
    #if __VERSION__ == 100 || !defined(dVaryingGroup)
        #ifdef requiredDrawBuffers
            varying vec4 vObject;
            varying vec4 vInstance;
            varying vec4 vGroup;
        #else
            varying vec4 vColor;
        #endif
    #else
        #ifdef requiredDrawBuffers
            flat out vec4 vObject;
            flat out vec4 vInstance;
            flat out vec4 vGroup;
        #else
            flat out vec4 vColor;
        #endif
    #endif
#endif

#ifdef dTransparency
    #if defined(dTransparencyType_instance) || defined(dTransparencyType_groupInstance) || defined(dTransparencyType_vertexInstance) || defined(dTransparencyType_volumeInstance)
        varying float vTransparency;
    #endif
#endif
`;
