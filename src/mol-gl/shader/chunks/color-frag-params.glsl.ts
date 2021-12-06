export const color_frag_params = `
uniform float uMetalness;
uniform float uRoughness;
uniform float uBumpiness;
#ifdef bumpEnabled
    uniform float uBumpFrequency;
    uniform float uBumpAmplitude;
#endif

#if defined(dRenderVariant_color)
    #if defined(dColorType_uniform)
        uniform vec3 uColor;
    #elif defined(dColorType_varying)
        varying vec4 vColor;
    #endif

    #ifdef dOverpaint
        varying vec4 vOverpaint;
    #endif

    #ifdef dSubstance
        varying vec4 vSubstance;
    #endif
#elif defined(dRenderVariant_pick)
    #if __VERSION__ == 100
        varying vec4 vColor;
    #else
        flat in vec4 vColor;
    #endif
#endif

#ifdef dTransparency
    varying float vGroup;
    varying float vTransparency;
#endif

#ifdef dUsePalette
    uniform sampler2D tPalette;
    varying float vPaletteV;
#endif
`;