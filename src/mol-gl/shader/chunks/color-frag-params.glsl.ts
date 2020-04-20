export default `
#if defined(dRenderVariant_color)
    #if defined(dColorType_uniform)
        uniform vec3 uColor;
    #elif defined(dColorType_varying)
        varying vec4 vColor;
    #endif

    #ifdef dOverpaint
        varying vec4 vOverpaint;
    #endif

    #ifdef dTransparency
        varying float vGroup;
        varying float vTransparency;
    #endif
#elif defined(dRenderVariant_pick)
    #if __VERSION__ != 300
        varying vec4 vColor;
    #else
        flat in vec4 vColor;
    #endif
#endif
`;