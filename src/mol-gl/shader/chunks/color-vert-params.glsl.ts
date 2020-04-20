export default `
#if defined(dRenderVariant_color)
    #if defined(dColorType_uniform)
        uniform vec3 uColor;
    #elif defined(dColorType_attribute)
        varying vec4 vColor;
        attribute vec3 aColor;
    #elif defined(dColorType_texture)
        varying vec4 vColor;
        uniform vec2 uColorTexDim;
        uniform sampler2D tColor;
    #endif

    #ifdef dOverpaint
        varying vec4 vOverpaint;
        uniform vec2 uOverpaintTexDim;
        uniform sampler2D tOverpaint;
    #endif

    #ifdef dTransparency
        varying float vGroup;
        varying float vTransparency;
        uniform vec2 uTransparencyTexDim;
        uniform sampler2D tTransparency;
    #endif
#elif defined(dRenderVariant_pick)
    #if __VERSION__ != 300
        varying vec4 vColor;
    #else
        flat out vec4 vColor;
    #endif
#endif
`;