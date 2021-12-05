export const color_vert_params = `
uniform float uMetalness;
uniform float uRoughness;
uniform float uBumpiness;

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
    #elif defined(dColorType_grid)
        varying vec4 vColor;
        uniform vec2 uColorTexDim;
        uniform vec3 uColorGridDim;
        uniform vec4 uColorGridTransform;
        uniform sampler2D tColorGrid;
    #endif

    #ifdef dOverpaint
        #if defined(dOverpaintType_groupInstance) || defined(dOverpaintType_vertexInstance)
            varying vec4 vOverpaint;
            uniform vec2 uOverpaintTexDim;
            uniform sampler2D tOverpaint;
        #elif defined(dOverpaintType_volumeInstance)
            varying vec4 vOverpaint;
            uniform vec2 uOverpaintTexDim;
            uniform vec3 uOverpaintGridDim;
            uniform vec4 uOverpaintGridTransform;
            uniform sampler2D tOverpaintGrid;
        #endif
    #endif

    #ifdef dSubstance
        #if defined(dSubstanceType_groupInstance) || defined(dSubstanceType_vertexInstance)
            varying vec4 vSubstance;
            uniform vec2 uSubstanceTexDim;
            uniform sampler2D tSubstance;
        #elif defined(dSubstanceType_volumeInstance)
            varying vec4 vSubstance;
            uniform vec2 uSubstanceTexDim;
            uniform vec3 uSubstanceGridDim;
            uniform vec4 uSubstanceGridTransform;
            uniform sampler2D tSubstanceGrid;
        #endif
    #endif
#elif defined(dRenderVariant_pick)
    #if __VERSION__ == 100
        varying vec4 vColor;
    #else
        flat out vec4 vColor;
    #endif
#endif

#ifdef dTransparency
    varying float vGroup;
    #if defined(dTransparencyType_groupInstance) || defined(dTransparencyType_vertexInstance)
        varying float vTransparency;
        uniform vec2 uTransparencyTexDim;
        uniform sampler2D tTransparency;
    #elif defined(dTransparencyType_volumeInstance)
        varying float vTransparency;
        uniform vec2 uTransparencyTexDim;
        uniform vec3 uTransparencyGridDim;
        uniform vec4 uTransparencyGridTransform;
        uniform sampler2D tTransparencyGrid;
    #endif
#endif

#ifdef dUsePalette
    varying float vPaletteV;
#endif
`;