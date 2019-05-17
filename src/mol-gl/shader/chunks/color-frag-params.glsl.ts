export default `
#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance)
    varying vec4 vColor;
#elif defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    #if __VERSION__ != 300
        varying vec4 vColor;
    #else
        flat in vec4 vColor;
    #endif
#endif

#ifdef dOverpaint
    varying vec4 vOverpaint;
#endif

#ifdef dTransparency
    varying float vGroup;
    varying float vTransparency;
#endif
`