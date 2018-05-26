#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_element) || defined(dColorType_elementInstance)
    varying vec3 vColor;
#endif