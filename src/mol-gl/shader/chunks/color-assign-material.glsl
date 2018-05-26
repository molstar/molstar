#if defined(dColorType_uniform)
    vec3 material = uColor;
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_element) || defined(dColorType_elementInstance)
    vec3 material = vColor;
#endif