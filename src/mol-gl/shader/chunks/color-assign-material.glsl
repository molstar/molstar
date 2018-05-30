#if defined(dColorType_uniform)
    vec4 material.rgb = uColor;
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_element) || defined(dColorType_elementInstance) || defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_elementPicking)
    vec4 material = vColor;
#endif