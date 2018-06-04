#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_element) || defined(dColorType_elementInstance) || defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_elementPicking)
    varying vec4 vColor;
#endif