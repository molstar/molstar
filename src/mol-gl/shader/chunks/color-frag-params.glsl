#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance) || defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    varying vec4 vColor;
#endif