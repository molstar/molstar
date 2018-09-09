#if defined(dColorType_uniform)
    vec4 material = vec4(uColor, uAlpha);
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance) || defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    vec4 material = vColor;
#endif