#if defined(dColorType_uniform)
    vec4 material = vec4(uColor, uAlpha);
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance)
    vec4 material = vec4(vColor.rgb, uAlpha);
#elif defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    vec4 material = uPickable == 1 ? vColor : vec4(0.0, 0.0, 0.0, 1.0); // set to empty picking id
#endif

// mix material with overpaint
#if defined(dOverpaint) && (defined(dColorType_uniform) || defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance))
    material.rgb = mix(material.rgb, vOverpaint.rgb, vOverpaint.a);
#endif