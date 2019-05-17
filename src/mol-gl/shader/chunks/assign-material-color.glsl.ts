export default `
#if defined(dColorType_uniform)
    vec4 material = vec4(uColor, uAlpha);
#elif defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance)
    vec4 material = vec4(vColor.rgb, uAlpha);
#elif defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    vec4 material = uPickable == 1 ? vColor : vec4(0.0, 0.0, 0.0, 1.0); // set to empty picking id
#elif defined(dColorType_depth)
    #ifdef enabledFragDepth
        vec4 material = packDepthToRGBA(gl_FragDepthEXT);
    #else
        vec4 material = packDepthToRGBA(gl_FragCoord.z);
    #endif
#endif

// mix material with overpaint
#if defined(dOverpaint) && (defined(dColorType_uniform) || defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance))
    material.rgb = mix(material.rgb, vOverpaint.rgb, vOverpaint.a);
#endif

// apply screendoor transparency
#if defined(dTransparency) && (defined(dColorType_uniform) || defined(dColorType_attribute) || defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance))
    float ta = 1.0 - vTransparency;
    float at = 0.0;

    #if defined(dTransparencyVariant_single)
        const mat4 thresholdMatrix = mat4(
            1.0 / 17.0,  9.0 / 17.0,  3.0 / 17.0, 11.0 / 17.0,
            13.0 / 17.0,  5.0 / 17.0, 15.0 / 17.0,  7.0 / 17.0,
            4.0 / 17.0, 12.0 / 17.0,  2.0 / 17.0, 10.0 / 17.0,
            16.0 / 17.0,  8.0 / 17.0, 14.0 / 17.0,  6.0 / 17.0
        );
        at = thresholdMatrix[int(intMod(gl_FragCoord.x, 4.0))][int(intMod(gl_FragCoord.y, 4.0))];
    #elif defined(dTransparencyVariant_multi)
        at = fract(dot(vec3(gl_FragCoord.xy, vGroup + 0.5), vec3(2.0, 7.0, 23.0) / 17.0f));
    #endif

    if (ta < 0.99 && (ta < 0.01 || ta < at)) discard;
#endif
`