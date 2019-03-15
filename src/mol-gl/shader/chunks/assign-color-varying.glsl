#if defined(dColorType_attribute)
    vColor.rgb = aColor;
#elif defined(dColorType_instance)
    vColor.rgb = readFromTexture(tColor, aInstance, uColorTexDim).rgb;
#elif defined(dColorType_group)
    vColor.rgb = readFromTexture(tColor, aGroup, uColorTexDim).rgb;
#elif defined(dColorType_groupInstance)
    vColor.rgb = readFromTexture(tColor, aInstance * float(uGroupCount) + aGroup, uColorTexDim).rgb;
#elif defined(dColorType_objectPicking)
    vColor = vec4(encodeFloatRGB(float(uObjectId)), 1.0);
#elif defined(dColorType_instancePicking)
    vColor = vec4(encodeFloatRGB(aInstance), 1.0);
#elif defined(dColorType_groupPicking)
    vColor = vec4(encodeFloatRGB(aGroup), 1.0);
#endif

#ifdef dOverpaint
    vOverpaint = readFromTexture(tOverpaint, aInstance * float(uGroupCount) + aGroup, uOverpaintTexDim);
#endif

#ifdef dTransparency
    vTransparency = readFromTexture(tTransparency, aInstance * float(uGroupCount) + aGroup, uTransparencyTexDim).a;
#endif