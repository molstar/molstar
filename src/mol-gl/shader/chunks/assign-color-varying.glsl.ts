export default `
#if defined(dRenderVariant_color)
    #if defined(dColorType_attribute)
        vColor.rgb = aColor;
    #elif defined(dColorType_instance)
        vColor.rgb = readFromTexture(tColor, aInstance, uColorTexDim).rgb;
    #elif defined(dColorType_group)
        vColor.rgb = readFromTexture(tColor, group, uColorTexDim).rgb;
    #elif defined(dColorType_groupInstance)
        vColor.rgb = readFromTexture(tColor, aInstance * float(uGroupCount) + group, uColorTexDim).rgb;
    #endif

    #ifdef dOverpaint
        vOverpaint = readFromTexture(tOverpaint, aInstance * float(uGroupCount) + group, uOverpaintTexDim);
    #endif

    #ifdef dTransparency
        vGroup = group;
        vTransparency = readFromTexture(tTransparency, aInstance * float(uGroupCount) + group, uTransparencyTexDim).a;
    #endif
#elif defined(dRenderVariant_pick)
    #if defined(dRenderVariant_pickObject)
        vColor = vec4(encodeFloatRGB(float(uObjectId)), 1.0);
    #elif defined(dRenderVariant_pickInstance)
        vColor = vec4(encodeFloatRGB(aInstance), 1.0);
    #elif defined(dRenderVariant_pickGroup)
        vColor = vec4(encodeFloatRGB(group), 1.0);
    #endif
#endif
`;