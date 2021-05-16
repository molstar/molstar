export const assign_color_varying = `
#if defined(dRenderVariant_color)
    #if defined(dColorType_attribute)
        vColor.rgb = aColor;
    #elif defined(dColorType_instance)
        vColor.rgb = readFromTexture(tColor, aInstance, uColorTexDim).rgb;
    #elif defined(dColorType_group)
        vColor.rgb = readFromTexture(tColor, group, uColorTexDim).rgb;
    #elif defined(dColorType_groupInstance)
        vColor.rgb = readFromTexture(tColor, aInstance * float(uGroupCount) + group, uColorTexDim).rgb;
    #elif defined(dColorType_vertex)
        vColor.rgb = readFromTexture(tColor, VertexID, uColorTexDim).rgb;
    #elif defined(dColorType_vertexInstance)
        vColor.rgb = readFromTexture(tColor, int(aInstance) * uVertexCount + VertexID, uColorTexDim).rgb;
    #elif defined(dColorType_volume)
        vec3 gridPos = (uColorGridTransform.w * (position - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, gridPos, uColorGridDim, uColorTexDim).rgb;
    #elif defined(dColorType_volumeInstance)
        vec3 gridPos = (uColorGridTransform.w * (vModelPosition - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, gridPos, uColorGridDim, uColorTexDim).rgb;
    #endif

    #ifdef dUsePalette
        vPaletteV = ((vColor.r * 256.0 * 256.0 * 255.0 + vColor.g * 256.0 * 255.0 + vColor.b * 255.0) - 1.0) / 16777215.0;
    #endif

    #ifdef dOverpaint
        vOverpaint = readFromTexture(tOverpaint, aInstance * float(uGroupCount) + group, uOverpaintTexDim);
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

#ifdef dTransparency
    vGroup = group;
    vTransparency = readFromTexture(tTransparency, aInstance * float(uGroupCount) + group, uTransparencyTexDim).a;
#endif
`;