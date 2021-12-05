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
        vec3 cgridPos = (uColorGridTransform.w * (position - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, cgridPos, uColorGridDim, uColorTexDim).rgb;
    #elif defined(dColorType_volumeInstance)
        vec3 cgridPos = (uColorGridTransform.w * (vModelPosition - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, cgridPos, uColorGridDim, uColorTexDim).rgb;
    #endif

    #ifdef dUsePalette
        vPaletteV = ((vColor.r * 256.0 * 256.0 * 255.0 + vColor.g * 256.0 * 255.0 + vColor.b * 255.0) - 1.0) / 16777215.0;
    #endif

    #ifdef dOverpaint
        #if defined(dOverpaintType_groupInstance)
            vOverpaint = readFromTexture(tOverpaint, aInstance * float(uGroupCount) + group, uOverpaintTexDim);
        #elif defined(dOverpaintType_vertexInstance)
            vOverpaint = readFromTexture(tOverpaint, int(aInstance) * uVertexCount + VertexID, uOverpaintTexDim);
        #elif defined(dOverpaintType_volumeInstance)
            vec3 ogridPos = (uOverpaintGridTransform.w * (vModelPosition - uOverpaintGridTransform.xyz)) / uOverpaintGridDim;
            vOverpaint = texture3dFrom2dLinear(tOverpaintGrid, ogridPos, uOverpaintGridDim, uOverpaintTexDim);
        #endif

        // pre-mix to avoid darkening due to empty overpaint
        #ifdef dColorType_uniform
            vOverpaint.rgb = mix(uColor.rgb, vOverpaint.rgb, vOverpaint.a);
        #else
            vOverpaint.rgb = mix(vColor.rgb, vOverpaint.rgb, vOverpaint.a);
        #endif
    #endif

    #ifdef dSubstance
        #if defined(dSubstanceType_groupInstance)
            vSubstance = readFromTexture(tSubstance, aInstance * float(uGroupCount) + group, uSubstanceTexDim);
        #elif defined(dSubstanceType_vertexInstance)
            vSubstance = readFromTexture(tSubstance, int(aInstance) * uVertexCount + VertexID, uSubstanceTexDim);
        #elif defined(dSubstanceType_volumeInstance)
            vec3 sgridPos = (uSubstanceGridTransform.w * (vModelPosition - uSubstanceGridTransform.xyz)) / uSubstanceGridDim;
            vSubstance = texture3dFrom2dLinear(tSubstanceGrid, sgridPos, uSubstanceGridDim, uSubstanceTexDim);
        #endif

        // pre-mix to avoid artifacts due to empty substance
        vSubstance.rgb = mix(vec3(uMetalness, uRoughness, uBumpiness), vSubstance.rgb, vSubstance.a);
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

    #if defined(dTransparencyType_groupInstance)
        vTransparency = readFromTexture(tTransparency, aInstance * float(uGroupCount) + group, uTransparencyTexDim).a;
    #elif defined(dTransparencyType_vertexInstance)
        vTransparency = readFromTexture(tTransparency, int(aInstance) * uVertexCount + VertexID, uTransparencyTexDim).a;
    #elif defined(dTransparencyType_volumeInstance)
        vec3 tgridPos = (uTransparencyGridTransform.w * (vModelPosition - uTransparencyGridTransform.xyz)) / uTransparencyGridDim;
        vTransparency = texture3dFrom2dLinear(tTransparencyGrid, tgridPos, uTransparencyGridDim, uTransparencyTexDim).a;
    #endif
#endif
`;