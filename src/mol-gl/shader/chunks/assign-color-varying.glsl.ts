export const assign_color_varying = `
#if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
    #if defined(dColorType_attribute)
        vColor.rgb = aColor;
    #elif defined(dColorType_instance)
        vColor.rgb = readFromTexture(tColor, aInstance, uColorTexDim).rgb;
    #elif defined(dColorType_group)
        #if defined(dDualColor)
            vec4 color2;
            if (aColorMode == 2.0) {
                vColor.rgb = readFromTexture(tColor, group, uColorTexDim).rgb;
            } else {
                vColor.rgb = readFromTexture(tColor, group * 2.0, uColorTexDim).rgb;
                color2.rgb = readFromTexture(tColor, group * 2.0 + 1.0, uColorTexDim).rgb;
            }
        #else
            vColor.rgb = readFromTexture(tColor, group, uColorTexDim).rgb;
        #endif
    #elif defined(dColorType_groupInstance)
        #if defined(dDualColor)
            vec4 color2;
            if (aColorMode == 2.0) {
                vColor.rgb = readFromTexture(tColor, aInstance * float(uGroupCount) + group, uColorTexDim).rgb;
            } else {
                vColor.rgb = readFromTexture(tColor, (aInstance * float(uGroupCount) + group) * 2.0, uColorTexDim).rgb;
                color2.rgb = readFromTexture(tColor, (aInstance * float(uGroupCount) + group) * 2.0 + 1.0, uColorTexDim).rgb;
            }
        #else
            vColor.rgb = readFromTexture(tColor, aInstance * float(uGroupCount) + group, uColorTexDim).rgb;
        #endif
    #elif defined(dColorType_vertex)
        vColor.rgb = readFromTexture(tColor, vertexId, uColorTexDim).rgb;
    #elif defined(dColorType_vertexInstance)
        vColor.rgb = readFromTexture(tColor, int(aInstance) * uVertexCount + vertexId, uColorTexDim).rgb;
    #elif defined(dColorType_volume)
        vec3 cgridPos = (uColorGridTransform.w * (position - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, cgridPos, uColorGridDim, uColorTexDim).rgb;
    #elif defined(dColorType_volumeInstance)
        vec3 cgridPos = (uColorGridTransform.w * (vModelPosition / uModelScale - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, cgridPos, uColorGridDim, uColorTexDim).rgb;
    #endif

    #ifdef dUsePalette
        vPaletteV = ((vColor.r * 256.0 * 256.0 * 255.0 + vColor.g * 256.0 * 255.0 + vColor.b * 255.0) - 1.0) / PALETTE_SCALE;
    #endif

    #ifdef dOverpaint
        #if defined(dOverpaintType_instance)
            vOverpaint = readFromTexture(tOverpaint, aInstance, uOverpaintTexDim);
        #elif defined(dOverpaintType_groupInstance)
            vOverpaint = readFromTexture(tOverpaint, aInstance * float(uGroupCount) + group, uOverpaintTexDim);
        #elif defined(dOverpaintType_vertexInstance)
            vOverpaint = readFromTexture(tOverpaint, int(aInstance) * uVertexCount + vertexId, uOverpaintTexDim);
        #elif defined(dOverpaintType_volumeInstance)
            vec3 ogridPos = (uOverpaintGridTransform.w * (vModelPosition / uModelScale - uOverpaintGridTransform.xyz)) / uOverpaintGridDim;
            vOverpaint = texture3dFrom2dLinear(tOverpaintGrid, ogridPos, uOverpaintGridDim, uOverpaintTexDim);
        #endif

        // pre-mix to avoid darkening due to empty overpaint
        #ifdef dColorType_uniform
            vOverpaint.rgb = mix(uColor.rgb, vOverpaint.rgb, vOverpaint.a);
        #else
            vOverpaint.rgb = mix(vColor.rgb, vOverpaint.rgb, vOverpaint.a);
        #endif
        vOverpaint *= uOverpaintStrength;
    #endif

    #ifdef dEmissive
        #if defined(dEmissiveType_instance)
            vEmissive = readFromTexture(tEmissive, aInstance, uEmissiveTexDim).a;
        #elif defined(dEmissiveType_groupInstance)
            vEmissive = readFromTexture(tEmissive, aInstance * float(uGroupCount) + group, uEmissiveTexDim).a;
        #elif defined(dEmissiveType_vertexInstance)
            vEmissive = readFromTexture(tEmissive, int(aInstance) * uVertexCount + vertexId, uEmissiveTexDim).a;
        #elif defined(dEmissiveType_volumeInstance)
            vec3 egridPos = (uEmissiveGridTransform.w * (vModelPosition / uModelScale - uEmissiveGridTransform.xyz)) / uEmissiveGridDim;
            vEmissive = texture3dFrom2dLinear(tEmissiveGrid, egridPos, uEmissiveGridDim, uEmissiveTexDim).a;
        #endif
        vEmissive *= uEmissiveStrength;
    #endif

    #ifdef dSubstance
        #if defined(dSubstanceType_instance)
            vSubstance = readFromTexture(tSubstance, aInstance, uSubstanceTexDim);
        #elif defined(dSubstanceType_groupInstance)
            vSubstance = readFromTexture(tSubstance, aInstance * float(uGroupCount) + group, uSubstanceTexDim);
        #elif defined(dSubstanceType_vertexInstance)
            vSubstance = readFromTexture(tSubstance, int(aInstance) * uVertexCount + vertexId, uSubstanceTexDim);
        #elif defined(dSubstanceType_volumeInstance)
            vec3 sgridPos = (uSubstanceGridTransform.w * (vModelPosition / uModelScale - uSubstanceGridTransform.xyz)) / uSubstanceGridDim;
            vSubstance = texture3dFrom2dLinear(tSubstanceGrid, sgridPos, uSubstanceGridDim, uSubstanceTexDim);
        #endif

        // pre-mix to avoid artifacts due to empty substance
        vSubstance.rgb = mix(vec3(uMetalness, uRoughness, uBumpiness), vSubstance.rgb, vSubstance.a);
        vSubstance *= uSubstanceStrength;
    #endif
#elif defined(dRenderVariant_emissive)
    #ifdef dEmissive
        #if defined(dEmissiveType_instance)
            vEmissive = readFromTexture(tEmissive, aInstance, uEmissiveTexDim).a;
        #elif defined(dEmissiveType_groupInstance)
            vEmissive = readFromTexture(tEmissive, aInstance * float(uGroupCount) + group, uEmissiveTexDim).a;
        #elif defined(dEmissiveType_vertexInstance)
            vEmissive = readFromTexture(tEmissive, int(aInstance) * uVertexCount + vertexId, uEmissiveTexDim).a;
        #elif defined(dEmissiveType_volumeInstance)
            vec3 egridPos = (uEmissiveGridTransform.w * (vModelPosition / uModelScale - uEmissiveGridTransform.xyz)) / uEmissiveGridDim;
            vEmissive = texture3dFrom2dLinear(tEmissiveGrid, egridPos, uEmissiveGridDim, uEmissiveTexDim).a;
        #endif
        vEmissive *= uEmissiveStrength;
    #endif
#elif defined(dRenderVariant_pick)
    #ifdef requiredDrawBuffers
        vObject = vec4(packIntToRGB(float(uObjectId)), 1.0);
        vInstance = vec4(packIntToRGB(aInstance), 1.0);
        vGroup = vec4(packIntToRGB(group), 1.0);
    #else
        if (uPickType == 1) {
            vColor = vec4(packIntToRGB(float(uObjectId)), 1.0);
        } else if (uPickType == 2) {
            vColor = vec4(packIntToRGB(aInstance), 1.0);
        } else {
            vColor = vec4(packIntToRGB(group), 1.0);
        }
    #endif
#endif

#ifdef dTransparency
    #if defined(dTransparencyType_instance)
        vTransparency = readFromTexture(tTransparency, aInstance, uTransparencyTexDim).a;
    #elif defined(dTransparencyType_groupInstance)
        vTransparency = readFromTexture(tTransparency, aInstance * float(uGroupCount) + group, uTransparencyTexDim).a;
    #elif defined(dTransparencyType_vertexInstance)
        vTransparency = readFromTexture(tTransparency, int(aInstance) * uVertexCount + vertexId, uTransparencyTexDim).a;
    #elif defined(dTransparencyType_volumeInstance)
        vec3 tgridPos = (uTransparencyGridTransform.w * (vModelPosition / uModelScale - uTransparencyGridTransform.xyz)) / uTransparencyGridDim;
        vTransparency = texture3dFrom2dLinear(tTransparencyGrid, tgridPos, uTransparencyGridDim, uTransparencyTexDim).a;
    #endif
    vTransparency *= uTransparencyStrength;
#endif
`;
