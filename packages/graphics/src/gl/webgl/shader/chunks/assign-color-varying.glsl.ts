/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

export const assign_color_varying = `
#if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
    #if defined(dColorType_texture)
        // the texture's index by the color type (\`ColorTypeIndex\`): instance, group, groupInstance, vertex,
        // vertexInstance
        float colorIndex = aInstance;
        if (uColorIndexType == 1) {
            colorIndex = groupIndex(group);
        } else if (uColorIndexType == 2) {
            colorIndex = groupInstanceIndex(aInstance, group);
        } else if (uColorIndexType == 3) {
            colorIndex = float(vertexId);
        } else if (uColorIndexType == 4) {
            colorIndex = float(int(aInstance) * uVertexCount + vertexId);
        }
        #if defined(dDualColor)
            // two colors per group, except for color mode 2 (cylinders)
            bool dualColor = uColorIndexType == 1 || uColorIndexType == 2;
            vec4 color2;
            if (dualColor && aColorMode != 2.0) {
                vColor.rgb = read_tColor(colorIndex * 2.0).rgb;
                color2.rgb = read_tColor(colorIndex * 2.0 + 1.0).rgb;
            } else {
                vColor.rgb = read_tColor(colorIndex).rgb;
            }
        #else
            vColor.rgb = read_tColor(colorIndex).rgb;
        #endif
    #elif defined(dColorType_grid)
        // the local position (volume) or the model position (volumeInstance), see \`ColorTypeIndex\`
        vec3 cgridSource = uColorIndexType == 1 ? vModelPosition / uModelScale : position;
        vec3 cgridPos = (uColorGridTransform.w * (cgridSource - uColorGridTransform.xyz)) / uColorGridDim;
        vColor.rgb = texture3dFrom2dLinear(tColorGrid, cgridPos, uColorGridDim, uColorGridTexDim).rgb;
    #endif

    #ifdef dUsePalette
        vPaletteV = ((vColor.r * 256.0 * 256.0 * 255.0 + vColor.g * 256.0 * 255.0 + vColor.b * 255.0) - 1.0) / PALETTE_SCALE;
    #endif

    #ifdef dOverpaint
        #if defined(dOverpaintType_texture)
            vOverpaint = read_tOverpaint(granularityIndex(aInstance, group));
        #elif defined(dOverpaintType_vertexInstance)
            vOverpaint = read_tOverpaint(int(aInstance) * uVertexCount + vertexId);
        #elif defined(dOverpaintType_volumeInstance)
            vec3 ogridPos = (uOverpaintGridTransform.w * (vModelPosition / uModelScale - uOverpaintGridTransform.xyz)) / uOverpaintGridDim;
            vOverpaint = texture3dFrom2dLinear(tOverpaintGrid, ogridPos, uOverpaintGridDim, uOverpaintGridTexDim);
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
        #if defined(dEmissiveType_texture)
            vEmissive = read_tEmissive(granularityIndex(aInstance, group)).r;
        #elif defined(dEmissiveType_vertexInstance)
            vEmissive = read_tEmissive(int(aInstance) * uVertexCount + vertexId).r;
        #elif defined(dEmissiveType_volumeInstance)
            vec3 egridPos = (uEmissiveGridTransform.w * (vModelPosition / uModelScale - uEmissiveGridTransform.xyz)) / uEmissiveGridDim;
            vEmissive = texture3dFrom2dLinear(tEmissiveGrid, egridPos, uEmissiveGridDim, uEmissiveGridTexDim).r;
        #endif
        vEmissive *= uEmissiveStrength;
    #endif

    #ifdef dSubstance
        #if defined(dSubstanceType_texture)
            vSubstance = read_tSubstance(granularityIndex(aInstance, group));
        #elif defined(dSubstanceType_vertexInstance)
            vSubstance = read_tSubstance(int(aInstance) * uVertexCount + vertexId);
        #elif defined(dSubstanceType_volumeInstance)
            vec3 sgridPos = (uSubstanceGridTransform.w * (vModelPosition / uModelScale - uSubstanceGridTransform.xyz)) / uSubstanceGridDim;
            vSubstance = texture3dFrom2dLinear(tSubstanceGrid, sgridPos, uSubstanceGridDim, uSubstanceGridTexDim);
        #endif

        // pre-mix to avoid artifacts due to empty substance
        vSubstance.rgb = mix(vec3(uMetalness, uRoughness, uBumpiness), vSubstance.rgb, vSubstance.a);
        vSubstance *= uSubstanceStrength;
    #endif
#elif defined(dRenderVariant_emissive)
    #ifdef dEmissive
        #if defined(dEmissiveType_texture)
            vEmissive = read_tEmissive(granularityIndex(aInstance, group)).r;
        #elif defined(dEmissiveType_vertexInstance)
            vEmissive = read_tEmissive(int(aInstance) * uVertexCount + vertexId).r;
        #elif defined(dEmissiveType_volumeInstance)
            vec3 egridPos = (uEmissiveGridTransform.w * (vModelPosition / uModelScale - uEmissiveGridTransform.xyz)) / uEmissiveGridDim;
            vEmissive = texture3dFrom2dLinear(tEmissiveGrid, egridPos, uEmissiveGridDim, uEmissiveGridTexDim).r;
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
    #if defined(dTransparencyType_texture)
        vTransparency = read_tTransparency(granularityIndex(aInstance, group)).r;
    #elif defined(dTransparencyType_vertexInstance)
        vTransparency = read_tTransparency(int(aInstance) * uVertexCount + vertexId).r;
    #elif defined(dTransparencyType_volumeInstance)
        vec3 tgridPos = (uTransparencyGridTransform.w * (vModelPosition / uModelScale - uTransparencyGridTransform.xyz)) / uTransparencyGridDim;
        vTransparency = texture3dFrom2dLinear(tTransparencyGrid, tgridPos, uTransparencyGridDim, uTransparencyGridTexDim).r;
    #endif
    vTransparency *= uTransparencyStrength;
#endif
`;
