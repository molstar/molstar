export default `
#if defined(dRenderVariant_color)
    #if defined(dColorType_uniform)
        vec4 material = vec4(uColor, uAlpha);
    #elif defined(dColorType_varying)
        vec4 material = vec4(vColor.rgb, uAlpha);
    #endif

    // mix material with overpaint
    #if defined(dOverpaint)
        material.rgb = mix(material.rgb, vOverpaint.rgb, vOverpaint.a);
    #endif

    // apply screendoor transparency
    #if defined(dTransparency)
        float ta = 1.0 - vTransparency;
        float at = 0.0;

        // shift by view-offset during multi-sample rendering to allow for blending
        vec2 coord = gl_FragCoord.xy + uViewOffset * 0.25;

        #if defined(dTransparencyVariant_single)
            const mat4 thresholdMatrix = mat4(
                1.0 / 17.0,  9.0 / 17.0,  3.0 / 17.0, 11.0 / 17.0,
                13.0 / 17.0,  5.0 / 17.0, 15.0 / 17.0,  7.0 / 17.0,
                4.0 / 17.0, 12.0 / 17.0,  2.0 / 17.0, 10.0 / 17.0,
                16.0 / 17.0,  8.0 / 17.0, 14.0 / 17.0,  6.0 / 17.0
            );
            at = thresholdMatrix[int(intMod(coord.x, 4.0))][int(intMod(coord.y, 4.0))];
        #elif defined(dTransparencyVariant_multi)
            at = fract(dot(vec3(coord, vGroup + 0.5), vec3(2.0, 7.0, 23.0) / 17.0f));
        #endif

        if (ta < 0.99 && (ta < 0.01 || ta < at)) discard;
    #endif
#elif defined(dRenderVariant_pick)
    vec4 material = vColor;
#elif defined(dRenderVariant_depth)
    #ifdef enabledFragDepth
        vec4 material = packDepthToRGBA(gl_FragDepthEXT);
    #else
        vec4 material = packDepthToRGBA(gl_FragCoord.z);
    #endif
#endif
`;