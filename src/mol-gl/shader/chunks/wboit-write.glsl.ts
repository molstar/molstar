export default `
#if defined(dRenderVariant_colorWboit)
    if (!uRenderWboit) {
        if (preFogAlpha < 1.0) {
            discard;
        }
    } else if (uRenderWboit) {
        if (preFogAlpha != 1.0 && !interior && fragmentDepth < getDepth(gl_FragCoord.xy / uDrawingBufferSize)) {
            float alpha = gl_FragColor.a;
            float wboitWeight = alpha * clamp(pow(1.0 - fragmentDepth, 2.0), 0.01, 1.0);
            gl_FragColor = vec4(gl_FragColor.rgb * alpha * wboitWeight, alpha);
            // extra alpha is to handle pre-multiplied alpha
            #if !defined(dRenderMode_volume) && !defined(dRenderMode_isosurface)
                gl_FragData[1] = vec4((uTransparentBackground ? alpha : 1.0) * alpha * wboitWeight);
            #else
                gl_FragData[1] = vec4(alpha * alpha * wboitWeight);
            #endif
        } else {
            discard;
        }
    }
#endif
`;