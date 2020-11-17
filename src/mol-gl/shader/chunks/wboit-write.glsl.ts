export default `
#if defined(dRenderVariant_colorWboit)
    if (uRenderWboit == 0) {
        if (preFogAlpha < 1.0) {
            discard;
        }
    } else if (uRenderWboit == 1) {
        if (preFogAlpha != 1.0 && !interior && fragmentDepth < getDepth(gl_FragCoord.xy / uDrawingBufferSize)) {
            float alpha = preFogAlpha;
            float wboitWeight = alpha * clamp(pow(1.0 - fragmentDepth, 2.0), 0.01, 1.0);
            gl_FragColor = vec4(gl_FragColor.rgb * alpha * wboitWeight, alpha);
            gl_FragData[1] = vec4(alpha * wboitWeight);
        } else {
            discard;
        }
    }
#endif
`;