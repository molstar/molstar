export default `
#if defined(dRenderVariant_colorWboit)
    if (uRenderWboit == 0) {
        if (gl_FragColor.a != 1.0) {
            discard;
        }
    } else if (uRenderWboit == 1) {
        if (gl_FragColor.a != 1.0 && absFragDepth < getDepth(gl_FragCoord.xy / uDrawingBufferSize)) {
            float alpha = gl_FragColor.a;
            float wboitWeight = alpha * clamp(pow(1.0 - absFragDepth, 2.0), 0.01, 1.0);
            gl_FragColor = vec4(gl_FragColor.rgb * alpha * wboitWeight, alpha);
            gl_FragData[1] = vec4(alpha * wboitWeight);
        } else {
            discard;
        }
    }
#endif
`;