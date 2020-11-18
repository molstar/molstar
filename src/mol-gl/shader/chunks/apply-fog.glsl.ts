export default `
float fogDepth = length(vViewPosition);
float fogFactor = smoothstep(uFogNear, uFogFar, fogDepth);
float fogAlpha = (1.0 - fogFactor) * gl_FragColor.a;
float preFogAlpha = gl_FragColor.a < 1.0 ? fogAlpha : 1.0;
if (!uTransparentBackground) {
    gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
    if (gl_FragColor.a < 1.0)
        gl_FragColor.a = fogAlpha;
} else {
    gl_FragColor.a = fogAlpha;
}
`;