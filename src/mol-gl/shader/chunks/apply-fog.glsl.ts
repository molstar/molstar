export default `
float depth = length(vViewPosition);
float fogFactor = smoothstep(uFogNear, uFogFar, depth);
float fogAlpha = (1.0 - fogFactor) * gl_FragColor.a;
if (uTransparentBackground == 0) {
    gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
    if (gl_FragColor.a < 1.0)
        gl_FragColor.a = fogAlpha;
} else {
    float fogAlpha = (1.0 - fogFactor) * gl_FragColor.a;
    gl_FragColor.a = fogAlpha;
}
`;