export default `
float depth = length(vViewPosition);
float fogFactor = smoothstep(uFogNear, uFogFar, depth);
if (uTransparentBackground == 0) {
    gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
} else {
    float fogAlpha = (1.0 - fogFactor) * gl_FragColor.a;
    gl_FragColor.a = fogAlpha;
}
`