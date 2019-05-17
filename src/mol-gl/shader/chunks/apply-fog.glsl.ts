export default `
#ifdef dUseFog
	float depth = length(vViewPosition);
    float fogFactor = smoothstep(uFogNear, uFogFar, depth);
	gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
    float fogAlpha = (1.0 - fogFactor) * gl_FragColor.a;
    if (fogAlpha < 0.01)
        discard;
    gl_FragColor = vec4(gl_FragColor.rgb, fogAlpha);
#endif
`