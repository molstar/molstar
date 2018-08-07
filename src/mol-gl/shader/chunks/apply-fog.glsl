#ifdef dUseFog
	float depth = length(vViewPosition);
    // float depth = gl_FragCoord.z / gl_FragCoord.w;
    float fogFactor = smoothstep(uFogNear, uFogFar, depth);
	gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
#endif