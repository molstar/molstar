#ifdef dUseFog
	float depth = length(vViewPosition);
    // float depth = gl_FragCoord.z / gl_FragCoord.w;
    float fogFactor = smoothstep(uFogNear, uFogFar, depth);
	gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
    float alpha = (1.0 - fogFactor) * gl_FragColor.a;
    if (alpha < 0.01)
        discard;
    gl_FragColor = vec4( gl_FragColor.rgb, alpha );
#endif