export const apply_fog = `
float viewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
float fogFactor = smoothstep(uFogNear, uFogFar, abs(viewZ));
float fogAlpha = (1.0 - fogFactor) * gl_FragColor.a;
float preFogAlpha = gl_FragColor.a;
if (!uTransparentBackground) {
    if (gl_FragColor.a < 1.0) {
        // transparent objects are blended with background color
        gl_FragColor.a = fogAlpha;
    } else {
        // mix opaque objects with background color
        gl_FragColor.rgb = mix(gl_FragColor.rgb, uFogColor, fogFactor);
    }
} else {
    // pre-multiplied alpha expected for transparent background
    gl_FragColor.rgb *= fogAlpha;
    gl_FragColor.a = fogAlpha;
}
`;