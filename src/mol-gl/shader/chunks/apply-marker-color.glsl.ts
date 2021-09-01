export const apply_marker_color = `
if (marker > 0.1) {
    if (intMod(marker, 2.0) > 0.1) {
        gl_FragColor.rgb = mix(gl_FragColor.rgb, uHighlightColor, uHighlightStrength);
        gl_FragColor.a = max(gl_FragColor.a, uHighlightStrength * 0.002); // for direct-volume rendering
    } else {
        gl_FragColor.rgb = mix(gl_FragColor.rgb, uSelectColor, uSelectStrength);
        gl_FragColor.a = max(gl_FragColor.a, uSelectStrength * 0.002); // for direct-volume rendering
    }
}
`;