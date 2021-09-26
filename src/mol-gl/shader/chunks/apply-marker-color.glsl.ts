export const apply_marker_color = `
if (marker > 0.0) {
    if ((uMarkerPriority == 1 && marker != 2.0) || (uMarkerPriority != 1 && marker == 1.0)) {
        gl_FragColor.rgb = mix(gl_FragColor.rgb, uHighlightColor, uHighlightStrength);
        gl_FragColor.a = max(gl_FragColor.a, uHighlightStrength * 0.002); // for direct-volume rendering
    } else {
        gl_FragColor.rgb = mix(gl_FragColor.rgb, uSelectColor, uSelectStrength);
        gl_FragColor.a = max(gl_FragColor.a, uSelectStrength * 0.002); // for direct-volume rendering
    }
}
`;