export const apply_marker_color = `
float marker = floor(vMarker * 255.0 + 0.5); // rounding required to work on some cards on win
if (marker > 0.1) {
    if (intMod(marker, 2.0) > 0.1) {
        gl_FragColor.rgb = mix(uHighlightColor, gl_FragColor.rgb, 0.3);
        gl_FragColor.a = max(0.02, gl_FragColor.a); // for direct-volume rendering
    } else {
        gl_FragColor.rgb = mix(uSelectColor, gl_FragColor.rgb, 0.3);
    }
}
`;