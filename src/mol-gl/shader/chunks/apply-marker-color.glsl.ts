export default `
float marker = floor(vMarker * 255.0 + 0.5); // rounding required to work on some cards on win
if (marker > 0.1) {
    if (intMod(marker, 2.0) > 0.1) {
        gl_FragColor.rgb = mix(uHighlightColor, gl_FragColor.rgb, 0.3);
    } else {
        gl_FragColor.rgb = mix(uSelectColor, gl_FragColor.rgb, 0.3);
    }
}
`;