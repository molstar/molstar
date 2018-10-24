float marker = vMarker * 255.0;
if (marker > 0.1) {
    if (mod(marker, 2.0) > 0.1) {
        gl_FragColor.rgb = mix(uHighlightColor, gl_FragColor.rgb, 0.3);
    } else {
        gl_FragColor.rgb = mix(uSelectColor, gl_FragColor.rgb, 0.3);
    }
}