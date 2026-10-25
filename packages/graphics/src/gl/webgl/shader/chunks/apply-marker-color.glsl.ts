export const apply_marker_color = `

#if defined(dColorMarker)
    if (marker > 0.0) {
        if ((uMarkerPriority == 1 && marker != 2.0) || (uMarkerPriority != 1 && marker == 1.0)) {
            gl_FragColor.rgb = mix(gl_FragColor.rgb, uHighlightColor, uHighlightStrength);
        } else {
            gl_FragColor.rgb = mix(gl_FragColor.rgb, uSelectColor, uSelectStrength);
        }
    } else if (uMarkerAverage > 0.0) {
        gl_FragColor.rgb = mix(gl_FragColor.rgb, uDimColor, uDimStrength);
    }
#endif
`;
