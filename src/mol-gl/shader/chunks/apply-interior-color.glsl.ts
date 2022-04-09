export const apply_interior_color = `
if (interior) {
    if (uInteriorColorFlag) {
        gl_FragColor.rgb = uInteriorColor;
    } else {
        gl_FragColor.rgb *= 1.0 - uInteriorDarkening;
    }

    #ifdef dOpaqueBackfaces
        gl_FragColor.a = 1.0;
    #endif
}
`;