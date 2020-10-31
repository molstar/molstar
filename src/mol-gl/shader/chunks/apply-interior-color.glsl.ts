export default `
if (interior) {
    if (uInteriorColorFlag) {
        gl_FragColor.rgb = uInteriorColor;
    } else {
        gl_FragColor.rgb *= 1.0 - uInteriorDarkening;
    }
}
`;