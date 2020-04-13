export default `
if (interior) {
    if (uInteriorColorFlag == 1) {
        gl_FragColor.rgb = uInteriorColor;
    } else {
        gl_FragColor.rgb *= 1.0 - uInteriorDarkening;
    }
}
`;