export const evaluateDpoit_frag = `
precision highp float;

uniform sampler2D tDpoitFrontColor;
uniform sampler2D tDpoitBlendBackColor;
uniform vec2 uTexSize;

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec4 frontColor = texture2D(tDpoitFrontColor, coords);
    vec4 backColor = texture2D(tDpoitBlendBackColor, coords);
    float alphaMultiplier = 1.0 - frontColor.a;

    gl_FragColor = vec4(
        frontColor.rgb + alphaMultiplier * backColor.rgb,
        frontColor.a + backColor.a
    );
}
`;
