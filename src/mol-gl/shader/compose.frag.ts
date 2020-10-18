export default `
precision highp float;
precision highp sampler2D;

uniform vec2 uQuadShift;

uniform sampler2D tColor;
uniform vec2 uTexSize;
uniform float uWeight;

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize - uQuadShift;
    gl_FragColor = texture2D(tColor, coords) * uWeight;
}
`;