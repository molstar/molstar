export default `
precision highp float;
precision highp sampler2D;

uniform sampler2D tColor;
uniform sampler2D tDepth;
uniform vec2 uTexSize;

#include common

float getDepth(const in vec2 coords) {
    return unpackRGBAToDepth(texture2D(tDepth, coords));
}

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    gl_FragColor = texture2D(tColor, coords);
    gl_FragDepthEXT = getDepth(coords);
}
`;