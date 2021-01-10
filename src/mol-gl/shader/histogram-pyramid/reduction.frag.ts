export default `
precision highp float;
precision highp sampler2D;

// input texture (previous level used to evaluate the new level)
uniform sampler2D tPreviousLevel;

// inverted size of the previous level texture.
uniform float uSize;
uniform float uTexSize;
uniform bool uFirst;

void main(void) {
    float k = 0.5 * uSize;
    vec2 position = floor((gl_FragCoord.xy / uTexSize) / uSize) * uSize;
    float a, b, c, d;

    if (uFirst) {
        a = texture2D(tPreviousLevel, position).r * 255.0;
        b = texture2D(tPreviousLevel, position + vec2(k, 0.0)).r * 255.0;
        c = texture2D(tPreviousLevel, position + vec2(0.0, k)).r * 255.0;
        d = texture2D(tPreviousLevel, position + vec2(k, k)).r * 255.0;
    } else {
        a = texture2D(tPreviousLevel, position).r;
        b = texture2D(tPreviousLevel, position + vec2(k, 0.0)).r;
        c = texture2D(tPreviousLevel, position + vec2(0.0, k)).r;
        d = texture2D(tPreviousLevel, position + vec2(k, k)).r;
    }

    gl_FragColor.a = a;
    gl_FragColor.b = a + b;
    gl_FragColor.g = gl_FragColor.b + c;
    gl_FragColor.r = gl_FragColor.g + d;
}
`;