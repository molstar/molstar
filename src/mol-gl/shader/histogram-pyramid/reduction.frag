precision highp float;
precision highp sampler2D;

// input texture (previous level used to evaluate the new level)
uniform sampler2D tPreviousLevel;

// 1/size of the previous level texture.
uniform float uSize;

varying vec2 vCoordinate;

void main(void) {
    float k = 0.5 * uSize;
    vec2 position = floor(vCoordinate / uSize) * uSize;
    float a = texture2D(tPreviousLevel, position + vec2(0., 0.)).r;
    float b = texture2D(tPreviousLevel, position + vec2(k, 0.)).r;
    float c = texture2D(tPreviousLevel, position + vec2(0., k)).r;
    float d = texture2D(tPreviousLevel, position + vec2(k, k)).r;
    gl_FragColor.a = a;
    gl_FragColor.b = a + b;
    gl_FragColor.g = gl_FragColor.b + c;
    gl_FragColor.r = gl_FragColor.g + d;
}