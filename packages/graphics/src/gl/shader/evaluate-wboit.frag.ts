export const evaluateWboit_frag = `
precision highp float;

uniform sampler2D tWboitA;
uniform sampler2D tWboitB;
uniform vec2 uTexSize;

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;

    vec4 accum = texture2D(tWboitA, coords);
    float r = 1.0 - accum.a;

    accum.a = texture2D(tWboitB, coords).r;
    // divisor needs to allow very small values for nice fading
    gl_FragColor = vec4(accum.rgb / clamp(accum.a, 0.00000001, 50000.0), r);
}
`;
