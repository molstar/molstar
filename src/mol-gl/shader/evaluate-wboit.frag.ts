export default `
precision highp float;

uniform sampler2D tWboitA;
uniform sampler2D tWboitB;

void main() {
    ivec2 coords = ivec2(gl_FragCoord.xy);
    
    vec4 accum = texelFetch(tWboitA, coords, 0);
    float r = 1.0 - accum.a;

    accum.a = texelFetch(tWboitB, coords, 0).r;
    gl_FragColor = vec4(accum.rgb / clamp(accum.a, 0.0001, 50000.0), r);
}
`;