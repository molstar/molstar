export const edge_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tMaskTexture;
uniform vec2 uTexSizeInv;

void main() {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;
    vec4 offset = vec4(float(dEdgeScale), 0.0, 0.0, float(dEdgeScale)) * vec4(uTexSizeInv, uTexSizeInv);
    vec4 c0 = texture2D(tMaskTexture, coords);
    vec4 c1 = texture2D(tMaskTexture, coords + offset.xy);
    vec4 c2 = texture2D(tMaskTexture, coords - offset.xy);
    vec4 c3 = texture2D(tMaskTexture, coords + offset.yw);
    vec4 c4 = texture2D(tMaskTexture, coords - offset.yw);
    float diff1 = (c1.r - c2.r) * 0.5;
    float diff2 = (c3.r - c4.r) * 0.5;
    float d = length(vec2(diff1, diff2));
    if (d <= 0.0)
        discard;
    float a1 = min(c1.g, c2.g);
    float a2 = min(c3.g, c4.g);
    float visibility = min(a1, a2) > 0.001 ? 1.0 : 0.0;
    float mask = c0.r;
    float marker = min(c1.b, min(c2.b, min(c3.b, c4.b)));
    float fogAlpha = min(c1.a, min(c2.a, min(c3.a, c4.a)));
    gl_FragColor = vec4(visibility, mask, marker, fogAlpha);
}
`;