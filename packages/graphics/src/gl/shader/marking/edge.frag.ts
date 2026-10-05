export const edge_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tMaskTexture;
uniform vec2 uTexSizeInv;
uniform float uEdgeScale;

// the mask is antialiased, so the coverage can be taken from a single texel
vec4 getMask(vec2 coords) {
    return texture2D(tMaskTexture, coords);
}

void main() {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;
    vec4 offset = vec4(uEdgeScale, 0.0, 0.0, uEdgeScale) * vec4(uTexSizeInv, uTexSizeInv);

    vec4 m1 = getMask(coords + offset.xy);
    vec4 m2 = getMask(coords - offset.xy);
    vec4 m3 = getMask(coords + offset.yw);
    vec4 m4 = getMask(coords - offset.yw);

    float edge = max(abs(m1.r - m2.r), abs(m3.r - m4.r));
    if (edge <= 0.0)
        discard;

    // the sample covering the most of the marked region, mask values are 0.0 where marked
    vec4 mx = m1.r < m2.r ? m1 : m2;
    vec4 my = m3.r < m4.r ? m3 : m4;
    vec4 m = mx.r < my.r ? mx : my;

    // the samples average marked and unmarked texels and unmarked ones are 1.0 in every channel,
    // so the values of the marked texels can be recovered from how much they are covered
    vec3 marked = clamp((m.gba - m.r) / max(1.0 - m.r, 0.001), 0.0, 1.0);

    float visibility = marked.x > 0.5 ? 1.0 : 0.0;
    float mask = texture2D(tMaskTexture, coords).r;
    float marker = marked.y > 0.5 ? 1.0 : 0.0;
    float fogAlpha = marked.z;
    // premultiplied, so that a blurred sample of this can be decoded the same way
    gl_FragColor = vec4(visibility, mask, marker, 1.0) * (edge * fogAlpha);
}
`;
