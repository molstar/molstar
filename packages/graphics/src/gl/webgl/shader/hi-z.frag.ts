export const hiZ_frag = `
precision highp float;

void main(void) {
    vec2 position = gl_FragCoord.xy * uInvSize + uRect.xy;

    float x = texture(tPreviousLevel, position).r;
    float y = textureOffset(tPreviousLevel, position, ivec2(-1, 0)).r;
    float z = textureOffset(tPreviousLevel, position, ivec2(-1, -1)).r;
    float w = textureOffset(tPreviousLevel, position, ivec2(0, -1)).r;

    // the farthest depth (depth is reversed)
    gl_FragColor = vec4(min(min(x, y), min(z, w)));
}
`;
