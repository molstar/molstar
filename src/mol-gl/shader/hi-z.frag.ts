export const hiZ_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tPreviousLevel;
uniform vec2 uInvSize;

void main(void) {
    vec2 position = gl_FragCoord.xy * uInvSize;

    float x = texture(tPreviousLevel, position).r;
    float y = textureOffset(tPreviousLevel, position, ivec2(-1, 0)).r;
    float z = textureOffset(tPreviousLevel, position, ivec2(-1, -1)).r;
    float w = textureOffset(tPreviousLevel, position, ivec2(0, -1)).r;

    gl_FragColor = vec4(max(max(x, y), max(z, w)));
}
`;