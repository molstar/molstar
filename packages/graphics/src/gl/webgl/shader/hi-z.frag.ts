/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const hiZ_frag = `
precision highp float;

// the farthest depth (the min, depth is reversed) of the texels of the previous level (or the depth texture) this texel
// covers: 2x2 of a level, up to 3x3 of the depth texture (level 0 maps the viewport onto a power-of-two square, e.g.
// 1280 pixels onto 1024 texels: a texel's span can touch three pixels)
void main(void) {
    ivec2 size = textureSize(tPreviousLevel, 0);
    vec2 s = vec2(size);
    vec2 start = ((gl_FragCoord.xy - 0.5) * uInvSize + uRect.xy) * s;
    vec2 end = ((gl_FragCoord.xy + 0.5) * uInvSize + uRect.xy) * s;
    ivec2 lo = clamp(ivec2(floor(start)), ivec2(0), size - 1);
    ivec2 hi = clamp(ivec2(ceil(end)) - 1, lo, min(lo + 2, size - 1));

    float depth = 1.0;
    for (int dy = 0; dy < 3; ++dy) {
        for (int dx = 0; dx < 3; ++dx) {
            ivec2 p = lo + ivec2(dx, dy);
            if (p.x <= hi.x && p.y <= hi.y) depth = min(depth, texelFetch(tPreviousLevel, p, 0).r);
        }
    }
    gl_FragColor = vec4(depth);
}
`;
