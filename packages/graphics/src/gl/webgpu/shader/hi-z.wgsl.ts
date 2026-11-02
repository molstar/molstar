/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl } from './chunks/quad.wgsl.js';

/**
 * A level of the hierarchical depth buffer: the min (the farthest, depth is reversed) of the texels of the previous
 * level (or the depth texture) the texel covers: 2x2 of a level, up to 3x3 of the depth texture (level 0 maps the
 * viewport onto a power-of-two square, e.g. 1280 pixels onto 1024 texels: a texel's span can touch three pixels)
 */
export const hi_z_wgsl = /* wgsl */ `
${quad_vert_wgsl}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let size = vec2<i32>(textureDimensions(tPreviousLevel));
    let s = vec2<f32>(size);
    // the texel's span in the previous level's texels: those it touches (at most three a side)
    let start = ((in.position.xy - 0.5) * uObject.uInvSize + uObject.uRect.xy) * s;
    let end = ((in.position.xy + 0.5) * uObject.uInvSize + uObject.uRect.xy) * s;
    let lo = clamp(vec2<i32>(floor(start)), vec2<i32>(0), size - 1);
    let hi = clamp(vec2<i32>(ceil(end)) - 1, lo, min(lo + 2, size - 1));

    var depth = 1.0;
    for (var y = lo.y; y <= hi.y; y++) {
        for (var x = lo.x; x <= hi.x; x++) {
            depth = min(depth, textureLoad(tPreviousLevel, vec2<i32>(x, y), 0).r);
        }
    }
    return vec4<f32>(depth);
}
`;
