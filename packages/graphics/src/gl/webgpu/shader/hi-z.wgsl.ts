/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl } from './chunks/quad.wgsl.js';

/** A level of the hierarchical depth buffer: the min of 2x2 texels (the farthest, depth is reversed) of the previous level (or the depth texture) */
export const hi_z_wgsl = /* wgsl */ `
${quad_vert_wgsl}

/** depth of the texel at \`p\`, clamped to the edge as the GLSL \`textureOffset\` reads */
fn load_tPreviousLevel(p: vec2<i32>, size: vec2<i32>) -> f32 {
    return textureLoad(tPreviousLevel, clamp(p, vec2<i32>(0), size - 1), 0).r;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let position = in.position.xy * uObject.uInvSize + uObject.uRect.xy;

    let size = vec2<i32>(textureDimensions(tPreviousLevel));
    let p = vec2<i32>(floor(position * vec2<f32>(size)));

    let x = load_tPreviousLevel(p, size);
    let y = load_tPreviousLevel(p + vec2<i32>(-1, 0), size);
    // the texels left and below, as the GLSL reads
    let z = load_tPreviousLevel(p + vec2<i32>(-1, 1), size);
    let w = load_tPreviousLevel(p + vec2<i32>(0, 1), size);

    return vec4<f32>(min(min(x, y), min(z, w)));
}
`;
