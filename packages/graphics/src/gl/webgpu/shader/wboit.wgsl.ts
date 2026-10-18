/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** Resolves the weighted blended transparency accumulation */
export const evaluate_wboit_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tWboitA')}
${texture_read_wgsl('tWboitB')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, uObject.uTexSize);

    let accum = read_tWboitA(coords);
    let r = 1.0 - accum.a;

    let weight = read_tWboitB(coords).r;
    // divisor needs to allow very small values for nice fading
    return vec4<f32>(accum.rgb / clamp(weight, 0.00000001, 50000.0), r);
}
`;
