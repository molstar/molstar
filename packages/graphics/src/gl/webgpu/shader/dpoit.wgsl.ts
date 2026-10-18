/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** The front color of the dual depth peeling layers */
export const evaluate_dpoit_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tDpoitFrontColor')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    return read_tDpoitFrontColor(fragmentUv(in.position, uObject.uTexSize));
}
`;

/** Blends the back color of a dual depth peeling layer */
export const blend_back_dpoit_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tDpoitBackColor')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let color = read_tDpoitBackColor(fragmentUv(in.position, uObject.uTexSize));
    if (color.a == 0.0) {
        discard;
    }
    return color;
}
`;
