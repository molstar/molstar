/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** Weighted sample of the multi-sample accumulation */
export const multi_sample_compose_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tColor')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, uObject.uTexSize);
    return read_tColor(coords) * uObject.uWeight;
}
`;
