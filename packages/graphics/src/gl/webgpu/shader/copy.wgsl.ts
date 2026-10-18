/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** Copies `tColor`, filtered as on WebGL when the source texture is linear (e.g. depth downsampling) */
export const copy_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tColor')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    return read_tColor(fragmentUv(in.position, uObject.uTexSize));
}
`;
