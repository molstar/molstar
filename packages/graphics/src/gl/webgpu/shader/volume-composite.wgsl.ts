/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/**
 * Direct volumes rendered at a lower resolution (`tVolume`, premultiplied) upsampled to the target: the four nearest
 * texels weighted bilinearly and by how close the opaque depth their rays stopped at is to the pixel's, so that the
 * volume does not bleed across the edges of opaque objects
 */
export const volume_composite_wgsl = /* wgsl */ `
${depth_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tDepth')}
${texture_read_wgsl('tVolume')}

fn viewZAt(pixel: vec2<f32>) -> f32 {
    let depth = read_tDepth((pixel + 0.5) / uObject.uTexSize).r;
    return depthToViewZ(uObject.uIsOrtho, depth, uObject.uNear, uObject.uFar);
}

fn volumeTexel(texelIn: vec2<f32>, bilinear: f32, viewZ: f32, weightSum: ptr<function, f32>) -> vec4<f32> {
    let volumeSize = uObject.uVolumeSize;
    let texel = clamp(texelIn, vec2<f32>(0.0), volumeSize - 1.0);
    // the pixel whose opaque depth the texel's ray read (its position over the drawing buffer size at the lower size)
    let texelViewZ = viewZAt(floor((texel + 0.5) / volumeSize * uObject.uTexSize));
    let weight = bilinear / (0.001 + abs(texelViewZ - viewZ));
    *weightSum += weight;
    return read_tVolume((texel + 0.5) / volumeSize) * weight;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let viewZ = viewZAt(floor(in.position.xy));
    let p = in.position.xy / uObject.uTexSize * uObject.uVolumeSize - 0.5;
    let i = floor(p);
    let f = p - i;
    var weightSum = 0.0;
    var sum = volumeTexel(i, (1.0 - f.x) * (1.0 - f.y), viewZ, &weightSum);
    sum += volumeTexel(i + vec2<f32>(1.0, 0.0), f.x * (1.0 - f.y), viewZ, &weightSum);
    sum += volumeTexel(i + vec2<f32>(0.0, 1.0), (1.0 - f.x) * f.y, viewZ, &weightSum);
    sum += volumeTexel(i + vec2<f32>(1.0, 1.0), f.x * f.y, viewZ, &weightSum);
    if (weightSum > 0.0) {
        return sum / weightSum;
    }
    return vec4<f32>(0.0);
}
`;
