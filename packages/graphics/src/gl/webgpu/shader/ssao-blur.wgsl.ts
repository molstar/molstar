/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** Depth-aware separable blur of the packed SSAO occlusion */
export const ssao_blur_wgsl = /* wgsl */ `
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tSsaoDepth')}

override oOcclusionKernelSize: i32 = 15;

// the packed depth is linear and reversed (see the SSAO pass)
fn getViewZ(depth: f32) -> f32 {
    return depth * (uObject.uFar - uObject.uNear) - uObject.uFar;
}

// the background (far plane) at 0, packed as (0, 0); the near plane at 1, which
// packUnitIntervalToRG stores as (0, 1) in the uint8 target and unpacks to 255/256, never 1.0, so the packed bytes
// are tested
fn isBackgroundPacked(packedDepth: vec2<f32>) -> bool {
    return packedDepth.x == 0.0 && packedDepth.y == 0.0;
}

fn isNearClipPacked(packedDepth: vec2<f32>) -> bool {
    return packedDepth.x == 0.0 && packedDepth.y == 1.0;
}

fn outsideBounds(p: vec2<f32>) -> bool {
    let b = uObject.uBounds;
    return p.x < b.x || p.y < b.y || p.x > b.z || p.y > b.w;
}

/** view space size of a pixel at view z */
fn getPixelSize(viewZ: f32) -> f32 {
    #if dOrthographic == 1
        return 2.0 * uObject.uInvProjection[0][0] / uObject.uTexSize.x;
    #else
        return 2.0 * uObject.uInvProjection[0][0] * -viewZ / uObject.uTexSize.x;
    #endif
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, uObject.uTexSize);
    let selfSsaoDepth = read_tSsaoDepth(coords);
    let packedDepth = selfSsaoDepth.zw;

    if (outsideBounds(coords)) {
        return vec4<f32>(packUnitIntervalToRG(1.0), packedDepth);
    }

    let selfDepth = unpackRGToUnitInterval(packedDepth);
    // (if background and if second pass) or if near clip
    if ((isBackgroundPacked(packedDepth) && uObject.uBlurDirectionY != 0.0) || isNearClipPacked(packedDepth)) {
        return vec4<f32>(packUnitIntervalToRG(1.0), packedDepth);
    }

    let selfViewZ = getViewZ(selfDepth);
    let pixelSize = getPixelSize(selfViewZ);

    let offset = vec2<f32>(uObject.uBlurDirectionX, uObject.uBlurDirectionY) / uObject.uTexSize;

    var sum = 0.0;
    var kernelSum = 0.0;
    let halfKernelSize = oOcclusionKernelSize / 2;
    // only if kernelSize is odd
    for (var i = -halfKernelSize; i <= halfKernelSize; i++) {
        let fi = f32(i);
        if (abs(fi) > 1.0 && abs(fi) * pixelSize > 0.8) {
            continue;
        }

        let sampleCoords = coords + fi * offset;
        if (outsideBounds(sampleCoords)) {
            continue;
        }

        let sampleSsaoDepth = read_tSsaoDepth(sampleCoords);
        let sampleDepth = unpackRGToUnitInterval(sampleSsaoDepth.zw);
        if (isBackgroundPacked(sampleSsaoDepth.zw) || isNearClipPacked(sampleSsaoDepth.zw)) {
            continue;
        }

        let sampleViewZ = getViewZ(sampleDepth);
        if (abs(selfViewZ - sampleViewZ) >= uObject.uBlurDepthBias) {
            continue;
        }

        let kernel = uObject.uKernel[abs(i)].x;
        let sampleValue = unpackRGToUnitInterval(sampleSsaoDepth.xy);

        sum += kernel * sampleValue;
        kernelSum += kernel;
    }
    // no sample within the depth bias (a background pixel of the first pass): the pixel unblurred
    if (kernelSum == 0.0) {
        return selfSsaoDepth;
    }
    return vec4<f32>(packUnitIntervalToRG(sum / kernelSum), packedDepth);
}
`;
