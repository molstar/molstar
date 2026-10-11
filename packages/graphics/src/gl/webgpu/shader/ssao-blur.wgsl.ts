/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** Depth-aware separable blur of the packed SSAO occlusion, in GL window space as in GLSL */
export const ssao_blur_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tSsaoDepth')}

override oOcclusionKernelSize: i32 = 15;

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

// packUnitIntervalToRG(1.0) stores (0, 1) in the uint8 target and unpacks to 255/256,
// never 1.0, so the packed bytes have to be tested instead
fn isBackgroundPacked(packedDepth: vec2<f32>) -> bool {
    return packedDepth.x == 0.0 && packedDepth.y == 1.0;
}

fn isNearClip(depth: f32) -> bool {
    return depth == 0.0;
}

fn outsideBounds(p: vec2<f32>) -> bool {
    let b = uObject.uBounds;
    return p.x < b.x || p.y < b.y || p.x > b.z || p.y > b.w;
}

fn getPixelSize(coords: vec2<f32>, depth: f32) -> f32 {
    let viewPos0 = screenSpaceToViewSpace(vec3<f32>(coords, depth), uObject.uInvProjection);
    let viewPos1 = screenSpaceToViewSpace(vec3<f32>(coords + vec2<f32>(1.0, 0.0) / uObject.uTexSize, depth), uObject.uInvProjection);
    return distance(viewPos0, viewPos1);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = glCoords(in.position, uObject.uTexSize);
    let packedDepth = read_tSsaoDepth(coords).zw;

    if (outsideBounds(coords)) {
        return vec4<f32>(packUnitIntervalToRG(1.0), packedDepth);
    }

    let selfDepth = unpackRGToUnitInterval(packedDepth);
    // (if background and if second pass) or if near clip
    if ((isBackgroundPacked(packedDepth) && uObject.uBlurDirectionY != 0.0) || isNearClip(selfDepth)) {
        return vec4<f32>(packUnitIntervalToRG(1.0), packedDepth);
    }

    let selfViewZ = getViewZ(selfDepth);
    let pixelSize = getPixelSize(coords, selfDepth);

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
        if (isBackgroundPacked(sampleSsaoDepth.zw) || isNearClip(sampleDepth)) {
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
    return vec4<f32>(packUnitIntervalToRG(sum / kernelSum), packedDepth);
}
`;
