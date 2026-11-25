/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** The ghost of what the clip objects remove over the image, as `clipGhost_frag` */
export const clip_ghost_wgsl = /* wgsl */ `
${packing_wgsl}
${depth_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tGhostDepth')}
${texture_read_wgsl('tDepth')}
${texture_read_wgsl('tTransparentDepth')}

fn getGhostDepth(coords: vec2<f32>) -> f32 {
    return unpackRGBAToDepth(read_tGhostDepth(coords));
}

fn getSceneDepth(coords: vec2<f32>) -> f32 {
    #if defined(dPackedDepth)
        var depth = unpackRGBAToDepth(read_tDepth(coords));
    #else
        var depth = read_tDepth(coords).r;
    #endif
    if (uObject.uHasTransparentDepth != 0u) {
        let layer = unpackRGBAToDepthWithAlpha(read_tTransparentDepth(coords));
        // depth is reversed: nearer is larger
        if (layer.x > depth && 1.0 - layer.y < uObject.uTransparentThreshold) {
            depth = layer.x;
        }
    }
    return depth;
}

fn isEdge(coords: vec2<f32>, viewZ: f32) -> bool {
    let n = getGhostDepth(coords);
    if (n <= 0.0) {
        return true;
    }
    let nZ = depthToViewZ(uObject.uIsOrtho, n, uObject.uNear, uObject.uFar);
    return abs(nZ - viewZ) > uObject.uEdgeDistance * abs(viewZ);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = in.position.xy * uObject.uTexSizeInv;
    let depth = getGhostDepth(coords);
    // nothing removed here (cleared to 0.0, the far plane), or behind the scene (depth is reversed)
    if (depth <= 0.0 || depth < getSceneDepth(coords)) {
        discard;
    }
    let viewZ = depthToViewZ(uObject.uIsOrtho, depth, uObject.uNear, uObject.uFar);
    let step = uObject.uGhostTexSizeInv;
    let edge = isEdge(coords + vec2<f32>(step.x, 0.0), viewZ) || isEdge(coords - vec2<f32>(step.x, 0.0), viewZ)
        || isEdge(coords + vec2<f32>(0.0, step.y), viewZ) || isEdge(coords - vec2<f32>(0.0, step.y), viewZ);
    let alpha = select(uObject.uAlpha, uObject.uEdgeAlpha, edge);
    return vec4<f32>(uObject.uColor * alpha, alpha);
}
`;
