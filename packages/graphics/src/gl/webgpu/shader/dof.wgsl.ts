/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/**
 * Depth of field: the color blurred by the circle of confusion of the nearest depth, planar or spherical around the
 * focus, in GL window space; bounded by the \`oBlurSize\` override constant. Depth textures are always available on
 * WebGPU.
 */
export const dof_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tColor')}
${texture_read_wgsl('tDepthOpaque')}
${texture_read_wgsl('tDepthTransparent')}

override oBlurSize: i32 = 9;

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn getDepth(coords: vec2<f32>) -> f32 {
    let opaque = read_tDepthOpaque(coords).r;
    let transparent = unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords)).x;
    return min(opaque, transparent);
}

/** circle of confusion */
fn getCOC(uv: vec2<f32>) -> f32 {
    let depth = getDepth(uv);
    let viewDist = getViewZ(depth);
    let aposition = screenSpaceToViewSpace(vec3<f32>(uv, depth), uObject.uInvProjection);
    let focusDist = length(aposition - uObject.uCenter);
    var coc = 0.0;
    if (uObject.uMode == 0) {
        // planar depth of field
        coc = (abs(viewDist) - uObject.uInFocus) / uObject.uPPM;
    } else if (uObject.uMode == 1) {
        // spherical depth of field
        coc = focusDist / uObject.uPPM;
    }
    return clamp(coc, -1.0, 1.0);
}

/** box blur weighted by the circle of confusion */
fn getBlurredImage(coords: vec2<f32>) -> vec3<f32> {
    var blurColor = vec3<f32>(0.0);
    let texelSize = 1.0 / uObject.uTexSize;
    let blurSize = f32(oBlurSize);
    var count = 0.0;
    for (var x = 0; x < oBlurSize; x++) {
        for (var y = 0; y < oBlurSize; y++) {
            let offset = vec2<f32>(f32(x) - blurSize / 2.0, f32(y) - blurSize / 2.0);
            let uvPixel = coords + offset * texelSize * uObject.uBlurSpread;
            let coc = smoothstep(0.0, 1.0, abs(getCOC(uvPixel)));
            blurColor += read_tColor(uvPixel).rgb * coc;
            count += coc;
        }
    }
    return blurColor / count;
}

// simplification from https://catlikecoding.com/unity/tutorials/advanced-rendering/depth-of-field/
@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let uv = glCoords(in.position, uObject.uTexSize);
    let color = read_tColor(uv);
    let blurColor = getBlurredImage(uv);
    let coc = getCOC(uv);
    return vec4<f32>(mix(color.rgb, blurColor, smoothstep(0.0, 1.0, abs(coc))), color.a);
}
`;
