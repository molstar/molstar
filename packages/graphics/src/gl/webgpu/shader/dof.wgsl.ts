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
 * focus; bounded by the \`oBlurSize\` override constant. Depth textures are always available on
 * WebGPU.
 *
 * The blur is a box weighted per sample, so it separates exactly: \`dDofPass\` 'horizontal' writes the weighted sums of
 * a row (rgb) and their weights (a), 'vertical' adds those of the column and mixes; 'single' does both in one pass
 */
export const dof_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tColor')}
${texture_read_wgsl('tDepthOpaque')}
${texture_read_wgsl('tDepthTransparent')}
#if defined(dDofPass_vertical)
    ${texture_read_wgsl('tBlur')}
#endif

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
    return max(opaque, transparent); // the nearer (depth is reversed)
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

/** the color weighted by the circle of confusion at \`uv\` (rgb), and the weight (a) */
fn getWeightedColor(uv: vec2<f32>) -> vec4<f32> {
    let coc = smoothstep(0.0, 1.0, abs(getCOC(uv)));
    return vec4<f32>(read_tColor(uv).rgb * coc, coc);
}

/**
 * the offset of sample \`i\` of a row or column, in texels, centered (texel centers for odd sizes, not their edges,
 * where rounding picks a neighbour); y up, as in the GLSL program
 */
fn sampleOffset(i: i32) -> f32 {
    return f32(i) - f32(oBlurSize - 1) / 2.0;
}

/** box blur weighted by the circle of confusion: weighted color sums (rgb) and weights (a) */
fn getBlurSums(coords: vec2<f32>) -> vec4<f32> {
    let texelSize = 1.0 / uObject.uTexSize;
    var sums = vec4<f32>(0.0);
    #if defined(dDofPass_horizontal)
        for (var x = 0; x < oBlurSize; x++) {
            sums += getWeightedColor(coords + vec2<f32>(sampleOffset(x), 0.0) * texelSize * uObject.uBlurSpread);
        }
    #elif defined(dDofPass_vertical)
        for (var y = 0; y < oBlurSize; y++) {
            sums += read_tBlur(coords + vec2<f32>(0.0, -sampleOffset(y)) * texelSize * uObject.uBlurSpread);
        }
    #else
        for (var x = 0; x < oBlurSize; x++) {
            for (var y = 0; y < oBlurSize; y++) {
                sums += getWeightedColor(coords + vec2<f32>(sampleOffset(x), -sampleOffset(y)) * texelSize * uObject.uBlurSpread);
            }
        }
    #endif
    return sums;
}

// simplification from https://catlikecoding.com/unity/tutorials/advanced-rendering/depth-of-field/
@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let uv = fragmentUv(in.position, uObject.uTexSize);
    #if defined(dDofPass_horizontal)
        return getBlurSums(uv);
    #else
        let color = read_tColor(uv);
        let sums = getBlurSums(uv);
        // no sample out of focus: the pixel itself
        let blurColor = select(color.rgb, sums.rgb / sums.a, sums.a > 0.0);
        let coc = getCOC(uv);
        return vec4<f32>(mix(color.rgb, blurColor, smoothstep(0.0, 1.0, abs(coc))), color.a);
    #endif
}
`;
