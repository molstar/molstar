/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** Marking mask edges, in GL window space as in GLSL */
export const edge_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tMaskTexture')}

// the mask is antialiased, so the coverage can be taken from a single texel
fn getMask(coords: vec2<f32>) -> vec4<f32> {
    return read_tMaskTexture(coords);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let texSizeInv = uObject.uTexSizeInv;
    let coords = glCoords(in.position, 1.0 / texSizeInv);
    let edgeScale = uObject.uEdgeScale;
    let offset = vec4<f32>(edgeScale, 0.0, 0.0, edgeScale) * vec4<f32>(texSizeInv, texSizeInv);

    let m1 = getMask(coords + offset.xy);
    let m2 = getMask(coords - offset.xy);
    let m3 = getMask(coords + offset.yw);
    let m4 = getMask(coords - offset.yw);

    let edge = max(abs(m1.r - m2.r), abs(m3.r - m4.r));
    if (edge <= 0.0) {
        discard;
    }

    // the sample covering the most of the marked region, mask values are 0.0 where marked
    let mx = select(m2, m1, m1.r < m2.r);
    let my = select(m4, m3, m3.r < m4.r);
    let m = select(my, mx, mx.r < my.r);

    // the samples average marked and unmarked texels and unmarked ones are 1.0 in every channel,
    // so the values of the marked texels can be recovered from how much they are covered
    let marked = clamp((m.gba - m.r) / max(1.0 - m.r, 0.001), vec3<f32>(0.0), vec3<f32>(1.0));

    let visibility = select(0.0, 1.0, marked.x > 0.5);
    let mask = read_tMaskTexture(coords).r;
    let marker = select(0.0, 1.0, marked.y > 0.5);
    let fogAlpha = marked.z;
    // premultiplied, so that a blurred sample of this can be decoded the same way
    return vec4<f32>(visibility, mask, marker, 1.0) * (edge * fogAlpha);
}
`;

/** Coverage of the visible unmarked geometry, in GL window space as in GLSL */
export const dim_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tDepthTexture')}
${texture_read_wgsl('tMaskTexture')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = glCoords(in.position, 1.0 / uObject.uTexSizeInv);

    // depth of the front-most unmarked fragment, cleared to 1.0 in every channel where there is none
    let packedDepth = read_tDepthTexture(coords);
    if (all(packedDepth == vec4<f32>(1.0))) {
        discard;
    }

    // mask values are 0.0 where marked, green is 0.0 where the marked fragment is in front
    let m = read_tMaskTexture(coords);
    if (m.r < 0.5 && m.g < 0.5) {
        discard;
    }

    let depthWithAlpha = unpackRGBAToDepthWithAlpha(packedDepth);
    let viewZ = depthToViewZ(uObject.uIsOrtho, depthWithAlpha.x, uObject.uNear, uObject.uFar);
    let fogAlpha = 1.0 - smoothstep(uObject.uFogNear, uObject.uFogFar, abs(viewZ));
    // after antialiasing: r = coverage * opacity * fogAlpha, g = coverage * opacity
    return vec4<f32>(fogAlpha * depthWithAlpha.y, depthWithAlpha.y, 0.0, 1.0);
}
`;

/** Marking edges, fill and dimming of one sample, in GL window space as in GLSL */
export const overlay_wgsl = /* wgsl */ `
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tEdgeTexture')}
${texture_read_wgsl('tMaskTexture')}
${texture_read_wgsl('tDimTexture')}
${texture_read_wgsl('tSsaoDepth')}
${texture_read_wgsl('tBase')}
${texture_read_wgsl('tShaded')}

fn unpackMarkingOpacity(g: f32) -> f32 {
    return clamp(1.0 - 4.0 * min(g, 1.0 - g), 0.0, 1.0);
}

// mixes toward the SSAO occlusion color, like the postprocessing pass
fn applyOcclusion(color: vec3<f32>, occlusion: f32, fogAlpha: f32) -> vec3<f32> {
    var occlusionColor = mix(uObject.uFogColor, uObject.uOcclusionColor, fogAlpha);
    if (uObject.uTransparentBackground != 0u) {
        occlusionColor = uObject.uOcclusionColor * fogAlpha;
    }
    return mix(occlusionColor, color, occlusion);
}

fn luminance(c: vec3<f32>) -> f32 {
    // https://www.w3.org/TR/WCAG21/#dfn-relative-luminance
    return dot(c, vec3<f32>(0.2125, 0.7154, 0.0721));
}

// the path traced occlusion and shadows, i.e. the base image relative to the direct shading
fn getTracedShading(baseLuminance: f32, shadedLuminance: f32, fogAlpha: f32) -> f32 {
    let a = max(fogAlpha, 0.001);
    // undo the fog of the base image
    var l = (baseLuminance - luminance(uObject.uFogColor) * (1.0 - fogAlpha)) / a;
    if (uObject.uTransparentBackground != 0u) {
        l = baseLuminance / a;
    }
    return clamp(l / max(shadedLuminance, 0.001), 0.0, 1.0);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = glCoords(in.position, 1.0 / uObject.uTexSizeInv);

    // coverage of the visible unmarked geometry, see dim
    var dimAlpha = 0.0;
    var dimFogAlpha = 1.0;
    if (uObject.uDimStrength > 0.0) {
        let d = read_tDimTexture(coords);
        dimAlpha = d.r * uObject.uDimStrength;
        dimFogAlpha = d.r / max(d.g, 0.001);
    }

    // solid fill covering the interior of marked regions, sampled directly from the mask
    var fillRgb = vec3<f32>(0.0);
    var fillAlpha = 0.0;
    var fillFogAlpha = 1.0;
    if (uObject.uHighlightFillStrength > 0.0 || uObject.uSelectFillStrength > 0.0) {
        let m = read_tMaskTexture(coords);
        let coverage = 1.0 - m.r;
        if (coverage > 0.0) {
            let marked = clamp((m.gba - m.r) / max(coverage, 0.001), vec3<f32>(0.0), vec3<f32>(1.0));
            let isHighlight = marked.y > 0.5;
            fillRgb = select(uObject.uSelectFillColor, uObject.uHighlightFillColor, isHighlight);
            let fillStrength = select(uObject.uSelectFillStrength, uObject.uHighlightFillStrength, isHighlight);
            fillFogAlpha = marked.z;
            // marked.x: 1.0 = hidden, 0.0 = visible; without depth test every marked texel reads as hidden
            let visibility = select(1.0, 0.0, uObject.uDepthTest != 0u && marked.x > 0.5);
            fillAlpha = coverage * fillStrength * visibility * fillFogAlpha * unpackMarkingOpacity(marked.x);
        }
    }

    // keep the shading of the image, so that blending toward the color does not wash it out
    var dimRgb = uObject.uDimColor;
    #if defined(dMarkingShading_ssao)
        if (dimAlpha > 0.0 || fillAlpha > 0.0) {
            let occlusion = unpackSsao(read_tSsaoDepth(coords + uObject.uOcclusionOffset));
            dimRgb = applyOcclusion(dimRgb, occlusion, dimFogAlpha);
            fillRgb = applyOcclusion(fillRgb, occlusion, fillFogAlpha);
        }
    #elif defined(dMarkingShading_traced)
        if (dimAlpha > 0.0 || fillAlpha > 0.0) {
            let baseLuminance = luminance(read_tBase(coords).rgb);
            let shadedLuminance = luminance(read_tShaded(coords).rgb);
            dimRgb *= getTracedShading(baseLuminance, shadedLuminance, dimFogAlpha);
            fillRgb *= getTracedShading(baseLuminance, shadedLuminance, fillFogAlpha);
        }
    #endif

    let e = read_tEdgeTexture(coords);
    var edgeRgb = vec3<f32>(0.0);
    var edgeAlpha = 0.0;
    if (e.a > 0.0) {
        // the edge is stored premultiplied, see edge
        let edgeValue = e.rgb / e.a;
        let isHighlight = edgeValue.z > 0.5;
        edgeRgb = select(uObject.uSelectEdgeColor, uObject.uHighlightEdgeColor, isHighlight);
        if (edgeValue.y <= 0.5) {
            edgeRgb *= uObject.uInnerEdgeFactor;
        }
        edgeAlpha = select(1.0, uObject.uGhostEdgeStrength, edgeValue.x > 0.5) * e.a;
        edgeAlpha *= select(uObject.uSelectEdgeStrength, uObject.uHighlightEdgeStrength, isHighlight);
    }

    // premultiplied, so that partially covered pixels are not weighted by alpha twice; dim, fill, edge from bottom to top
    let rgb = fillRgb * fillAlpha + dimRgb * dimAlpha * (1.0 - fillAlpha);
    let alpha = fillAlpha + dimAlpha * (1.0 - fillAlpha);
    return vec4<f32>(edgeRgb * edgeAlpha + rgb * (1.0 - edgeAlpha), edgeAlpha + alpha * (1.0 - edgeAlpha));
}
`;
