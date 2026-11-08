/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';
import type { DefineLookup } from '../../renderable/schema.js';

/** Coverage of the visible unmarked geometry */
export const dim_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tDepthTexture')}
${texture_read_wgsl('tMaskTexture')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, 1.0 / uObject.uTexSizeInv);

    // depth of the front-most unmarked fragment, cleared to 0.0 (the far plane) in every channel where there is none
    let packedDepth = read_tDepthTexture(coords);
    if (all(packedDepth == vec4<f32>(0.0))) {
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

/**
 * The marking depth from the scene's: its opaque depth, or its front transparent layer where that is nearer and less
 * transparent than \`uTransparentThreshold\`; packed with the opacity as the marking depth pass writes it, cleared to 0.0
 * (the far plane) where there is nothing
 */
export const marking_depth_wgsl = /* wgsl */ `
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tDepth')}
${texture_read_wgsl('tTransparentDepth')}

fn getDepth(coords: vec2<f32>) -> vec2<f32> {
    var depth = read_tDepth(coords).r;
    var alpha = 1.0;
    if (uObject.uHasTransparentDepth != 0u) {
        let layer = unpackRGBAToDepthWithAlpha(read_tTransparentDepth(coords));
        // depth is reversed: nearer is larger
        if (layer.x > depth && 1.0 - layer.y < uObject.uTransparentThreshold) {
            depth = layer.x;
            alpha = layer.y;
        }
    }
    return vec2<f32>(depth, alpha);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    var d: vec2<f32>;
    if (all(uObject.uShift == vec2<f32>(0.0))) {
        d = getDepth(in.position.xy * uObject.uTexSizeInv);
    } else {
        // the centers of the four texels around the point
        let p = floor(in.position.xy + uObject.uShift - 0.5) + 0.5;
        d = getDepth(p * uObject.uTexSizeInv);
        let d1 = getDepth((p + vec2<f32>(1.0, 0.0)) * uObject.uTexSizeInv);
        if (d1.x < d.x) { d = d1; }
        let d2 = getDepth((p + vec2<f32>(0.0, 1.0)) * uObject.uTexSizeInv);
        if (d2.x < d.x) { d = d2; }
        let d3 = getDepth((p + vec2<f32>(1.0, 1.0)) * uObject.uTexSizeInv);
        if (d3.x < d.x) { d = d3; }
    }
    if (d.x <= 0.0) {
        discard;
    }
    return packDepthWithAlphaToRGBA(d.x, d.y);
}
`;

/** The SSAO shading of the overlay (`dMarkingShading_ssao`): the occlusion read with or without upsampling */
const overlay_ssao_wgsl = (d: DefineLookup) => /* wgsl */ `
${d('dOcclusionUpsample') ? texture_read_wgsl('tDepth') : `${texture_read_wgsl('tSsaoDepth')}\n${texture_read_wgsl('tSsaoDepthTransparent')}`}

/** occlusion at \`coords\`, depth-aware upsampled from a lower SSAO resolution with the scene's opaque depth */
fn getOcclusion(coords: vec2<f32>) -> f32 {
    #if defined(dOcclusionUpsample)
        let viewZ = depthToViewZ(uObject.uIsOrtho, read_tDepth(coords).r, uObject.uNear, uObject.uFar);
        return getSsaoUpsampled(tSsaoDepth, coords + uObject.uOcclusionOffset, viewZ, uObject.uNear, uObject.uFar);
    #else
        return unpackSsao(read_tSsaoDepth(coords + uObject.uOcclusionOffset));
    #endif
}

/** occlusion of the front transparent layer at \`coords\`, at its \`depth\` (0: none there) */
fn getOcclusionTransparent(coords: vec2<f32>, depth: f32) -> f32 {
    if (uObject.uHasSsaoTransparent == 0u || depth == 0.0) {
        return 1.0;
    }
    #if defined(dOcclusionUpsample)
        let viewZ = depthToViewZ(uObject.uIsOrtho, depth, uObject.uNear, uObject.uFar);
        return getSsaoUpsampled(tSsaoDepthTransparent, coords + uObject.uOcclusionOffset, viewZ, uObject.uNear, uObject.uFar);
    #else
        return unpackSsao(read_tSsaoDepthTransparent(coords + uObject.uOcclusionOffset));
    #endif
}

// mixes toward the SSAO occlusion color, like the postprocessing pass
fn applyOcclusion(color: vec3<f32>, occlusion: f32, fogAlpha: f32) -> vec3<f32> {
    var occlusionColor = mix(uObject.uFogColor, uObject.uOcclusionColor, fogAlpha);
    if (uObject.uTransparentBackground != 0u) {
        occlusionColor = uObject.uOcclusionColor * fogAlpha;
    }
    return mix(occlusionColor, color, occlusion);
}
`;

/** The traced shading of the overlay (`dMarkingShading_traced`): the traced image relative to the direct shading */
const overlay_traced_wgsl = /* wgsl */ `
${texture_read_wgsl('tBase')}
${texture_read_wgsl('tShaded')}

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
`;

/** Marking edges, fill and dimming of one sample; the shading's reads and helpers by the variant's defines */
export const overlay_wgsl = (d: DefineLookup) => /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tMaskTexture')}
${texture_read_wgsl('tDimTexture')}
${texture_read_wgsl('tTransparentDepth')}
${d('dMarkingShading_ssao') ? overlay_ssao_wgsl(d) : ''}
${d('dMarkingShading_traced') ? overlay_traced_wgsl : ''}

// edge of the marked regions at \`coords\` (alpha 0 if none): the mask is antialiased, so the coverage can be taken from
// single texels, \`uEdgeScale\` texels apart; rgb: hidden, outer (the mask), highlight
fn getEdge(coords: vec2<f32>) -> vec4<f32> {
    let texSizeInv = uObject.uTexSizeInv;
    let edgeScale = uObject.uEdgeScale;
    let offset = vec4<f32>(edgeScale, 0.0, 0.0, edgeScale) * vec4<f32>(texSizeInv, texSizeInv);

    let m1 = read_tMaskTexture(coords + offset.xy);
    let m2 = read_tMaskTexture(coords - offset.xy);
    let m3 = read_tMaskTexture(coords + offset.yw);
    let m4 = read_tMaskTexture(coords - offset.yw);

    let edge = max(abs(m1.r - m2.r), abs(m3.r - m4.r));
    if (edge <= 0.0) {
        return vec4<f32>(0.0);
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
    return vec4<f32>(visibility, mask, marker, edge * fogAlpha);
}

fn unpackMarkingOpacity(g: f32) -> f32 {
    return clamp(1.0 - 4.0 * min(g, 1.0 - g), 0.0, 1.0);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, 1.0 / uObject.uTexSizeInv);

    // coverage of the visible unmarked geometry, see dim
    var dimCoverage = 0.0;
    var dimFogAlpha = 1.0;
    if (uObject.uDimStrength > 0.0) {
        let d = read_tDimTexture(coords);
        dimCoverage = d.r;
        dimFogAlpha = d.r / max(d.g, 0.001);
    }

    // solid fill covering the interior of marked regions, sampled directly from the mask
    var fillRgb = vec3<f32>(0.0);
    var fillAlpha = 0.0;
    var fillOpacity = 1.0;
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
            fillOpacity = unpackMarkingOpacity(marked.x);
            fillAlpha = coverage * fillStrength * visibility * fillFogAlpha * fillOpacity;
        }
    }

    // the front transparent layer, cleared to 0.0 (no alpha) where there is none: a fill below full opacity is of a
    // transparent fragment, the fill of opaque geometry is seen through the layer; an unmarked layer is dimmed by its
    // alpha, the opaque share behind it by what it lets through
    var layerAlpha = 0.0;
    var layerDepth = 0.0;
    if (uObject.uHasTransparentDepth != 0u && (uObject.uDimStrength > 0.0 || fillAlpha > 0.0)) {
        let layer = unpackRGBAToDepthWithAlpha(read_tTransparentDepth(coords));
        layerDepth = layer.x;
        layerAlpha = layer.y;
    }
    let fillOfLayer = fillOpacity < 0.98;
    let dimLayer = select(0.0, layerAlpha, uObject.uDimStrength > 0.0 && !(fillAlpha > 0.0 && fillOfLayer));
    let dimOpaque = (1.0 - dimLayer) * dimCoverage;
    let dimWeight = dimLayer + dimOpaque;
    let dimAlpha = uObject.uDimStrength * dimWeight;
    if (!fillOfLayer) {
        fillAlpha *= 1.0 - layerAlpha;
    }

    // keep the shading of the image, so that blending toward the color does not wash it out
    var dimRgb = uObject.uDimColor;
    #if defined(dMarkingShading_ssao)
        if (dimAlpha > 0.0 || fillAlpha > 0.0) {
            let occlusion = getOcclusion(coords);
            var occlusionLayer = 1.0;
            if (dimLayer > 0.0 || fillOfLayer) {
                occlusionLayer = getOcclusionTransparent(coords, layerDepth);
            }
            var dimOcclusion = occlusion;
            if (dimWeight > 0.0) {
                dimOcclusion = (occlusionLayer * dimLayer + occlusion * dimOpaque) / dimWeight;
            }
            dimRgb = applyOcclusion(dimRgb, dimOcclusion, dimFogAlpha);
            fillRgb = applyOcclusion(fillRgb, select(occlusion, occlusionLayer, fillOfLayer), fillFogAlpha);
        }
    #elif defined(dMarkingShading_traced)
        if (dimAlpha > 0.0 || fillAlpha > 0.0) {
            let baseLuminance = luminance(read_tBase(coords).rgb);
            let shadedLuminance = luminance(read_tShaded(coords).rgb);
            dimRgb *= getTracedShading(baseLuminance, shadedLuminance, dimFogAlpha);
            fillRgb *= getTracedShading(baseLuminance, shadedLuminance, fillFogAlpha);
        }
    #endif

    let e = getEdge(coords);
    var edgeRgb = vec3<f32>(0.0);
    var edgeAlpha = 0.0;
    if (e.a > 0.0) {
        let edgeValue = e.rgb;
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
