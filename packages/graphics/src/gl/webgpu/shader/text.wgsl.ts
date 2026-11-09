/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { clip_flag_wgsl, clip_functions_wgsl, clip_instance_wgsl, clip_pixel_wgsl } from './chunks/clip.wgsl.js';
import { segment_wgsl } from './chunks/segment.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import {
  apply_fog_wgsl,
  apply_marker_wgsl,
  assign_material_wgsl,
  apply_occlusion_wgsl,
  check_transparency_wgsl,
  material_functions_wgsl,
} from './chunks/material.wgsl.js';
import {
  depth_output_wgsl,
  dpoit_output_wgsl,
  emissive_output_wgsl,
  fragment_output_wgsl,
  marking_output_wgsl,
  pick_output_wgsl,
  tracing_output_wgsl,
  wboit_output_wgsl,
} from './chunks/output.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { texture_sample_wgsl } from './chunks/quad.wgsl.js';
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';
import { assign_size_wgsl, assign_varyings_wgsl, vertex_output_wgsl } from './chunks/varyings.wgsl.js';

/** Text as glyph quads (and background quads) facing the camera, glyphs from the SDF font atlas `tFont`, unlit */
export const text_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}
${texture_sample_wgsl('tFont')}

${vertex_output_wgsl('    @location(9) vTexCoord: vec2<f32>,')}
${fragment_output_wgsl('    @builtin(frag_depth) fragDepth: f32,')}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    // the item of the vertex (4 texels of \`tText\`) and its corner (4 vertices per item)
    let item = i32(in.vertexIndex) / 4;
    let corner = i32(in.vertexIndex) - item * 4;
    // the data vertex ('vertex' color, size): one per item
    let vertexId = item;
    let textPosition = read_tText(item * 4); // position xyz, group
    let textFlags = read_tText(item * 4 + 1); // depth, tether
    let textA = read_tText(item * 4 + 2);
    let textB = read_tText(item * 4 + 3);
    let group = textPosition.w;
    var mapping: vec2<f32>;
    var texCoord: vec2<f32>;
    if (textFlags.y > 0.5) {
        // a tether: its four corners, no glyph
        mapping = select(select(select(textB.zw, textB.xy, corner == 2), textA.zw, corner == 1), textA.xy, corner == 0);
        texCoord = vec2<f32>(10.0);
    } else {
        // a glyph or background: top left, bottom left, top right, bottom right of its rectangles
        let right = corner >= 2;
        let top = corner == 0 || corner == 2;
        mapping = vec2<f32>(select(textA.x, textA.z, right), select(textA.w, textA.y, top));
        texCoord = vec2<f32>(select(textB.x, textB.z, right), select(textB.w, textB.y, top));
    }
    ${assign_varyings_wgsl}
    ${assign_size_wgsl}

    out.vTexCoord = texCoord;

    let scale = uGlobal.uModelScale;
    let offsetX = uMaterial.uOffsetX * scale;
    let offsetY = uMaterial.uOffsetY * scale;
    var offsetZ = (uMaterial.uOffsetZ + textFlags.x * 0.95) * scale;

    let aTransform = affineTransform(in.aTransform0, in.aTransform1, in.aTransform2);
    let transform = uObject.uObjectTransform * aTransform;
    let position4 = vec4<f32>(textPosition.xyz, 1.0);
    let hasEyeCamera = uGlobal.uHasEyeCamera != 0u;
    var mvPosition = uGlobal.uModelView * transform * position4;
    if (hasEyeCamera) {
        mvPosition = uGlobal.uModelViewEye * transform * position4;
    }

    out.vModelPosition = (uGlobal.uModel * transform * position4).xyz;

    var mvCenter = vec4<f32>(mvPosition.xyz, 1.0);

    if (texCoord.x == 10.0) { // indicates background plane
        // move a bit to the back, taking distance to camera into account to avoid z-fighting; the distance is to the
        // GL clip-space position (z of the [-1, 1] projection: 2 * z - w of the [0, 1] one)
        let p = uGlobal.uProjection * mvCenter;
        offsetZ -= 0.001 * distance(uGlobal.uCameraPosition, vec3<f32>(p.xy, 2.0 * p.z - p.w));
    }

    // apply Z offset in view space
    if (!hasEyeCamera) {
        if (uGlobal.uIsOrtho == 1.0) {
            mvCenter.z += offsetZ;
        } else {
            mvCenter = vec4<f32>(mvCenter.xyz + normalize(-mvCenter.xyz) * offsetZ, mvCenter.w);
        }
    }

    if (hasEyeCamera) {
        mvCenter = uGlobal.uModelView * uGlobal.uInvModelViewEye * mvCenter;
    }

    // project center to clip space
    var clip = uGlobal.uProjection * mvCenter;

    // compute corner offset in screen-space units
    var cornerOffset = mapping * size * scale;
    cornerOffset.x += offsetX;
    cornerOffset.y += offsetY;

    if (uGlobal.uHasHeadRotation != 0u) {
        let rotatedOffset = (uGlobal.uInvHeadRotation * vec4<f32>(cornerOffset, 0.0, 0.0)).xyz;
        clip += uGlobal.uProjection * vec4<f32>(rotatedOffset, 0.0);
    } else {
        // apply offset in clip space to avoid perspective distortion on the quad
        clip = vec4<f32>(clip.xy + vec2<f32>(uGlobal.uProjection[0][0], uGlobal.uProjection[1][1]) * cornerOffset, clip.zw);
    }

    out.position = clip;
    out.vViewPosition = -mvCenter.xyz;
    ${clip_instance_wgsl}
    return out;
}

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;
    ${clip_flag_wgsl('in.vMarkerClipping.y')}

    // anti-aliased SDF text (.v6/designs/labels.md §5)
    // sampled and differentiated before anything can discard: derivatives need uniform control flow
    let rawSdf = sample_tFont(in.vTexCoord).r;
    // how much the distance field changes across one pixel, at any label size
    let sdfPerPixel = max(length(vec2<f32>(dpdx(rawSdf), dpdy(rawSdf))), 1e-4);

    if (fadeLod(in.position, in.vModelPosition)) {
        discard;
    }
    ${clip_pixel_wgsl('in.vModelPosition')}

    let fragmentDepth = in.position.z;
    let interior = false;

    // determine if this is a background or glyph fragment
    let isBackground = in.vTexCoord.x > 1.0;

    // discard background for pick/marking/emissive, but keep it in color, tracing and depth
    #if !defined(dRenderVariant_color) && !defined(dRenderVariant_tracing) && !defined(dRenderVariant_depth)
        if (isBackground) {
            discard;
        }
    #endif

    // glyph + border coverage: in color and tracing, anti-aliased over one pixel (the fractional edge is drawn in the
    // transparent pass, check_transparency); in the other variants, the edge at the field's middle as before
    var edgeAlpha = 1.0;
    var fillAlpha = 1.0;
    if (!isBackground) {
        // clamp to avoid exceeding max SDF range
        let sdf = rawSdf + min(uMaterial.uBorderWidth, 0.49);
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            // edges a quarter pixel out: keeps one-pixel stems of small text solid, invisible on large text
            edgeAlpha = clamp((sdf - 0.5) / sdfPerPixel + 0.75, 0.0, 1.0);
            fillAlpha = clamp((rawSdf - 0.5) / sdfPerPixel + 0.75, 0.0, 1.0);
            if (edgeAlpha < 1.0 / 255.0) {
                discard;
            }
        #else
            if (sdf < 0.5) {
                discard;
            }
        #endif
    }

    out.fragDepth = fragmentDepth;

    ${assign_material_wgsl}

    if (isBackground) {
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            material = vec4<f32>(uMaterial.uBackgroundColor, uMaterial.uBackgroundOpacity * material.a);
        #endif
    } else {
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            material = vec4<f32>(material.rgb, material.a * edgeAlpha);
            if (uMaterial.uBorderWidth > 0.0) {
                material = vec4<f32>(mix(uMaterial.uBorderColor, material.rgb, fillAlpha), material.a);
            }
            if (uMaterial.uBorderWidth == 0.0 || fillAlpha >= 0.5) {
                // push text fragments forward in depth so they render in front of border (depth is reversed)
                out.fragDepth = fragmentDepth + 0.0001;
            }
        #endif
    }

    ${apply_occlusion_wgsl}
    ${check_transparency_wgsl}

    #if defined(dRenderVariant_pick)
        ${pick_output_wgsl}
    #elif defined(dRenderVariant_marking)
        ${marking_output_wgsl}
    #elif defined(dRenderVariant_emissive)
        ${emissive_output_wgsl}
    #elif defined(dRenderVariant_depth)
        ${depth_output_wgsl('select(1.0, uMaterial.uBackgroundOpacity, isBackground)')}
    #else
        var color = material;
        ${apply_marker_wgsl}
        #if defined(dRenderVariant_tracing)
            ${tracing_output_wgsl('-normalize(in.vViewPosition)')}
        #else
            ${apply_fog_wgsl}
            ${wboit_output_wgsl}
            ${dpoit_output_wgsl}
            out.color = color;
        #endif
    #endif
    return out;
}
`;
