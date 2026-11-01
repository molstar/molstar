/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { animation_wgsl, assign_wiggle_wgsl } from './chunks/animation.wgsl.js';
import { clip_functions_wgsl, clip_instance_wgsl } from './chunks/clip.wgsl.js';
import { segment_wgsl } from './chunks/segment.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import {
  apply_fog_wgsl,
  apply_marker_wgsl,
  assign_material_wgsl,
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
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';
import {
  assign_quad_position_group_wgsl,
  assign_size_wgsl,
  assign_varyings_wgsl,
  vertex_output_wgsl,
  assign_lod_factor_wgsl,
} from './chunks/varyings.wgsl.js';

/** Points as screen-aligned quads with a size in pixels, unlit */
export const points_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${animation_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}

${vertex_output_wgsl('    @location(9) vPointCoord: vec2<f32>,')}
${fragment_output_wgsl()}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    ${assign_quad_position_group_wgsl}
    ${assign_varyings_wgsl}

    let aTransform = affineTransform(in.aTransform0, in.aTransform1, in.aTransform2);
    let transform = applyTumble(aTransform, invariantBoundingSphere(in));
    let model = uGlobal.uModel * transform;
    let modelView = uGlobal.uView * model;
    ${assign_wiggle_wgsl}
    position = applyWiggle(position, group, wiggle);
    let position4 = vec4<f32>(position, 1.0);
    out.vModelPosition = (model * position4).xyz;
    let mvPosition = modelView * position4;
    out.vViewPosition = mvPosition.xyz;

    ${assign_size_wgsl}

    let viewport = uGlobal.uViewport;
    #if defined(dPointSizeAttenuation)
        var pointSize = size * uGlobal.uPixelRatio * ((viewport.w / 2.0) / -mvPosition.z) * 5.0;
    #else
        var pointSize = size * uGlobal.uPixelRatio;
    #endif

    ${assign_lod_factor_wgsl}
    pointSize = max(1.0, pointSize * lodFactor);

    let p = uGlobal.uProjection * mvPosition;
    // offset the corner by half the point size in pixels, in clip space
    out.position = vec4<f32>(p.xy + mapping * (pointSize / viewport.zw) * p.w, p.zw);
    out.vPointCoord = mapping * 0.5 + 0.5;

    if (lodCulled) {
        // move out of [0, w] to 'discard' in the vertex stage
        out.position.z = 2.0 * out.position.w;
    }
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping: every fragment of the point is at its center, the point tested once
        if (clipTest(out.vModelPosition / uGlobal.uModelScale, clipFlag)) {
            out.position.z = 2.0 * out.position.w;
        }
    #endif
    ${clip_instance_wgsl}
    return out;
}

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;

    if (fadeLod(in.position, in.vModelPosition)) {
        discard;
    }
    // pixel clipping in the vertex stage: the point's fragments are all at its center

    let fragmentDepth = in.position.z;
    let interior = false;

    ${assign_material_wgsl}

    #if defined(dPointStyle_circle)
        if (distance(in.vPointCoord, vec2<f32>(0.5)) > 0.5) {
            discard;
        }
    #elif defined(dPointStyle_fuzzy)
        let fuzzyAlpha = 1.0 - smoothstep(0.0, 0.5, distance(in.vPointCoord, vec2<f32>(0.5)));
        if (fuzzyAlpha < 0.0001) {
            discard;
        }
        material.a *= fuzzyAlpha;
    #endif

    ${check_transparency_wgsl}

    #if defined(dRenderVariant_pick)
        ${pick_output_wgsl}
    #elif defined(dRenderVariant_marking)
        ${marking_output_wgsl}
    #elif defined(dRenderVariant_emissive)
        ${emissive_output_wgsl}
    #elif defined(dRenderVariant_depth)
        ${depth_output_wgsl()}
    #else
        var color = material;
        ${apply_marker_wgsl}
        #if defined(dRenderVariant_tracing)
            ${tracing_output_wgsl('normalize(in.vViewPosition)')}
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
