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
    let transform = uObject.uObjectTransform * applyTumble(aTransform, invariantBoundingSphere(in));
    let model = uGlobal.uModel * transform;
    let modelView = uGlobal.uView * model;
    ${assign_wiggle_wgsl}
    position = applyWiggle(position, group, wiggle);
    let position4 = vec4<f32>(position, 1.0);
    out.vModelPosition = (model * position4).xyz;
    var mvPosition = modelView * position4;
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

    #if defined(dPointFront)
        // twice its radius toward the camera along the view ray (an orthographic projection's [2][3] is 0): a point's
        // radius in view units is its pixels over the projection's y scale and the viewport's height, times the depth
        let projection = uGlobal.uProjection;
        let ortho = projection[2][3] == 0.0;
        let radius = pointSize / (projection[1][1] * viewport.w) * select(-mvPosition.z, 1.0, ortho);
        if (ortho) {
            mvPosition.z += 2.0 * radius;
        } else {
            mvPosition = vec4<f32>(mvPosition.xyz * max(0.0, 1.0 - 2.0 * radius / length(mvPosition.xyz)), mvPosition.w);
        }
    #endif

    let p = uGlobal.uProjection * mvPosition;
    // offset the corner by half the point size in pixels, in clip space
    out.position = vec4<f32>(p.xy + mapping * (pointSize / viewport.zw) * p.w, p.zw);
    out.vPointCoord = mapping * 0.5 + 0.5;
    if (uObject.uPointDepthOffset != 0.0) {
        // in depth only (its pixels stay, as a mesh's \`uDepthOffset\`): that many CSS pixels at the point's depth, kept
        // in front of the camera and within the depth range (the planes fit the bounds: past them, at a plane, still in
        // front of or behind what is drawn)
        let projection = uGlobal.uProjection;
        let ortho = projection[2][3] == 0.0;
        let cssPixel = select(-2.0 * mvPosition.z, 2.0, ortho) / (projection[1][1] * uGlobal.uCssHeight);
        var offsetZ = mvPosition.z + uObject.uPointDepthOffset * cssPixel;
        if (!ortho) {
            offsetZ = min(offsetZ, mvPosition.z * 1e-3);
        }
        let offsetPosition = projection * vec4<f32>(mvPosition.xy, offsetZ, 1.0);
        out.position.z = clamp(offsetPosition.z / offsetPosition.w, 1e-6, 1.0 - 1e-6) * out.position.w;
    }

    if (lodCulled) {
        // move out of [0, w] to 'discard' in the vertex stage
        out.position.z = 2.0 * out.position.w;
    }
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping: every fragment of the point is at its center, the point tested once; the ghost pass
        // (\`uClipGhost\`) keeps what the clip objects remove
        if (clipTest(out.vModelPosition / uGlobal.uModelScale, clipFlag) != (uGlobal.uClipGhost != 0u)) {
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
    #elif defined(dPointStyle_ring)
        // a ring: the circle less an inner one, its width a fraction of the radius
        let ringDist = distance(in.vPointCoord, vec2<f32>(0.5));
        if (ringDist > 0.5 || ringDist < 0.5 * (1.0 - uObject.uPointRingWidth)) {
            discard;
        }
    #elif defined(dPointStyle_fuzzy)
        let fuzzyAlpha = 1.0 - smoothstep(0.0, 0.5, distance(in.vPointCoord, vec2<f32>(0.5)));
        if (fuzzyAlpha < 0.0001) {
            discard;
        }
        material.a *= fuzzyAlpha;
    #endif

    ${apply_occlusion_wgsl}
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
