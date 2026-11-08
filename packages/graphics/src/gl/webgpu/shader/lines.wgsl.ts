/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { animation_wgsl, assign_wiggle_wgsl } from './chunks/animation.wgsl.js';
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
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';
import { assign_size_wgsl, assign_varyings_wgsl, vertex_output_wgsl } from './chunks/varyings.wgsl.js';

/** Line segments as screen-space quads (4 vertices, 2 triangles per segment), unlit */
export const lines_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${animation_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}

${vertex_output_wgsl(`#if defined(dClipVariant_pixel) && defined(dClipObjects)
    // 1 where the clip objects may cut the line (\`clipSphereClass\` of its bounding sphere): its clip tests skipped
    // if 0
    @location(9) @interpolate(flat) vClipSome: u32,
#endif`)}
${fragment_output_wgsl()}

/**
 * Trims the end of a segment so it terminates between the camera plane and the near plane, with the conservative near
 * plane estimate of the GL projection (z row of the [-1, 1] projection: -(2 * reversed [0, 1] row + w row), the sign
 * cancels in the ratio)
 */
fn trimSegment(start: vec4<f32>, end: vec4<f32>) -> vec4<f32> {
    let a = 2.0 * uGlobal.uProjection[2][2] + 1.0;
    let b = 2.0 * uGlobal.uProjection[3][2];
    let nearEstimate = -0.5 * b / a;
    let alpha = (nearEstimate - start.z) / (end.z - start.z);
    return vec4<f32>(mix(start.xyz, end.xyz, alpha), end.w);
}


@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let viewport = uGlobal.uViewport;
    let aspect = viewport.z / viewport.w;

    // the line of the vertex (2 texels of \`tLine\`) and its corner (4 vertices per line: left and right of the start,
    // then of the end)
    let item = i32(in.vertexIndex) / 4;
    let corner = i32(in.vertexIndex) - item * 4;
    // the data vertex ('vertex' color, size): the start for the corners at the start, else the end
    let vertexId = item * 2 + select(1, 0, corner < 2);
    let lineStart = read_tLine(item * 2); // start xyz, group
    let lineEnd = read_tLine(item * 2 + 1); // end xyz, whether the line continues the one before
    let group = lineStart.w;
    // -1 or 1: left or right, start or end
    var mapping = vec2<f32>(select(1.0, -1.0, corner == 0 || corner == 2), select(1.0, -1.0, corner < 2));
    var segmentStart = lineStart.xyz;
    var segmentEnd = lineEnd.xyz;
    // the start corners of a line continuing the one before (strips): that line's end corners, across it
    if (lineEnd.w > 0.5 && mapping.y < 0.0) {
        segmentEnd = lineStart.xyz;
        segmentStart = read_tLine(item * 2 - 2).xyz;
        mapping.y = 1.0;
    }
    ${assign_varyings_wgsl}
    ${assign_size_wgsl}

    let aTransform = affineTransform(in.aTransform0, in.aTransform1, in.aTransform2);
    let transform = uObject.uObjectTransform * applyTumble(aTransform, invariantBoundingSphere(in));
    ${assign_wiggle_wgsl}
    let wigStart = applyWiggle(segmentStart, group, wiggle);
    let wigEnd = applyWiggle(segmentEnd, group, wiggle);
    let model = uGlobal.uModel * transform;
    let modelView = uGlobal.uView * model;

    // camera space
    var start = modelView * vec4<f32>(wigStart, 1.0);
    var end = modelView * vec4<f32>(wigEnd, 1.0);

    // assign position
    let isStart = mapping.y < 0.5;
    let position4 = vec4<f32>(select(wigEnd, wigStart, isStart), 1.0);
    out.vViewPosition = select(end.xyz, start.xyz, isStart);
    out.vModelPosition = (model * position4).xyz;

    // special case for perspective projection, and segments that terminate either in, or behind, the camera plane
    let perspective = uGlobal.uProjection[2][3] == -1.0;
    if (perspective) {
        if (start.z < 0.0 && end.z >= 0.0) {
            end = trimSegment(start, end);
        } else if (end.z < 0.0 && start.z >= 0.0) {
            start = trimSegment(end, start);
        }
    }

    // clip space
    let clipStart = uGlobal.uProjection * start;
    let clipEnd = uGlobal.uProjection * end;

    // ndc space
    let ndcStart = clipStart.xy / clipStart.w;
    let ndcEnd = clipEnd.xy / clipEnd.w;

    // direction, accounting for the clip-space aspect ratio
    var dir = ndcEnd - ndcStart;
    dir.x *= aspect;
    dir = normalize(dir);

    // perpendicular to dir, undoing the aspect ratio adjustment
    var offset = vec2<f32>(dir.y, -dir.x);
    offset.x /= aspect;

    // sign flip
    if (mapping.x < 0.0) {
        offset *= -1.0;
    }

    // calculate linewidth
    #if defined(dLineSizeAttenuation)
        var linewidth = size * uGlobal.uPixelRatio * ((viewport.w / 2.0) / -start.z) * 5.0;
    #else
        var linewidth = size * uGlobal.uPixelRatio;
    #endif
    linewidth = max(1.0, linewidth);

    // adjust for linewidth and the clip-space to screen-space conversion
    offset *= linewidth;
    offset /= viewport.w;

    // select end, back to clip space
    var clip = select(clipEnd, clipStart, isStart);
    offset *= clip.w;
    clip = vec4<f32>(clip.xy + offset, clip.zw);
    out.position = clip;
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping by the bounding sphere of the item's line (not the segment of a strip's start corners, which
        // continue the line before): a line clipped whole moved out of [0, w], one not cut not tested per fragment
        let modelTransform = uGlobal.uModel * transform;
        let lineA = (modelTransform * vec4<f32>(applyWiggle(lineStart.xyz, group, wiggle), 1.0)).xyz;
        let lineB = (modelTransform * vec4<f32>(applyWiggle(lineEnd.xyz, group, wiggle), 1.0)).xyz;
        // the ghost pass (\`uClipGhost\`): a line whole by its center, kept where clipped, not tested per fragment
        var clipClass: u32;
        if (uGlobal.uClipGhost != 0u) {
            clipClass = select(2u, 0u, clipTest((lineA + lineB) * 0.5 / uGlobal.uModelScale, clipFlag));
        } else {
            clipClass = clipSphereClass((lineA + lineB) * 0.5 / uGlobal.uModelScale, length(lineB - lineA) * 0.5 / uGlobal.uModelScale, clipFlag);
        }
        if (clipClass == 2u) {
            out.position.z = 2.0 * out.position.w;
        }
        out.vClipSome = select(0u, 1u, clipClass == 1u);
    #endif
    ${clip_instance_wgsl}
    return out;
}

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;
    ${clip_flag_wgsl('in.vMarkerClipping.y')}

    if (fadeLod(in.position, in.vModelPosition)) {
        discard;
    }
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        if (in.vClipSome != 0u) {
            ${clip_pixel_wgsl('in.vModelPosition')}
        }
    #endif

    let fragmentDepth = in.position.z;
    let interior = false;

    ${assign_material_wgsl}
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
