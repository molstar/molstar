/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * heavily based on code by WestLangley from https://github.com/WestLangley/three.js/blob/af28b2fb706ac109771ecad0a7447fad90ab3210/examples/js/lines/LineMaterial.js
 */

import { common } from './chunks/common.glsl.js';
import { read_from_texture } from './chunks/read-from-texture.glsl.js';
import { common_vert_params } from './chunks/common-vert-params.glsl.js';
import { color_vert_params } from './chunks/color-vert-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { common_animation } from './chunks/common-animation.glsl.js';
import { assign_color_varying } from './chunks/assign-color-varying.glsl.js';
import { assign_marker_varying } from './chunks/assign-marker-varying.glsl.js';
import { assign_clipping_varying } from './chunks/assign-clipping-varying.glsl.js';
import { assign_size } from './chunks/assign-size.glsl.js';
import { clip_instance } from './chunks/clip-instance.glsl.js';

export const lines_vert = `
precision highp float;
precision highp int;

${common}
${read_from_texture}
${common_vert_params}
${color_vert_params}
${common_clip}
${common_animation}

#if defined(dClipVariant_pixel) && defined(dClipObjects)
    // 1 where the clip objects may cut the line (\`clipSphereClass\` of its bounding sphere): its clip tests skipped if 0
    #if __VERSION__ == 100
        varying float vClipSome;
    #else
        flat out float vClipSome;
    #endif
#endif

void trimSegment(const in vec4 start, inout vec4 end) {
    // trim end segment so it terminates between the camera plane and the near plane
    // conservative estimate of the near plane, from the z row of the GL projection (the sign of the reversed rows
    // cancels in the ratio)
    #ifdef clipZeroToOne
        float a = 2.0 * uProjection[2][2] + 1.0;
        float b = 2.0 * uProjection[3][2];
    #else
        float a = uProjection[2][2];  // 3rd entry in 3rd column
        float b = uProjection[3][2];  // 3rd entry in 4th column
    #endif
    float nearEstimate = -0.5 * b / a;
    float alpha = (nearEstimate - start.z) / (end.z - start.z);
    end.xyz = mix(start.xyz, end.xyz, alpha);
}

void main(){
    float aspect = uViewport.z / uViewport.w;
    #if defined(dInstanceData)
        readInstanceData();
    #endif

    // the line of the vertex (2 texels of \`tLine\`) and its corner (4 vertices per line: left and right of the start,
    // then of the end)
    int item = VertexID / 4;
    int corner = VertexID - item * 4;
    // the data vertex ('vertex' color, size): the start for the corners at the start, else the end
    int vertexId = item * 2 + (corner < 2 ? 0 : 1);
    vec4 lineStart, lineEnd; // start xyz and group, end xyz and whether the line continues the one before
    read_tLine2(item * 2, lineStart, lineEnd);
    float group = lineStart.w;
    // -1 or 1: left or right, start or end
    vec2 mapping = vec2(corner == 0 || corner == 2 ? -1.0 : 1.0, corner < 2 ? -1.0 : 1.0);
    vec3 segmentStart = lineStart.xyz;
    vec3 segmentEnd = lineEnd.xyz;
    // the start corners of a line continuing the one before (strips): that line's end corners, across it
    if (lineEnd.w > 0.5 && mapping.y < 0.0) {
        segmentEnd = lineStart.xyz;
        segmentStart = read_tLine(item * 2 - 2).xyz;
        mapping.y = 1.0;
    }

    ${assign_color_varying}
    ${assign_marker_varying}
    ${assign_clipping_varying}
    ${assign_size}

    mat4 transform = applyTumble(aTransform);
    vec3 wigStart = applyWiggle(segmentStart, group, aInstance);
    vec3 wigEnd = applyWiggle(segmentEnd, group, aInstance);
    mat4 modelView = uView * uModel * transform;

    // camera space
    vec4 start = modelView * vec4(wigStart, 1.0);
    vec4 end = modelView * vec4(wigEnd, 1.0);

    // assign position
    vec4 position4 = vec4((mapping.y < 0.5) ? wigStart : wigEnd, 1.0);
    vViewPosition = (mapping.y < 0.5) ? start.xyz : end.xyz;

    vModelPosition = (uModel * transform * position4).xyz; // for clipping in frag shader

    // special case for perspective projection, and segments that terminate either in, or behind, the camera plane
    // clearly the gpu firmware has a way of addressing this issue when projecting into ndc space
    // but we need to perform ndc-space calculations in the shader, so we must address this issue directly
    // perhaps there is a more elegant solution -- WestLangley
    bool perspective = (uProjection[2][3] == -1.0); // 4th entry in the 3rd column
    if (perspective) {
        if (start.z < 0.0 && end.z >= 0.0) {
            trimSegment(start, end);
        } else if (end.z < 0.0 && start.z >= 0.0) {
            trimSegment(end, start);
        }
    }

    // clip space
    vec4 clipStart = uProjection * start;
    vec4 clipEnd = uProjection * end;

    // ndc space
    vec2 ndcStart = clipStart.xy / clipStart.w;
    vec2 ndcEnd = clipEnd.xy / clipEnd.w;

    // direction
    vec2 dir = ndcEnd - ndcStart;

    // account for clip-space aspect ratio
    dir.x *= aspect;
    dir = normalize(dir);

    // perpendicular to dir
    vec2 offset = vec2(dir.y, - dir.x);

    // undo aspect ratio adjustment
    dir.x /= aspect;
    offset.x /= aspect;

    // sign flip
    if (mapping.x < 0.0) offset *= -1.0;

    // calculate linewidth
    float linewidth;
    #ifdef dLineSizeAttenuation
        linewidth = size * uPixelRatio * ((uViewport.w / 2.0) / -start.z) * 5.0;
    #else
        linewidth = size * uPixelRatio;
    #endif
    linewidth = max(1.0, linewidth);

    // adjust for linewidth
    offset *= linewidth;

    // adjust for clip-space to screen-space conversion
    offset /= uViewport.w;

    // select end
    vec4 clip = (mapping.y < 0.5) ? clipStart : clipEnd;

    // back to clip space
    offset *= clip.w;
    clip.xy += offset;
    gl_Position = clip;

    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping by the bounding sphere of the item's line (not the segment of a strip's start corners, which
        // continue the line before): a line clipped whole moved out of [ -w, +w ], one not cut not tested per fragment
        mat4 modelTransform = uModel * transform;
        vec3 lineA = (modelTransform * vec4(applyWiggle(lineStart.xyz, group, aInstance), 1.0)).xyz;
        vec3 lineB = (modelTransform * vec4(applyWiggle(lineEnd.xyz, group, aInstance), 1.0)).xyz;
        int clipClass = clipSphereClass((lineA + lineB) * 0.5 / uModelScale, length(lineB - lineA) * 0.5 / uModelScale);
        if (clipClass == 2) {
            gl_Position.z = 2.0 * gl_Position.w;
        }
        vClipSome = clipClass == 1 ? 1.0 : 0.0;
    #endif

    ${clip_instance}
}
`;
