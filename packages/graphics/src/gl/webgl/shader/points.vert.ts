/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';
import { assign_quad_position_group } from './chunks/assign-quad-position-group.glsl.js';
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
import { assign_lod_factor } from './chunks/assign-lod-factor.glsl.js';

/** Points are drawn as screen-aligned quads of two triangles (6 vertices) with the point data read per quad */
export const points_vert = `
precision highp float;
precision highp int;

${common}
${read_from_texture}
${common_vert_params}
${color_vert_params}
${common_clip}
${common_animation}

/** position of the fragment within the point, from (0, 0) to (1, 1) */
varying vec2 vPointCoord;

void main(){
    ${assign_quad_position_group}

    ${assign_color_varying}
    ${assign_marker_varying}
    ${assign_clipping_varying}

    mat4 transform = uObjectTransform * applyTumble(aTransform);
    mat4 model = uModel * transform;
    mat4 modelView = uView * model;
    position = applyWiggle(position, group, aInstance);
    vec4 position4 = vec4(position, 1.0);
    // for clipping in frag shader
    vModelPosition = (model * position4).xyz;
    vec4 mvPosition = modelView * position4;
    vViewPosition = mvPosition.xyz;

    ${assign_size}

    #ifdef dPointSizeAttenuation
        float pointSize = size * uPixelRatio * ((uViewport.w / 2.0) / -mvPosition.z) * 5.0;
    #else
        float pointSize = size * uPixelRatio;
    #endif

    ${assign_lod_factor}
    pointSize = max(1.0, pointSize * lodFactor);

    #ifdef dPointFront
        // twice its radius toward the camera along the view ray (an orthographic projection's [2][3] is 0): a point's
        // radius in view units is its pixels over the projection's y scale and the viewport's height, times the depth
        bool ortho = uProjection[2][3] == 0.0;
        float radius = pointSize / (uProjection[1][1] * uViewport.w) * (ortho ? 1.0 : -mvPosition.z);
        if (ortho) {
            mvPosition.z += 2.0 * radius;
        } else {
            mvPosition.xyz *= max(0.0, 1.0 - 2.0 * radius / length(mvPosition.xyz));
        }
    #endif

    gl_Position = uProjection * mvPosition;
    // offset the corner by half the point size in pixels, in clip space
    gl_Position.xy += mapping * (pointSize / uViewport.zw) * gl_Position.w;
    vPointCoord = mapping * 0.5 + 0.5;
    if (uPointDepthOffset != 0.0) {
        // in depth only (its pixels stay, as a mesh's \`uDepthOffset\`): that many CSS pixels at the point's depth, kept
        // in front of the camera and within the depth range (the planes fit the bounds: past them, at a plane, still in
        // front of or behind what is drawn)
        bool orthoOffset = uProjection[2][3] == 0.0;
        float cssPixel = (orthoOffset ? 2.0 : -2.0 * mvPosition.z) / (uProjection[1][1] * uCssHeight);
        float offsetZ = mvPosition.z + uPointDepthOffset * cssPixel;
        if (!orthoOffset) offsetZ = min(offsetZ, mvPosition.z * 1e-3);
        vec4 offsetPosition = uProjection * vec4(mvPosition.xy, offsetZ, 1.0);
        gl_Position.z = clamp(offsetPosition.z / offsetPosition.w, -1.0 + 1e-6, 1.0 - 1e-6) * gl_Position.w;
    }

    if (lodCulled) {
        // move out of [ -w, +w ] to 'discard' in vert shader
        gl_Position.z = 2.0 * gl_Position.w;
    }

    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping: every fragment of the point is at its center, the point tested once; the ghost pass
        // (\`uClipGhost\`) keeps what the clip objects remove
        if (clipTest(vModelPosition / uModelScale) != uClipGhost) {
            gl_Position.z = 2.0 * gl_Position.w;
        }
    #endif

    ${clip_instance}
}
`;
