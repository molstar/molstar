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

    mat4 transform = applyTumble(aTransform, aInstance, float(uObjectId));
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
    pointSize = max(1.0, pointSize);

    gl_Position = uProjection * mvPosition;
    // offset the corner by half the point size in pixels, in clip space
    gl_Position.xy += mapping * (pointSize / uViewport.zw) * gl_Position.w;
    vPointCoord = mapping * 0.5 + 0.5;

    ${clip_instance}
}
`;
