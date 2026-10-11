/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { common } from './chunks/common.glsl.js';
import { read_from_texture } from './chunks/read-from-texture.glsl.js';
import { common_vert_params } from './chunks/common-vert-params.glsl.js';
import { color_vert_params } from './chunks/color-vert-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { common_animation } from './chunks/common-animation.glsl.js';
import { texture3d_from_2d_linear } from './chunks/texture3d-from-2d-linear.glsl.js';
import { assign_group } from './chunks/assign-group.glsl.js';
import { assign_marker_varying } from './chunks/assign-marker-varying.glsl.js';
import { assign_clipping_varying } from './chunks/assign-clipping-varying.glsl.js';
import { assign_position } from './chunks/assign-position.glsl.js';
import { assign_color_varying } from './chunks/assign-color-varying.glsl.js';
import { clip_instance } from './chunks/clip-instance.glsl.js';

export const mesh_vert = `
precision highp float;
precision highp int;
precision highp sampler2D;

${common}
${read_from_texture}
${common_vert_params}
${color_vert_params}
${common_clip}
${common_animation}
${texture3d_from_2d_linear}

varying vec3 vNormal;

void main(){
    int vertexId = VertexID;

    ${assign_group}
    ${assign_marker_varying}
    ${assign_clipping_varying}
    ${assign_position}

    #ifdef dSolidInterior
        if (uSolidInteriorPass != 0) {
            gl_Position.z = min(gl_Position.z, gl_Position.w * 0.9999);
        }
    #endif

    ${assign_color_varying}
    ${clip_instance}

    #ifdef dGeometryType_textureMesh
        vec3 normal = readFromTexture(tNormal, vertexId, uGeoTexDim).xyz;
    #else
        vec3 normal = aNormal;
    #endif
    mat3 normalMatrix = adjoint(modelView);
    vec3 transformedNormal = normalize(normalMatrix * normalize(normal));
    #if defined(dFlipSided)
        if (!uDoubleSided) { // TODO checking uDoubleSided should not be required, ASR
            transformedNormal = -transformedNormal;
        }
    #endif
    vNormal = transformedNormal;
}
`;
