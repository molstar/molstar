/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const mesh_vert = `
precision highp float;
precision highp int;
precision highp sampler2D;

#include common
#include read_from_texture
#include common_vert_params
#include color_vert_params
#include common_clip
#include texture3d_from_2d_linear

#ifdef dGeometryType_textureMesh
    uniform vec2 uGeoTexDim;
    uniform sampler2D tPosition;
    uniform sampler2D tGroup;
    uniform sampler2D tNormal;
#else
    attribute vec3 aPosition;
    attribute float aGroup;
    attribute vec3 aNormal;
#endif
attribute mat4 aTransform;
attribute float aInstance;

varying vec3 vNormal;

void main(){
    #include assign_group
    #include assign_marker_varying
    #include assign_clipping_varying
    #include assign_position
    #include assign_color_varying
    #include clip_instance

    #ifdef dGeometryType_textureMesh
        vec3 normal = readFromTexture(tNormal, VertexID, uGeoTexDim).xyz;
    #else
        vec3 normal = aNormal;
    #endif
    mat3 normalMatrix = transpose3(inverse3(mat3(modelView)));
    vec3 transformedNormal = normalize(normalMatrix * normalize(normal));
    #if defined(dFlipSided)
        if (!uDoubleSided) { // TODO checking uDoubleSided should not be required, ASR
            transformedNormal = -transformedNormal;
        }
    #endif
    vNormal = transformedNormal;
}
`;