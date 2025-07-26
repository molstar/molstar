/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const text_vert = `
precision highp float;
precision highp int;

#include common
#include read_from_texture
#include common_vert_params
#include color_vert_params
#include size_vert_params
#include common_clip

uniform mat4 uModelView;

attribute vec3 aPosition;
attribute vec2 aMapping;
attribute float aDepth;
attribute vec2 aTexCoord;
attribute mat4 aTransform;
attribute float aInstance;
attribute float aGroup;

uniform float uOffsetX;
uniform float uOffsetY;
uniform float uOffsetZ;

uniform float uIsOrtho;
uniform float uPixelRatio;
uniform vec4 uViewport;
uniform mat4 uInvHeadRotation;
uniform bool uHasHeadRotation;

varying vec2 vTexCoord;

void main(void){
    int vertexId = VertexID;

    #include assign_group
    #include assign_color_varying
    #include assign_marker_varying
    #include assign_clipping_varying
    #include assign_size

    vTexCoord = aTexCoord;

    float scale = uModelScale;

    float offsetX = uOffsetX * scale;
    float offsetY = uOffsetY * scale;
    float offsetZ = (uOffsetZ + aDepth * 0.95) * scale;

    vec4 position4 = vec4(aPosition, 1.0);
    vec4 mvPosition = uModelView * aTransform * position4;

    vModelPosition = (uModel * aTransform * position4).xyz; // for clipping in frag shader

    // TODO
    // #ifdef FIXED_SIZE
    //     if (ortho) {
    //         scale /= pixelRatio * ((uViewport.w / 2.0) / -uCameraPosition.z) * 0.1;
    //     } else {
    //         scale /= pixelRatio * ((uViewport.w / 2.0) / -mvPosition.z) * 0.1;
    //     }
    // #endif

    vec4 mvCorner = vec4(mvPosition.xyz, 1.0);

    if (vTexCoord.x == 10.0) { // indicates background plane
        // move a bit to the back, taking distance to camera into account to avoid z-fighting
        offsetZ -= 0.001 * distance(uCameraPosition, (uProjection * mvCorner).xyz);
    }

    vec3 cornerOffset = vec3(0.0);
    cornerOffset.xy += aMapping * size * scale;
    cornerOffset.x += offsetX;
    cornerOffset.y += offsetY;

    if (uHasHeadRotation) {
        mvCorner.xyz += (uInvHeadRotation * vec4(cornerOffset, 1.0)).xyz;
    } else {
        mvCorner.xyz += cornerOffset;
    }

    if (uIsOrtho == 1.0) {
        mvCorner.z += offsetZ;
    } else {
        mvCorner.xyz += normalize(-mvCorner.xyz) * offsetZ;
    }

    gl_Position = uProjection * mvCorner;

    vViewPosition = -mvCorner.xyz;

    #include clip_instance
}
`;