/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export default `
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

// uniform bool ortho;
uniform float uPixelRatio;
uniform float uViewportHeight;

varying vec2 vTexCoord;

#include matrix_scale

void main(void){
    #include assign_group
    #include assign_color_varying
    #include assign_marker_varying
    #include assign_clipping_varying
    #include assign_size

    vTexCoord = aTexCoord;

    float scale = matrixScale(uModelView);

    float offsetX = uOffsetX * scale;
    float offsetY = uOffsetY * scale;
    float offsetZ = (uOffsetZ + aDepth * 0.95) * scale;

    vec4 position4 = vec4(aPosition, 1.0);
    vec4 mvPosition = uModelView * aTransform * position4;

    vModelPosition = (uModel * aTransform * position4).xyz; // for clipping in frag shader

    // TODO
    // #ifdef FIXED_SIZE
    //     if (ortho) {
    //         scale /= pixelRatio * ((uViewportHeight / 2.0) / -uCameraPosition.z) * 0.1;
    //     } else {
    //         scale /= pixelRatio * ((uViewportHeight / 2.0) / -mvPosition.z) * 0.1;
    //     }
    // #endif

    vec4 mvCorner = vec4(mvPosition.xyz, 1.0);

    if (vTexCoord.x == 10.0) { // indicates background plane
        // move a bit to the back, taking distance to camera into account to avoid z-fighting
        offsetZ -= 0.001 * distance(uCameraPosition, (uProjection * mvCorner).xyz);
    }

    mvCorner.xy += aMapping * size * scale;
    mvCorner.x += offsetX;
    mvCorner.y += offsetY;

    // TODO
    // if(ortho){
    //     mvCorner.xyz += normalize(-uCameraPosition) * offsetZ;
    // } else {
    //     mvCorner.xyz += normalize(-mvCorner.xyz) * offsetZ;
    // }
    mvCorner.xyz += normalize(-mvCorner.xyz) * offsetZ;

    gl_Position = uProjection * mvCorner;

    vViewPosition = -mvCorner.xyz;

    #include clip_instance
}
`;